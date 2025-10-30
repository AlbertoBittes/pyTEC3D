#include "apec_m_eq_system_object.h"
#include "apec_m_system_objects.h"
#include "apec_m_mesh.h"
#include "element_builder.h"
#include "system_builder.h"
#include "time_step.h"

#include <algorithm>

void ApecMEquationSystems::apply_ramped_effects(const Substep &time_substep)
{
    auto ramp_sets = collect_ramp_vectors();

    for (auto &[old_vec, current_vec, ramped_vec] : ramp_sets)
        ApecMSystem<System>::ramp_vectors(time_substep, *old_vec, *current_vec, *ramped_vec);
}

void ApecMEquationSystems::configure_systems(const std::shared_ptr<FEMDataRepository> &model_data)
{
    libMesh::out << "Building system of equations..." << std::endl;

    this->clear_preserving_solver_settings();

    libMesh::out << "Checking system options..." << std::endl;

    auto system_builder_view = std::make_shared<
        TypedContextView<ElementTypesFactory>>(
        TypedContextView<ElementTypesFactory>::from(
            model_data));

    for (int i = 0; i < Globals::NUM_SYSTEMS; ++i)
    {
        Globals::SystemType type = static_cast<Globals::SystemType>(i);
        const std::string sys_name = Globals::get_info_by_key(Globals::SYSTEM_INFO_MAP, type).name;

        System *system = nullptr;
        if (type == Globals::DOF_SOLUTION_SYSTEM)
            system = &this->add_system<ApecMLinearSystem>(sys_name);
        else
            system = &this->add_system<ApecMExplicitSystem>(sys_name);

        SystemBuilder *builder = new SystemBuilder(system, system_builder_view);
        builder->init_variables();
        system->attach_init_object(*builder);
    }
}

void ApecMEquationSystems::initialize_physics()
{
    libMesh::out << "Distributing degrees of freedom..." << std::endl;
    this->init();
    this->print_info();
    this->get_element_builder()->init();
}

void ApecMEquationSystems::clear_preserving_solver_settings()
{
    Real tol = parameters.get<Real>("linear solver tolerance");
    unsigned int max_iter = parameters.get<unsigned int>("linear solver maximum iterations");

    this->clear();

    parameters.set<Real>("linear solver tolerance") = tol;
    parameters.set<unsigned int>("linear solver maximum iterations") = max_iter;
}

std::shared_ptr<ApecMEquationSystems> ApecMEquationSystems::build(ApecMMesh &mesh,
                                                                  const std::shared_ptr<FEMDataRepository> &model_data)
{
    auto equation_systems = std::make_shared<ApecMEquationSystems>(mesh);
    equation_systems->configure_systems(model_data);
    return equation_systems;
}

ApecMEquationSystems::AnalysisActivationState ApecMEquationSystems::activate_analysis(const std::string &analysis_name)
{
    if (analysis_name.empty())
        return AnalysisActivationState::AlreadyActive;

    const auto previous_active = active_analysis_name;

    auto it = std::find_if(analyses.begin(), analyses.end(),
                           [&analysis_name](const ApecMAnalysis &analysis)
                           { return analysis.get_name() == analysis_name; });

    AnalysisActivationState result = AnalysisActivationState::AlreadyActive;

    if (it == analyses.end())
    {
        analyses.emplace_back(analysis_name);
        result = AnalysisActivationState::Created;
    }
    else if (!previous_active || *previous_active != analysis_name)
    {
        result = AnalysisActivationState::ActivatedExisting;
    }

    active_analysis_name = analysis_name;

    if (result == AnalysisActivationState::AlreadyActive &&
        (!previous_active || *previous_active != analysis_name))
    {
        result = AnalysisActivationState::ActivatedExisting;
    }

    return result;
}

void ApecMEquationSystems::clear_active_analysis()
{
    active_analysis_name.reset();
}

bool ApecMEquationSystems::has_analysis(const std::string &analysis_name) const
{
    return std::any_of(analyses.begin(), analyses.end(),
                       [&analysis_name](const ApecMAnalysis &analysis)
                       { return analysis.get_name() == analysis_name; });
}

std::vector<ApecMEquationSystems::RampVectorTuple> ApecMEquationSystems::collect_ramp_vectors()
{
    std::vector<RampVectorTuple> ramp_sets;

    System &main_system = this->get_system<System>("ApecMLinearSystem");
    ramp_sets.emplace_back(
        main_system.request_vector("old_node_dof_constraints"),
        main_system.request_vector("node_dof_constraints"),
        main_system.request_vector("ramped_node_dof_constraints"));

    System &element_loads_system = this->get_system<System>("ElementLoadsSystem");
    ramp_sets.emplace_back(
        element_loads_system.request_vector("old_loads"),
        element_loads_system.current_local_solution.get(),
        element_loads_system.request_vector("ramped_loads"));

    System &nodal_loads_system = this->get_system<System>("NodalLoadsSystem");
    ramp_sets.emplace_back(
        nodal_loads_system.request_vector("old_loads"),
        nodal_loads_system.current_local_solution.get(),
        nodal_loads_system.request_vector("ramped_loads"));

    return ramp_sets;
}
