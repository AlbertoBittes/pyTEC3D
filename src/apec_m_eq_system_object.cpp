#include "apec_m_eq_system_object.h"
#include "apec_m_system_objects.h"
#include "apec_m_mesh.h"
#include "time_step.h"
#include <libmesh/solver_configuration.h>
#include "solver_configuration.h"

#include <utility>

void ApecMEquationSystems::apply_ramped_effects(const Substep &time_substep)
{
    auto ramp_sets = collect_ramp_vectors();

    for (auto &[old_vec, current_vec, ramped_vec] : ramp_sets)
        ApecMSystem<System>::ramp_vectors(time_substep, *old_vec, *current_vec, *ramped_vec);
}

void ApecMEquationSystems::configure_systems()
{
    libMesh::out << "Building system of equations..." << std::endl;

    this->clear_preserving_solver_settings();

    libMesh::out << "Checking system options..." << std::endl;

    for (int i = 0; i < Globals::NUM_SYSTEMS; ++i)
    {
        Globals::SystemType type = static_cast<Globals::SystemType>(i);
        const std::string sys_name = Globals::get_info_by_key(Globals::SYSTEM_INFO_MAP, type).name;

        if (type == Globals::DOF_SOLUTION_SYSTEM)
        {
            auto &linear_system = this->add_system<ApecMLinearSystem>(sys_name);
            if (auto assembler = linear_system.get_assembler())
                assembler->set_parent_object(this);
        }
        else
            this->add_system<ApecMExplicitSystem>(sys_name);
    }
}

void ApecMEquationSystems::initialize_system_contexts()
{
    auto system_builder_view = this->get_context_view();
    if (!system_builder_view)
        return;

    for (int i = 0; i < Globals::NUM_SYSTEMS; ++i)
    {
        const auto type = static_cast<Globals::SystemType>(i);
        const auto &info = Globals::get_info_by_key(Globals::SYSTEM_INFO_MAP, type);
        auto &base_system = this->get_system<System>(info.name);

        if (type == Globals::DOF_SOLUTION_SYSTEM)
        {
            if (auto *system = dynamic_cast<ApecMLinearSystem *>(&base_system))
            {
                system->set_context_view(system_builder_view);
                system->init_variables();
            }
        }
        else if (auto *system = dynamic_cast<ApecMExplicitSystem *>(&base_system))
        {
            system->set_context_view(system_builder_view);
            system->init_variables();
        }
    }
}

void ApecMEquationSystems::post_initialize_systems()
{
    for (int i = 0; i < Globals::NUM_SYSTEMS; ++i)
    {
        const auto type = static_cast<Globals::SystemType>(i);
        const auto &info = Globals::get_info_by_key(Globals::SYSTEM_INFO_MAP, type);
        auto &base_system = this->get_system<System>(info.name);

        if (type == Globals::DOF_SOLUTION_SYSTEM)
        {
            if (auto *system = dynamic_cast<ApecMLinearSystem *>(&base_system))
                system->post_init();
        }
        else if (auto *system = dynamic_cast<ApecMExplicitSystem *>(&base_system))
        {
            system->post_init();
        }
    }
}

void ApecMEquationSystems::initialize_systems(const std::shared_ptr<FEMDataRepository> &model_data)
{
    if (!model_data)
    {
        libMesh::out << "ApecMEquationSystems: model repository not provided for physics initialization." << std::endl;
        return;
    }

    auto system_view = std::make_shared<ContextView>(ContextView::from(model_data));
    this->set_context_view(system_view);

    libMesh::out << "Distributing degrees of freedom..." << std::endl;
    initialize_system_contexts();
}

void ApecMEquationSystems::initialize_physics()
{
    this->init();
    post_initialize_systems();
    this->print_info();
}

void ApecMEquationSystems::clear_preserving_solver_settings()
{
    Real tol = parameters.get<Real>("linear solver tolerance");
    unsigned int max_iter = parameters.get<unsigned int>("linear solver maximum iterations");

    this->clear();

    parameters.set<Real>("linear solver tolerance") = tol;
    parameters.set<unsigned int>("linear solver maximum iterations") = max_iter;
}

void ApecMEquationSystems::zero_all_systems()
{
    for (int i = 0; i < Globals::NUM_SYSTEMS; ++i)
    {
        const auto type = static_cast<Globals::SystemType>(i);
        const auto &info = Globals::get_info_by_key(Globals::SYSTEM_INFO_MAP, type);
        auto &base_system = this->get_system<System>(info.name);

        if (auto *zeroable = dynamic_cast<ApecMSystemBase *>(&base_system))
            zeroable->zero_all();
    }
}

std::shared_ptr<ApecMEquationSystems> ApecMEquationSystems::build(ApecMMesh &mesh, APEC_M *program)
{
    auto equation_systems = std::make_shared<ApecMEquationSystems>(mesh);
    if (program)
        equation_systems->set_parent_object(program);
    equation_systems->configure_systems();
    return equation_systems;
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
