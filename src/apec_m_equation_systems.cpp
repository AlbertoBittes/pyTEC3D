#include "apec_m_equation_systems.h"
#include "apec_m_system.h"
#include "apec_m_mesh.h"
#include "element_types_factory.h"
#include "time_step.h"
#include "vtkio_updated.h"
#include "session.h"
#include <libmesh/solver_configuration.h>
#include "solver_configuration.h"

#include <optional>
#include <set>
#include <utility>

void ApecMEquationSystems::apply_ramped_effects(const Substep &time_substep)
{
    auto ramp_sets = collect_ramp_vectors();

    for (auto &[system, old_vec, current_vec, ramped_vec] : ramp_sets)
        system->ramp_vectors(time_substep, *old_vec, *current_vec, *ramped_vec);
}

void ApecMEquationSystems::configure_services(const std::shared_ptr<AnalysisContextView> &analysis_ctx)
{
    analysis_context_view_ = analysis_ctx;

    auto &system = this->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    auto bc_view = analysis_ctx ? analysis_ctx->template project_view<TimeStep>() : nullptr;
    auto global_constraints_object = system.get_global_constraints_object();
    global_constraints_object->set_context_view(bc_view);

    auto assembler_view = analysis_ctx ? analysis_ctx->template project_view<TimeStep, Session>() : nullptr;
    auto assembler = system.get_assembler();
    assembler->set_context_view(assembler_view);
}

void ApecMEquationSystems::prepare_time_substep_solve(Substep &substep)
{
    auto &system = this->get_system<ApecMLinearSystem>("ApecMLinearSystem");

    this->apply_ramped_effects(substep);

    system.get_dof_map().swap_dof_constraints();

    system.reinit_constraints();
}

void ApecMEquationSystems::solve_equation_systems(Substep &substep, bool nodal_solution)
{
    auto &system = this->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    auto assembler = system.get_assembler();

    assembler->write_element_matrices_and_rhs();

    system.solve();
    system.get_dof_map().enforce_constraints_exactly(system);
    system.compute_solution_time_derivative(substep);
    assembler->write_solution_vectors(nodal_solution);
}

void ApecMEquationSystems::print_output(Substep &substep, bool picard_iteration_print)
{
    libMesh::out << "Writing output for substep " << substep.get_index() << std::endl;

    auto time_step = analysis_context_view_ ? analysis_context_view_->template get<TimeStep>() : nullptr;
    auto session = analysis_context_view_ ? analysis_context_view_->template get<Session>() : nullptr;

    if (!time_step || !session)
        return;

    const bool should_write = substep.get_write_results() ||
                              substep.get_index() == substep.get_parent_object()->get_n_substeps() - 1;

    if (!should_write)
        return;

    const auto base_substep_index = time_step->count_substeps_before_current();
    std::set<std::string> system_names = {"ApecMLinearSystem", "NodalSolutionSystem"};

    std::optional<unsigned int> picard_iteration;
    if (picard_iteration_print)
        picard_iteration = substep.get_number_picard_iterations();
    else
        picard_iteration = std::nullopt;

    session->set_output_filepath(base_substep_index + substep.get_index(), picard_iteration);
    VTKIO_Updated vtk_output(this->get_mesh(),
                             this);
    vtk_output.write_equation_systems(session->get_path("final_output_path"),
                                      *this,
                                      &system_names);
}

void ApecMEquationSystems::update_old_local_solution()
{
    auto &system = this->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    if (!system.old_local_solution || !system.current_local_solution)
        return;
    *system.old_local_solution = *system.current_local_solution;
}

void ApecMEquationSystems::configure_systems()
{
    libMesh::out << "Building system of equations..." << std::endl;

    this->clear_preserving_solver_settings();

    libMesh::out << "Checking system options..." << std::endl;

    for (int i = 0; i < ApecMSystemBase::NUM_SYSTEMS; ++i)
    {
        ApecMSystemBase::SystemType type = static_cast<ApecMSystemBase::SystemType>(i);
        const std::string sys_name = ApecMSystemBase::SYSTEM_INFO_MAP.at(type);

        if (type == ApecMSystemBase::DOF_SOLUTION_SYSTEM)
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

    for (int i = 0; i < ApecMSystemBase::NUM_SYSTEMS; ++i)
    {
        const auto type = static_cast<ApecMSystemBase::SystemType>(i);
        const auto &info = ApecMSystemBase::SYSTEM_INFO_MAP.at(type);
        auto &base_system = this->get_system<System>(info);

        if (type == ApecMSystemBase::DOF_SOLUTION_SYSTEM)
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
    for (int i = 0; i < ApecMSystemBase::NUM_SYSTEMS; ++i)
    {
        const auto type = static_cast<ApecMSystemBase::SystemType>(i);
        const auto &info = ApecMSystemBase::SYSTEM_INFO_MAP.at(type);
        auto &base_system = this->get_system<System>(info);

        if (type == ApecMSystemBase::DOF_SOLUTION_SYSTEM)
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
    for (int i = 0; i < ApecMSystemBase::NUM_SYSTEMS; ++i)
    {
        const auto type = static_cast<ApecMSystemBase::SystemType>(i);
        const auto &info = ApecMSystemBase::SYSTEM_INFO_MAP.at(type);
        auto &base_system = this->get_system<System>(info);

        if (auto *zeroable = dynamic_cast<ApecMSystemBase *>(&base_system))
            zeroable->zero_all();
    }
}

std::shared_ptr<ApecMEquationSystems> ApecMEquationSystems::build(ApecMMesh &mesh, APEC_M *program)
{
    if (program == nullptr)
        libmesh_error_msg("APEC_M program pointer not provided to ApecMEquationSystems::build.");
    auto equation_systems = std::make_shared<ApecMEquationSystems>(mesh);
    equation_systems->set_parent_object(program);
    equation_systems->configure_systems();
    return equation_systems;
}

std::vector<ApecMEquationSystems::RampVectorTuple> ApecMEquationSystems::collect_ramp_vectors()
{
    std::vector<RampVectorTuple> ramp_sets;

    System &main_system = this->get_system<System>("ApecMLinearSystem");
    auto *main_apec_system = dynamic_cast<ApecMSystemBase *>(&main_system);
    if (!main_apec_system)
        libmesh_error_msg("ApecMEquationSystems: ApecMLinearSystem is not an ApecMSystemBase.");
    ramp_sets.emplace_back(
        main_apec_system,
        main_system.request_vector("old_node_dof_constraints"),
        main_system.request_vector("node_dof_constraints"),
        main_system.request_vector("ramped_node_dof_constraints"));

    System &element_loads_system = this->get_system<System>("ElementLoadsSystem");
    auto *element_apec_system = dynamic_cast<ApecMSystemBase *>(&element_loads_system);
    if (!element_apec_system)
        libmesh_error_msg("ApecMEquationSystems: ElementLoadsSystem is not an ApecMSystemBase.");
    ramp_sets.emplace_back(
        element_apec_system,
        element_loads_system.request_vector("old_loads"),
        element_loads_system.request_vector("new_loads"),
        element_loads_system.solution.get());

    System &nodal_loads_system = this->get_system<System>("NodalLoadsSystem");
    auto *nodal_apec_system = dynamic_cast<ApecMSystemBase *>(&nodal_loads_system);
    if (!nodal_apec_system)
        libmesh_error_msg("ApecMEquationSystems: NodalLoadsSystem is not an ApecMSystemBase.");
    ramp_sets.emplace_back(
        nodal_apec_system,
        nodal_loads_system.request_vector("old_loads"),
        nodal_loads_system.request_vector("new_loads"),
        nodal_loads_system.solution.get());

    return ramp_sets;
}
