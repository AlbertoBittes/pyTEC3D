#include "apec_m_equation_systems.h"
#include "apec_m_system.h"
#include "apec_m_mesh.h"
#include "element_types_factory.h"
#include "vtkio_updated.h"
#include <libmesh/solver_configuration.h>
#include "solver_configuration.h"

#include <optional>
#include <set>
#include <utility>

void ApecMEquationSystems::configure_services(const std::shared_ptr<AnalysisContextView> &analysis_ctx)
{
    analysis_context_view_ = analysis_ctx;

    auto &system = this->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    auto bc_view = analysis_ctx ? analysis_ctx->template project_view<TimeStep>() : nullptr;
    auto global_constraints_object = system.get_global_constraints_object();
    global_constraints_object->set_context_view(bc_view);

    auto assembler_view = analysis_ctx ? analysis_ctx->template project_view<TimeStep>() : nullptr;
    auto assembler = system.get_assembler();
    assembler->set_context_view(assembler_view);
}

void ApecMEquationSystems::prepare_time_step_solve(TimeStep &time_step)
{
    auto &system = this->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    system.get_dof_map().swap_dof_constraints();
    system.reinit_constraints();
}

void ApecMEquationSystems::solve_equation_systems(TimeStep &time_step, bool nodal_solution)
{
    auto &system = this->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    auto assembler = system.get_assembler();

    assembler->write_element_matrices_and_rhs();

    system.solve();
    system.get_dof_map().enforce_constraints_exactly(system);
    system.update();
    system.compute_solution_time_derivative(time_step);
    if (nodal_solution)
        assembler->write_nodal_solution();
    assembler->write_element_solution();
}

void ApecMEquationSystems::print_output(TimeStep &time_step, bool picard_iteration_print)
{
    libMesh::out << "Writing output for time step " << time_step.get_index() << std::endl;

    std::set<std::string> system_names = {"ApecMLinearSystem", "NodalSolutionSystem"};

    VTKIO_Updated vtk_output(this->get_mesh(), this);
    vtk_output.write_equation_systems(,
                                      *this,
                                      &system_names);
}

void ApecMEquationSystems::update_old_local_solution()
{
    auto &system = this->get_system<ApecMLinearSystem>("ApecMLinearSystem");
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
    auto equation_systems = std::make_shared<ApecMEquationSystems>(mesh);
    equation_systems->set_parent_object(program);
    equation_systems->configure_systems();
    return equation_systems;
}
