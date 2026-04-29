#include "apec_m_equation_systems.h"
#include "apec_m_system.h"
#include "element_types_factory.h"
#include "vtkio_updated.h"
#include <libmesh/solver_configuration.h>
#include "solver_configuration.h"

#include <optional>
#include <set>
#include <utility>

void ApecMEquationSystems::solve_equation_systems(const pybind11::dict &time_step)
{
    auto &system = this->get_system<ApecMLinearSystem>("ApecMLinearSystem");

    system.prepare_time_step_solve(time_step);

    system.get_dof_map().swap_dof_constraints();
    system.reinit_constraints();

    auto assembler = system.get_assembler();
    assembler->write_element_matrices_and_rhs();

    system.solve();
    system.get_dof_map().enforce_constraints_exactly(system);
    system.update();
    system.compute_solution_time_derivative();
    assembler->write_nodal_solution();
    assembler->write_element_solution();
}

void ApecMEquationSystems::configure_systems()
{
    libMesh::out << "Building system of equations..." << std::endl;

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

void ApecMEquationSystems::initialize_systems()
{
    libMesh::out << "Distributing degrees of freedom..." << std::endl;

    for (int i = 0; i < ApecMSystemBase::NUM_SYSTEMS; ++i)
    {
        const auto type = static_cast<ApecMSystemBase::SystemType>(i);
        const auto &info = ApecMSystemBase::SYSTEM_INFO_MAP.at(type);
        auto &base_system = this->get_system<System>(info);

        if (type == ApecMSystemBase::DOF_SOLUTION_SYSTEM)
        {
            if (auto *system = dynamic_cast<ApecMLinearSystem *>(&base_system))
            {
                system->init_variables();
            }
        }
        else if (auto *system = dynamic_cast<ApecMExplicitSystem *>(&base_system))
        {
            system->init_variables();
        }
    }
}

void ApecMEquationSystems::initialize_physics()
{
    this->init();
    post_initialize_systems();
    this->print_info();
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
