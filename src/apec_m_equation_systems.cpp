#include "apec_m_equation_systems.h"
#include "apec_m_system.h"
#include "element_types_factory.h"
#include "vtkio_updated.h"
#include "apec_m_logging.h"
#include <libmesh/enum_to_string.h>
#include <libmesh/petsc_linear_solver.h>
#include <libmesh/solver_configuration.h>
#include "solver_configuration.h"

#include <iomanip>
#include <optional>
#include <sstream>
#include <set>
#include <utility>

namespace
{
    std::string fit_cell(const std::string &value, std::size_t width)
    {
        if (value.size() <= width)
            return value;
        if (width <= 3)
            return value.substr(0, width);
        return value.substr(0, width - 3) + "...";
    }

    std::string format_equation_systems_summary_table(ApecMEquationSystems &equation_systems)
    {
        constexpr int system_w = 22;
        constexpr int variable_w = 16;
        constexpr int fe_w = 16;
        constexpr int order_w = 10;
        constexpr int total_w = 12;
        constexpr int local_w = 70;

        const auto make_row = [](const std::string &system_name,
                                 const std::string &variable_name,
                                 const std::string &fe_name,
                                 const std::string &order_name,
                                 const std::string &total_dofs,
                                 const std::string &local_dofs) -> std::string
        {
            std::ostringstream row;
            row << "| " << std::left << std::setw(system_w) << fit_cell(system_name, system_w)
                << " | " << std::left << std::setw(variable_w) << fit_cell(variable_name, variable_w)
                << " | " << std::left << std::setw(fe_w) << fit_cell(fe_name, fe_w)
                << " | " << std::left << std::setw(order_w) << fit_cell(order_name, order_w)
                << " | " << std::right << std::setw(total_w) << fit_cell(total_dofs, total_w)
                << " | " << std::left << std::setw(local_w) << fit_cell(local_dofs, local_w)
                << " |";
            return row.str();
        };

        const std::string header = make_row("System", "Variable", "FE Type", "Order", "Total DoFs", "Local DoFs by rank");
        const std::string border = "+" + std::string(header.size() - 2, '-') + "+";

        std::ostringstream table;
        table << border << "\n";
        table << "| " << std::left << std::setw(static_cast<int>(header.size() - 4))
              << "Systems Summary"
              << " |\n";
        table << border << "\n";
        table << header << "\n";
        table << border << "\n";

        for (int i = 0; i < ApecMSystemBase::NUM_SYSTEMS; ++i)
        {
            const auto type = static_cast<ApecMSystemBase::SystemType>(i);
            const auto &system_name = ApecMSystemBase::SYSTEM_INFO_MAP.at(type);
            auto &system = equation_systems.get_system<System>(system_name);

            std::vector<dof_id_type> local_dofs_by_rank;
            system.comm().allgather(system.n_local_dofs(), local_dofs_by_rank);

            std::ostringstream local_dofs_stream;
            for (std::size_t rank = 0; rank < local_dofs_by_rank.size(); ++rank)
            {
                if (rank != 0)
                    local_dofs_stream << ", ";
                local_dofs_stream << "p" << rank << "=" << local_dofs_by_rank[rank];
            }
            const std::string local_dofs_text = local_dofs_stream.str();

            if (system.n_vars() == 0)
            {
                table << make_row(system_name,
                                  "-",
                                  "-",
                                  "-",
                                  std::to_string(system.n_dofs()),
                                  local_dofs_text)
                      << "\n";
                table << border << "\n";
                continue;
            }

            for (unsigned int var = 0; var < system.n_vars(); ++var)
            {
                const FEType fe_type = system.variable_type(var);
                const std::string fe_name = Utility::enum_to_string<FEFamily>(fe_type.family);
                const std::string var_order_name = Utility::enum_to_string<Order>(static_cast<Order>(fe_type.order));

                table << make_row(var == 0 ? system_name : "",
                                  system.variable_name(var),
                                  fe_name,
                                  var_order_name,
                                  var == 0 ? std::to_string(system.n_dofs()) : "",
                                  var == 0 ? local_dofs_text : "")
                      << "\n";
            }
            table << border << "\n";
        }

        return table.str();
    }

}

void ApecMEquationSystems::solve_equation_systems(const pybind11::dict &time_step)
{
    apec::logging::info("Preparing time-step solution");
    auto &system = this->get_system<ApecMLinearSystem>(ApecMSystemBase::SYSTEM_INFO_MAP.at(ApecMSystemBase::DOF_SOLUTION_SYSTEM));

    system.prepare_time_step_solve(time_step);

    system.get_dof_map().swap_dof_constraints();
    system.reinit_constraints();

    auto assembler = system.get_assembler();
    assembler->write_element_matrices_and_rhs();

    system.solve();
    system.get_dof_map().enforce_constraints_exactly(system);
    system.update();
    system.compute_solution_time_derivative();
    apec::logging::info("Finished system solve");
    assembler->write_nodal_solution();
    apec::logging::info("Nodal solution updated");
    assembler->write_element_solution();
    apec::logging::info("Element solution updated");
}

void ApecMEquationSystems::configure_systems()
{
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
    apec::logging::info("Setting systems variables");

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
    apec::logging::table_root(format_equation_systems_summary_table(*this));
}

void ApecMEquationSystems::zero_all_systems()
{
    apec::logging::info("Zeroing all equation systems");
    for (int i = 0; i < ApecMSystemBase::NUM_SYSTEMS; ++i)
    {
        const auto type = static_cast<ApecMSystemBase::SystemType>(i);
        const auto &info = ApecMSystemBase::SYSTEM_INFO_MAP.at(type);
        auto &base_system = this->get_system<System>(info);

        if (auto *zeroable = dynamic_cast<ApecMSystemBase *>(&base_system))
            zeroable->zero_all();
    }
}
