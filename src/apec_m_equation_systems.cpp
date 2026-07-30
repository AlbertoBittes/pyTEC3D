/**
 * @file apec_m_equation_systems.cpp
 * @brief ApecMEquationSystems: system registration, variable initialisation, and solve orchestration.
 */
#include "apec_m_equation_systems.h"
#include "apec_m_system.h"
#include "apec_m_solver.h"
#include "element_types_factory.h"
#include "vtkio_updated.h"
#include "apec_m_logging.h"
#include <libmesh/enum_to_string.h>
#include <libmesh/petsc_linear_solver.h>
#include <libmesh/solver_configuration.h>
#include "solver_configuration.h"

#include <iomanip>
#include <numeric>
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
        constexpr int var_total_w = 14;
        constexpr int local_w = 53;

        const auto make_row = [](const std::string &system_name,
                                 const std::string &variable_name,
                                 const std::string &fe_name,
                                 const std::string &order_name,
                                 const std::string &total_dofs,
                                 const std::string &variable_dofs,
                                 const std::string &local_dofs) -> std::string
        {
            std::ostringstream row;
            row << "| " << std::left << std::setw(system_w) << fit_cell(system_name, system_w)
                << " | " << std::left << std::setw(variable_w) << fit_cell(variable_name, variable_w)
                << " | " << std::left << std::setw(fe_w) << fit_cell(fe_name, fe_w)
                << " | " << std::left << std::setw(order_w) << fit_cell(order_name, order_w)
                << " | " << std::right << std::setw(total_w) << fit_cell(total_dofs, total_w)
                << " | " << std::right << std::setw(var_total_w) << fit_cell(variable_dofs, var_total_w)
                << " | " << std::left << std::setw(local_w) << fit_cell(local_dofs, local_w)
                << " |";
            return row.str();
        };

        const std::string header = make_row("System", "Variable", "FE Type", "Order", "Total DoFs", "Variable DoFs", "Local DoFs by rank");
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
            auto &system = equation_systems.get_system<libMesh::System>(system_name);

            std::vector<libMesh::dof_id_type> local_dofs_by_rank;
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
                                  "-",
                                  local_dofs_text)
                      << "\n";
                table << border << "\n";
                continue;
            }

            for (unsigned int var = 0; var < system.n_vars(); ++var)
            {
                const libMesh::FEType fe_type = system.variable_type(var);
                const std::string fe_name = libMesh::Utility::enum_to_string<libMesh::FEFamily>(fe_type.family);
                const std::string var_order_name = libMesh::Utility::enum_to_string<libMesh::Order>(static_cast<libMesh::Order>(fe_type.order));

                std::vector<libMesh::dof_id_type> var_local_dof_indices;
                system.get_dof_map().local_variable_indices(var_local_dof_indices, var);

                const libMesh::dof_id_type var_local_dofs = static_cast<libMesh::dof_id_type>(var_local_dof_indices.size());
                std::vector<libMesh::dof_id_type> var_local_dofs_by_rank;
                system.comm().allgather(var_local_dofs, var_local_dofs_by_rank);
                const libMesh::dof_id_type var_total_dofs = std::accumulate(var_local_dofs_by_rank.begin(),
                                                                            var_local_dofs_by_rank.end(),
                                                                            libMesh::dof_id_type(0));

                table << make_row(var == 0 ? system_name : "",
                                  system.variable_name(var),
                                  fe_name,
                                  var_order_name,
                                  var == 0 ? std::to_string(system.n_dofs()) : "",
                                  std::to_string(var_total_dofs),
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
    auto &system = get_dof_solution_system();

    system.init_time_step(time_step);

    system.write_element_data();

    this->solve();
}

void ApecMEquationSystems::configure_systems()
{
    auto *solver = this->get_parent_object();
    use_nonlinear_system_ = (solver != nullptr) && solver->has_bh_curves();

    for (int i = 0; i < ApecMSystemBase::NUM_SYSTEMS; ++i)
    {
        ApecMSystemBase::SystemType type = static_cast<ApecMSystemBase::SystemType>(i);
        const std::string &sys_name = ApecMSystemBase::system_name_from_type(type);

        if (type == ApecMSystemBase::DOF_SOLUTION_SYSTEM)
        {
            if (use_nonlinear_system_)
                this->add_system<ApecMNonLinearSystem>(sys_name);
            else
                this->add_system<ApecMLinearSystem>(sys_name);
        }
        else if (type == ApecMSystemBase::NODAL_SOLUTION)
            this->add_system<ApecMNodalSolutionSystem>(sys_name);
        else if (type == ApecMSystemBase::ELEMENT_SOLUTION)
            this->add_system<ApecMElementSolutionSystem>(sys_name);
    }
}

void ApecMEquationSystems::post_initialize_systems()
{
    for (int i = 0; i < ApecMSystemBase::NUM_SYSTEMS; ++i)
    {
        const auto type = static_cast<ApecMSystemBase::SystemType>(i);
        const auto &info = ApecMSystemBase::system_name_from_type(type);

        if (type == ApecMSystemBase::DOF_SOLUTION_SYSTEM)
            get_dof_solution_system().post_init();
        else if (type == ApecMSystemBase::NODAL_SOLUTION)
        {
            auto &system = this->get_system<ApecMNodalSolutionSystem>(info);
            system.post_init();
        }
        else if (type == ApecMSystemBase::ELEMENT_SOLUTION)
        {
            auto &system = this->get_system<ApecMElementSolutionSystem>(info);
            system.post_init();
        }
    }
}

void ApecMEquationSystems::initialize_systems()
{
    apec::logging::info("Setting systems variables");

    for (int i = 0; i < ApecMSystemBase::NUM_SYSTEMS; ++i)
    {
        const auto type = static_cast<ApecMSystemBase::SystemType>(i);
        const auto &info = ApecMSystemBase::system_name_from_type(type);

        if (type == ApecMSystemBase::DOF_SOLUTION_SYSTEM)
            get_dof_solution_system().init_variables();
        else if (type == ApecMSystemBase::NODAL_SOLUTION)
        {
            auto &system = this->get_system<ApecMNodalSolutionSystem>(info);
            system.init_variables();
        }
        else if (type == ApecMSystemBase::ELEMENT_SOLUTION)
        {
            auto &system = this->get_system<ApecMElementSolutionSystem>(info);
            system.init_variables();
        }
    }
}

void ApecMEquationSystems::initialize_physics()
{
    this->init();
    post_initialize_systems();
    apec::logging::table_root(format_equation_systems_summary_table(*this));
}

bool ApecMEquationSystems::uses_nonlinear_system() const
{
    return use_nonlinear_system_;
}

ApecMDOFSystemBase &ApecMEquationSystems::get_dof_solution_system()
{
    const auto &name = ApecMSystemBase::SYSTEM_INFO_MAP.at(ApecMSystemBase::DOF_SOLUTION_SYSTEM);
    if (use_nonlinear_system_)
        return this->get_system<ApecMNonLinearSystem>(name);
    return this->get_system<ApecMLinearSystem>(name);
}

const ApecMDOFSystemBase &ApecMEquationSystems::get_dof_solution_system() const
{
    const auto &name = ApecMSystemBase::SYSTEM_INFO_MAP.at(ApecMSystemBase::DOF_SOLUTION_SYSTEM);
    if (use_nonlinear_system_)
        return this->get_system<ApecMNonLinearSystem>(name);
    return this->get_system<ApecMLinearSystem>(name);
}
