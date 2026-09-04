/**
 * @file apec_m_logging.cpp
 * @brief Process-global spdlog initialisation, MPI-rank filtering, and structured output helpers.
 *
 * The logger is initialised once (lazy singleton pattern protected by a mutex) and
 * uses a custom spdlog flag (%C) that substitutes elapsed CPU time in seconds.
 * Rank filtering is applied at the call site so that non-root ranks emit nothing
 * unless log_all_ranks() returns true.
 */
#include "apec_m_solver/apec_m_logging.h"

#include "apec_m_solver/apec_m_equation_systems.h"
#include "apec_m_solver/apec_m_system.h"

#include <libmesh/enum_to_string.h>

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/fmt/fmt.h>
#include <spdlog/pattern_formatter.h>

#include <ctime>
#include <iomanip>
#include <memory>
#include <mutex>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <libmesh/mesh_base.h>

namespace
{
    std::mutex g_logger_mutex;
    std::shared_ptr<spdlog::logger> g_logger;
    int g_rank = 0;
    int g_size = 1;
    std::clock_t g_cpu_clock_start = 0;

    class cpu_elapsed_flag_formatter final : public spdlog::custom_flag_formatter
    {
    public:
        void format(const spdlog::details::log_msg &, const std::tm &, spdlog::memory_buf_t &dest) override
        {
            const std::clock_t now = std::clock();
            const std::clock_t delta = (now >= g_cpu_clock_start) ? (now - g_cpu_clock_start) : 0;
            const double elapsed_seconds =
                static_cast<double>(delta) / static_cast<double>(CLOCKS_PER_SEC);
            fmt::format_to(std::back_inserter(dest), "{:.3f}s", elapsed_seconds);
        }

        std::unique_ptr<spdlog::custom_flag_formatter> clone() const override
        {
            return std::make_unique<cpu_elapsed_flag_formatter>();
        }
    };

    std::string to_scientific_string(libMesh::Real value)
    {
        std::ostringstream out;
        out << std::scientific << std::setprecision(6) << value;
        return out.str();
    }

    std::string fit_cell(const std::string &value, std::size_t width)
    {
        if (value.size() <= width)
            return value;
        if (width <= 3)
            return value.substr(0, width);
        return value.substr(0, width - 3) + "...";
    }

    std::string centered_line(std::string_view text, int width)
    {
        const std::string clipped =
            (static_cast<int>(text.size()) > width) ? std::string(text.substr(0, width)) : std::string(text);
        const int pad_total = width - static_cast<int>(clipped.size());
        const int pad_left = pad_total / 2;
        const int pad_right = pad_total - pad_left;
        return "|" + std::string(pad_left, ' ') + clipped + std::string(pad_right, ' ') + "|";
    }

    std::string kv_line(std::string_view key, std::string_view value, int width, int key_width)
    {
        std::ostringstream out;
        out << "| " << std::left << std::setw(key_width) << key << " : ";
        std::string row = out.str() + std::string(value);
        if (static_cast<int>(row.size()) > width + 1)
            row = row.substr(0, width + 1);
        if (static_cast<int>(row.size()) < width + 1)
            row += std::string(width + 1 - row.size(), ' ');
        row += "|";
        return row;
    }

    std::string join_elem_dimensions(const std::set<unsigned char> &dims)
    {
        std::ostringstream out;
        bool first = true;
        for (const auto dim : dims)
        {
            if (!first)
                out << ",";
            out << static_cast<unsigned int>(dim);
            first = false;
        }
        return out.str();
    }

    std::string join_elem_orders(const std::set<libMesh::Order> &orders)
    {
        std::ostringstream out;
        bool first = true;
        for (const auto order : orders)
        {
            if (!first)
                out << ",";
            out << libMesh::Utility::enum_to_string<libMesh::Order>(order);
            first = false;
        }
        return out.str();
    }

    std::string join_counts_by_rank(const std::vector<libMesh::dof_id_type> &counts)
    {
        std::ostringstream out;
        for (std::size_t rank = 0; rank < counts.size(); ++rank)
        {
            if (rank != 0)
                out << ", ";
            out << "p" << rank << "=" << counts[rank];
        }
        return out.str();
    }
} // namespace

namespace apec::logging
{
    void initialize(MPI_Comm comm, std::string_view component)
    {
        std::lock_guard<std::mutex> lock(g_logger_mutex);
        if (g_logger)
            return;

        MPI_Comm_rank(comm, &g_rank);
        MPI_Comm_size(comm, &g_size);

        std::vector<spdlog::sink_ptr> sinks;
        sinks.emplace_back(std::make_shared<spdlog::sinks::stdout_color_sink_mt>());

        const std::string logger_name =
            fmt::format("{}[r{}/{}]", component, g_rank, g_size);

        g_cpu_clock_start = std::clock();
        g_logger = std::make_shared<spdlog::logger>(logger_name, sinks.begin(), sinks.end());
        auto formatter = std::make_unique<spdlog::pattern_formatter>();
        formatter->add_flag<cpu_elapsed_flag_formatter>('C');
        formatter->set_pattern("%Y-%m-%d %H:%M:%S.%e [%^%l%$] [cpu=%C] %v");
        g_logger->set_formatter(std::move(formatter));
        g_logger->set_level(spdlog::level::info);

        g_logger->flush_on(spdlog::level::warn);
        spdlog::set_default_logger(g_logger);
    }

    void shutdown()
    {
        std::lock_guard<std::mutex> lock(g_logger_mutex);
        g_logger.reset();
        spdlog::shutdown();
    }

    int rank()
    {
        return g_rank;
    }

    int size()
    {
        return g_size;
    }

    bool log_all_ranks()
    {
        return false;
    }

    void info_root(const std::string &message)
    {
        if (g_rank == 0)
            spdlog::info("{}", message);
    }

    void print_root(const std::string &message)
    {
        if (g_rank != 0)
            return;

        std::lock_guard<std::mutex> lock(g_logger_mutex);
        fmt::print("{}\n", message);
    }

    void table_root(const std::string &table)
    {
        if (g_rank != 0)
            return;

        std::lock_guard<std::mutex> lock(g_logger_mutex);

        std::istringstream lines(table);
        std::string line;
        while (std::getline(lines, line))
            fmt::print("{}\n", line);
    }

    std::string format_solver_banner(std::string_view version, std::string_view build_mode, std::string_view build_stamp)
    {
        constexpr int inner_width = 76;
        const std::string horizontal_rule = "+" + std::string(inner_width, '-') + "+";

        std::ostringstream banner;
        banner << horizontal_rule << "\n";
        banner << centered_line("A P E C - M  S O L V E R", inner_width) << "\n";
        banner << centered_line("Low-Frequency Electromagnetics Finite Elements solver", inner_width) << "\n";
        banner << horizontal_rule << "\n";
        banner << kv_line("Version", version, inner_width, 20) << "\n";
        banner << kv_line("Build mode", build_mode, inner_width, 20) << "\n";
        banner << kv_line("Build stamp", build_stamp, inner_width, 20) << "\n";
        banner << horizontal_rule;
        return banner.str();
    }

    std::string format_mesh_summary_table(const libMesh::MeshBase &mesh)
    {
        constexpr int inner_width = 76;
        const std::string horizontal_rule = "+" + std::string(inner_width, '-') + "+";

        std::vector<libMesh::dof_id_type> local_nodes_by_rank;
        std::vector<libMesh::dof_id_type> local_elems_by_rank;
        mesh.comm().allgather(mesh.n_local_nodes(), local_nodes_by_rank);
        mesh.comm().allgather(mesh.n_local_elem(), local_elems_by_rank);

        std::ostringstream out;
        out << horizontal_rule << "\n";
        out << "| " << std::left << std::setw(inner_width) << "Mesh Summary" << "|\n";
        out << horizontal_rule << "\n";
        out << kv_line("elem dimensions", "{" + join_elem_dimensions(mesh.elem_dimensions()) + "}", inner_width, 26) << "\n";
        out << kv_line("elem default orders", "{" + join_elem_orders(mesh.elem_default_orders()) + "}", inner_width, 26) << "\n";
        out << kv_line("supported nodal order", libMesh::Utility::enum_to_string<libMesh::Order>(mesh.supported_nodal_order()), inner_width, 26) << "\n";
        out << kv_line("spatial dimension", std::to_string(mesh.spatial_dimension()), inner_width, 26) << "\n";
        out << kv_line("total nodes", std::to_string(mesh.n_nodes()), inner_width, 26) << "\n";
        out << kv_line("total elements", std::to_string(mesh.n_elem()), inner_width, 26) << "\n";
        out << kv_line("nodes by rank", join_counts_by_rank(local_nodes_by_rank), inner_width, 26) << "\n";
        out << kv_line("elements by rank", join_counts_by_rank(local_elems_by_rank), inner_width, 26) << "\n";
        out << horizontal_rule;
        return out.str();
    }

    std::string format_time_step_table(const ApecMTimeStep &time_step, std::string_view title)
    {
        std::ostringstream out;
        const auto append_row = [&out](const std::string &label, const std::string &value)
        {
            out << "| " << std::left << std::setw(20) << label
                << " | " << std::left << std::setw(20) << value << " |\n";
        };

        out << "+---------------------------------------------+\n";
        out << "| " << std::left << std::setw(44) << title << "|\n";
        out << "+----------------------+----------------------+\n";
        append_row("index", std::to_string(time_step.time_step_idx));
        append_row("time", to_scientific_string(time_step.time));
        append_row("dt", to_scientific_string(time_step.dt));
        append_row("theta", to_scientific_string(time_step.theta));
        append_row("transient_effects", time_step.transient_effects ? "true" : "false");
        out << "+----------------------+----------------------+";
        return out.str();
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

    std::string format_constraints_summary_table(std::size_t time_integration_count,
                                                 std::size_t dof_constraint_count,
                                                 std::size_t coupling_count,
                                                 std::size_t total_count,
                                                 std::size_t applied_count,
                                                 const std::vector<ConstraintVariableStat> &variable_stats,
                                                 std::size_t unknown_total_count,
                                                 std::size_t unknown_applied_count)
    {
        std::ostringstream out;
        const auto fit_label = [](const std::string &label) -> std::string
        {
            if (label.size() <= 29)
                return label;
            return label.substr(0, 26) + "...";
        };

        const auto append_row = [&out](const std::string &label, std::size_t value)
        {
            out << "| " << std::left << std::setw(29) << label
                << " | " << std::left << std::setw(19) << value << " |\n";
        };

        const auto append_value_row = [&out, &fit_label](const std::string &label, const std::string &value)
        {
            out << "| " << std::left << std::setw(29) << fit_label(label)
                << " | " << std::left << std::setw(19) << value << " |\n";
        };

        out << "+-----------------------------------------------------+\n";
        out << "| " << std::left << std::setw(52) << "Constraints Summary" << "|\n";
        out << "+-------------------------------+---------------------+\n";
        append_row("time integration constraints", time_integration_count);
        append_row("dof constraints", dof_constraint_count);
        append_row("coupling constraints", coupling_count);
        append_row("total constraints", total_count);
        append_row("applied constraints", applied_count);

        out << "+-------------------------------+---------------------+\n";
        append_value_row("per-variable (total/applied)", "");
        for (const auto &var_stat : variable_stats)
        {
            append_value_row("  " + var_stat.variable_name,
                             std::to_string(var_stat.total_count) + " / " + std::to_string(var_stat.applied_count));
        }
        if (unknown_total_count > 0 || unknown_applied_count > 0)
        {
            append_value_row("  <unknown>",
                             std::to_string(unknown_total_count) + " / " + std::to_string(unknown_applied_count));
        }
        out << "+-------------------------------+---------------------+";
        return out.str();
    }

    std::string format_petsc_solver_configuration_table(std::string_view backend_name,
                                                        std::string_view solver_class,
                                                        std::string_view ksp_name,
                                                        std::string_view pc_name)
    {
        std::ostringstream table;
        table << "+------------------------------------------------------------+\n";
        table << "| PETSc Solver Configuration                                 |\n";
        table << "+----------------+-------------------------------------------+\n";
        table << fmt::format("| {:<14} | {:<41} |\n", "Backend", backend_name);
        table << fmt::format("| {:<14} | {:<41} |\n", "Type", solver_class);
        table << fmt::format("| {:<14} | {:<41} |\n", "KSP", ksp_name);
        table << fmt::format("| {:<14} | {:<41} |\n", "Preconditioner", pc_name);
        table << "+------------------------------------------------------------+";
        return table.str();
    }
} // namespace apec::logging
