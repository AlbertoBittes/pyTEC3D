/**
 * @file apec_m_logging.h
 * @brief MPI-aware structured logging interface for the APEC-M solver.
 *
 * Wraps spdlog with rank-filtering helpers so that, by default, only the
 * root rank (rank 0) emits messages.  All free functions forward to the
 * process-wide spdlog default logger initialised by apec::logging::initialize().
 * Use apec::logging::error_and_abort() for unrecoverable errors; it logs the
 * message on all ranks and then throws in debug builds so Python/debuggers can
 * show the original failure. In release builds, set APEC_M_HARD_ABORT=1 to
 * force libmesh_error_msg() hard-abort semantics.
 */
#pragma once

#include <mpi.h>
#include <spdlog/spdlog.h>
#include <libmesh/libmesh_common.h>

#include <cstdlib>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "apec_m_time_step.h"

class ApecMEquationSystems;
namespace libMesh
{
    class MeshBase;
}

namespace apec::logging
{
    struct ConstraintVariableStat
    {
        std::string variable_name;
        std::size_t total_count = 0;
        std::size_t applied_count = 0;
    };

    void initialize(MPI_Comm comm, std::string_view component = "apec_m");
    void shutdown();

    int rank();
    int size();
    bool log_all_ranks();

    void info_root(const std::string &message);
    void print_root(const std::string &message);
    void table_root(const std::string &table);
    std::string format_solver_banner(std::string_view version, std::string_view build_mode, std::string_view build_stamp);
    std::string format_mesh_summary_table(const libMesh::MeshBase &mesh);
    std::string format_time_step_table(const ApecMTimeStep &time_step, std::string_view title = "Time Step Information");
    std::string format_equation_systems_summary_table(ApecMEquationSystems &equation_systems);
    std::string format_constraints_summary_table(std::size_t time_integration_count,
                                                 std::size_t dof_constraint_count,
                                                 std::size_t coupling_count,
                                                 std::size_t total_count,
                                                 std::size_t applied_count,
                                                 const std::vector<ConstraintVariableStat> &variable_stats,
                                                 std::size_t unknown_total_count,
                                                 std::size_t unknown_applied_count);
    std::string format_petsc_solver_configuration_table(std::string_view backend_name,
                                                        std::string_view solver_class,
                                                        std::string_view ksp_name,
                                                        std::string_view pc_name);

    inline bool should_log_this_rank()
    {
        return log_all_ranks() || rank() == 0;
    }

    template <typename... Args>
    inline void trace(spdlog::format_string_t<Args...> fmt, Args &&...args)
    {
        (void)fmt;
        (void)sizeof...(args);
    }

    template <typename... Args>
    inline void info(spdlog::format_string_t<Args...> fmt, Args &&...args)
    {
        if (should_log_this_rank())
            spdlog::info(fmt, std::forward<Args>(args)...);
    }

    template <typename... Args>
    inline void info_all_ranks(spdlog::format_string_t<Args...> fmt, Args &&...args)
    {
        spdlog::info(fmt, std::forward<Args>(args)...);
    }

    template <typename... Args>
    inline void warn(spdlog::format_string_t<Args...> fmt, Args &&...args)
    {
        if (should_log_this_rank())
            spdlog::warn(fmt, std::forward<Args>(args)...);
    }

    template <typename... Args>
    inline void warn_all_ranks(spdlog::format_string_t<Args...> fmt, Args &&...args)
    {
        spdlog::warn(fmt, std::forward<Args>(args)...);
    }

    template <typename... Args>
    inline void error(spdlog::format_string_t<Args...> fmt, Args &&...args)
    {
        if (should_log_this_rank())
            spdlog::error(fmt, std::forward<Args>(args)...);
    }

    template <typename... Args>
    inline void error_all_ranks(spdlog::format_string_t<Args...> fmt, Args &&...args)
    {
        spdlog::error(fmt, std::forward<Args>(args)...);
    }

    template <typename... Args>
    inline void critical(spdlog::format_string_t<Args...> fmt, Args &&...args)
    {
        if (should_log_this_rank())
            spdlog::critical(fmt, std::forward<Args>(args)...);
    }

    template <typename... Args>
    inline void critical_all_ranks(spdlog::format_string_t<Args...> fmt, Args &&...args)
    {
        spdlog::critical(fmt, std::forward<Args>(args)...);
    }

    [[noreturn]] inline void error_and_abort(const std::string &message)
    {
        error_all_ranks("[rank {}/{}] {}", rank(), size(), message);
        throw std::runtime_error(message);
    }
} // namespace apec::logging
