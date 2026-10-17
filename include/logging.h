/**
 * @file logging.h
 * @brief MPI-aware logging for the pyTEC3D solver.
 *
 * Wraps spdlog with rank-filtering helpers so that, by default, only the root rank (rank 0) emits messages.
 * All free functions forward to the process-wide spdlog default logger initialised by initialize().
 * error_and_abort() is for unrecoverable errors: it logs the message on all ranks, then calls MPI_Abort when
 * running on several ranks (an error raised on some ranks only would otherwise leave the others waiting in the
 * next collective call) or throws std::runtime_error in a serial run, which reaches Python as an exception.
 */
#pragma once

#include <mpi.h>
#include <spdlog/spdlog.h>

#include <chrono>
#include <cstddef>
#include <string>
#include <string_view>
#include <utility>

namespace pytec3d::logging
{
    void initialize(MPI_Comm comm, std::string_view component = "pytec3d");
    void shutdown();

    int rank();
    int size();
    bool log_all_ranks();

    /**
     * @brief Also write the log (messages and tables, rank 0) to @p path, appended. Can be called before
     * initialize(); an empty path stops writing to a file.
     */
    void set_log_file(const std::string &path);

    /** @brief Warnings and errors logged on this rank since initialize(). */
    std::pair<std::size_t, std::size_t> message_counts();

    void info_root(const std::string &message);
    void print_root(const std::string &message);
    void table_root(const std::string &table);

    /** @brief Logs @p message on all ranks, then MPI_Abort on several ranks or throws std::runtime_error. */
    [[noreturn]] void error_and_abort(const std::string &message);

    /**
     * @brief Logs, on the root rank, the wall time of a phase when it goes out of scope: "<label>: 1.234 s", or a
     * warning that the phase was interrupted when an exception leaves it.
     *
     * For phases (load, assembly, solve), not for loops over elements.
     */
    class ScopedTimer
    {
    public:
        /** @param seconds If given, receives the elapsed time when the timer goes out of scope. */
        explicit ScopedTimer(std::string label, double *seconds = nullptr);
        ~ScopedTimer();
        ScopedTimer(const ScopedTimer &) = delete;
        ScopedTimer &operator=(const ScopedTimer &) = delete;

        /** @brief Seconds since the timer started. */
        double elapsed() const;

    private:
        std::string label_;
        std::chrono::steady_clock::time_point start_;
        int exceptions_ = 0; ///< std::uncaught_exceptions() at the start.
        double *seconds_ = nullptr;
    };

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
} // namespace pytec3d::logging
