#pragma once

#include <mpi.h>
#include <spdlog/spdlog.h>

#include <string>
#include <string_view>
#include <utility>

namespace apec::logging
{
    void initialize(MPI_Comm comm, std::string_view component = "apec_m");
    void shutdown();

    int rank();
    int size();
    bool log_all_ranks();

    void info_root(const std::string &message);
    void print_root(const std::string &message);
    void table_root(const std::string &table);

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
} // namespace apec::logging
