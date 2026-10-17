/**
 * @file logging.cpp
 * @brief Process-global spdlog initialisation and MPI-rank filtering.
 *
 * The logger is initialised once (protected by a mutex) and uses a custom spdlog flag (%C) that substitutes
 * elapsed CPU time in seconds. Rank filtering is applied at the call site so that non-root ranks emit nothing
 * unless log_all_ranks() returns true.
 */
#include "logging.h"

#include <spdlog/fmt/fmt.h>
#include <spdlog/pattern_formatter.h>
#include <spdlog/sinks/base_sink.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <algorithm>
#include <cstdio>
#include <ctime>
#include <exception>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

namespace
{
    /** @brief Counts the warnings and errors logged on this rank. */
    class CountingSink final : public spdlog::sinks::base_sink<std::mutex>
    {
    public:
        std::size_t warnings = 0;
        std::size_t errors = 0;

    protected:
        void sink_it_(const spdlog::details::log_msg &message) override
        {
            if (message.level == spdlog::level::warn)
                ++warnings;
            else if (message.level >= spdlog::level::err)
                ++errors;
        }
        void flush_() override {}
    };

    int g_rank = 0;
    int g_size = 1;
    std::clock_t g_cpu_clock_start = 0;
    std::mutex g_logger_mutex;
    std::shared_ptr<spdlog::logger> g_logger;
    std::shared_ptr<CountingSink> g_counter;
    std::shared_ptr<spdlog::sinks::basic_file_sink_mt> g_file_sink;
    std::string g_log_file;
    std::FILE *g_table_file = nullptr; ///< The log file, for the tables (written raw).

    /** @brief Attach the log file to the logger (rank 0); the caller holds g_logger_mutex. */
    void attach_log_file()
    {
        if (!g_logger || g_rank != 0)
            return;
        auto &sinks = g_logger->sinks();
        if (g_file_sink)
            sinks.erase(std::remove(sinks.begin(), sinks.end(), g_file_sink), sinks.end());
        g_file_sink.reset();
        if (g_table_file)
            std::fclose(g_table_file);
        g_table_file = nullptr;
        if (g_log_file.empty())
            return;
        g_file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(g_log_file, false);
        g_file_sink->set_pattern("%Y-%m-%d %H:%M:%S.%e [%l] %v");
        sinks.push_back(g_file_sink);
        g_table_file = std::fopen(g_log_file.c_str(), "a");
    }


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
} // namespace

namespace pytec3d::logging
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
        g_counter = std::make_shared<CountingSink>();
        sinks.emplace_back(g_counter);

        const std::string logger_name =
            fmt::format("{}[r{}/{}]", component, g_rank, g_size);

        g_cpu_clock_start = std::clock();
        g_logger = std::make_shared<spdlog::logger>(logger_name, sinks.begin(), sinks.end());
        auto formatter = std::make_unique<spdlog::pattern_formatter>();
        formatter->add_flag<cpu_elapsed_flag_formatter>('C');
        formatter->set_pattern("%Y-%m-%d %H:%M:%S.%e [%^%l%$] [cpu=%C] %v");
        g_logger->set_formatter(std::move(formatter));
        g_logger->set_level(spdlog::level::info);

        g_logger->flush_on(spdlog::level::info);
        spdlog::set_default_logger(g_logger);
        attach_log_file();
    }

    void set_log_file(const std::string &path)
    {
        std::lock_guard<std::mutex> lock(g_logger_mutex);
        g_log_file = path;
        attach_log_file();
    }

    std::pair<std::size_t, std::size_t> message_counts()
    {
        std::lock_guard<std::mutex> lock(g_logger_mutex);
        return g_counter ? std::pair{g_counter->warnings, g_counter->errors} : std::pair<std::size_t, std::size_t>{0, 0};
    }

    void shutdown()
    {
        std::lock_guard<std::mutex> lock(g_logger_mutex);
        if (g_table_file)
            std::fclose(g_table_file);
        g_table_file = nullptr;
        g_file_sink.reset();
        g_counter.reset();
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
        std::fflush(stdout);
        if (g_table_file)
        {
            fmt::print(g_table_file, "{}\n", message);
            std::fflush(g_table_file);
        }
    }

    void table_root(const std::string &table)
    {
        if (g_rank != 0)
            return;

        std::lock_guard<std::mutex> lock(g_logger_mutex);

        std::istringstream lines(table);
        std::string line;
        while (std::getline(lines, line))
        {
            fmt::print("{}\n", line);
            if (g_table_file)
                fmt::print(g_table_file, "{}\n", line);
        }
        // Flushed so that the Python log, written to the same streams, does not land inside the table.
        std::fflush(stdout);
        if (g_table_file)
            std::fflush(g_table_file);
    }

    ScopedTimer::ScopedTimer(std::string label, double *seconds)
        : label_(std::move(label)), start_(std::chrono::steady_clock::now()), exceptions_(std::uncaught_exceptions()),
          seconds_(seconds)
    {
    }

    ScopedTimer::~ScopedTimer()
    {
        if (seconds_)
            *seconds_ = elapsed();
        if (std::uncaught_exceptions() > exceptions_)
            warn("Interrupted by an error after {:.3f} s: {}", elapsed(), label_);
        else
            info("{}: {:.3f} s", label_, elapsed());
    }

    double ScopedTimer::elapsed() const
    {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - start_).count();
    }

    [[noreturn]] void error_and_abort(const std::string &message)
    {
        error_all_ranks("[rank {}/{}] {}", rank(), size(), message);
        if (size() > 1)
        {
            // A failure on some ranks only would leave the others waiting in the next collective call.
            spdlog::default_logger()->flush();
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        throw std::runtime_error(message);
    }
} // namespace pytec3d::logging
