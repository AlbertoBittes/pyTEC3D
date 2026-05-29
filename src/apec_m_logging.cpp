#include "apec_m_logging.h"

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/fmt/fmt.h>
#include <spdlog/pattern_formatter.h>

#include <ctime>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

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
} // namespace apec::logging
