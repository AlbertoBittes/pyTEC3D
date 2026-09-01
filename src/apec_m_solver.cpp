/**
 * @file apec_m_solver.cpp
 * @brief ApecMSolver implementation: runtime lifecycle, model loading, and solve dispatch.
 */
#include "apec_m_solver.h"

#include "apec_m_logging.h"
#include "apec_m_model.h"

#include <libmesh/libmesh.h>

#include <mpi.h>

#include <exception>
#include <iostream>
#include <mutex>
#include <string>
#include <utility>

namespace
{
    [[noreturn]] void log_and_rethrow(const char *context, const std::exception &ex)
    {
        try
        {
            apec::logging::error_all_ranks("[rank {}/{}] {}: {}",
                                           apec::logging::rank(),
                                           apec::logging::size(),
                                           context,
                                           ex.what());
        }
        catch (...)
        {
        }

        std::cerr << "[apec_m] " << context << ": " << ex.what() << std::endl;
        throw;
    }

    [[noreturn]] void log_unknown_and_rethrow(const char *context)
    {
        try
        {
            apec::logging::error_all_ranks("[rank {}/{}] {}: unknown exception",
                                           apec::logging::rank(),
                                           apec::logging::size(),
                                           context);
        }
        catch (...)
        {
        }

        std::cerr << "[apec_m] " << context << ": unknown exception" << std::endl;
        throw;
    }

    enum class RuntimeState
    {
        Uninitialized,
        Initialized,
        Closed
    };

    class RuntimeManager
    {
    public:
        static RuntimeManager &instance()
        {
            static RuntimeManager manager;
            return manager;
        }

        std::shared_ptr<libMesh::LibMeshInit> acquire(int argc, char **argv, int n_threads)
        {
            std::lock_guard<std::mutex> guard(mutex_);

            if (state_ == RuntimeState::Closed)
                apec::logging::error_and_abort("ApecMSolver runtime is permanently closed in this process.");

            if (state_ == RuntimeState::Uninitialized)
            {
                try
                {
                    static const char *default_argv[] = {"apec_m_solver"};
                    const bool has_user_args = (argc > 0 && argv != nullptr);
                    const int init_argc = has_user_args ? argc : 1;
                    const char *const *init_argv = has_user_args
                                                       ? const_cast<const char *const *>(argv)
                                                       : default_argv;

                    process_libmesh_init_ = std::make_shared<libMesh::LibMeshInit>(
                        init_argc,
                        init_argv,
                        MPI_COMM_WORLD,
                        n_threads);
                    apec::logging::initialize(MPI_COMM_WORLD, "apec_m");
                    logging_initialized_ = true;
                    state_ = RuntimeState::Initialized;
                }
                catch (const std::exception &ex)
                {
                    log_and_rethrow("RuntimeManager::acquire failed", ex);
                }
                catch (...)
                {
                    log_unknown_and_rethrow("RuntimeManager::acquire failed");
                }
            }

            ++active_instances_;
            return process_libmesh_init_;
        }

        void release_instance()
        {
            std::lock_guard<std::mutex> guard(mutex_);
            if (active_instances_ > 0)
                --active_instances_;
        }

        void shutdown_runtime()
        {
            std::lock_guard<std::mutex> guard(mutex_);

            if (state_ == RuntimeState::Closed)
                return;

            if (active_instances_ != 0)
                apec::logging::error_and_abort("Cannot shutdown ApecMSolver runtime while solver instances are still active.");

            if (logging_initialized_)
            {
                apec::logging::shutdown();
                logging_initialized_ = false;
            }

            process_libmesh_init_.reset();
            state_ = RuntimeState::Closed;
        }

        bool consume_banner_needed()
        {
            std::lock_guard<std::mutex> guard(mutex_);
            if (banner_printed_)
                return false;
            banner_printed_ = true;
            return true;
        }

    private:
        RuntimeManager() = default;

        std::mutex mutex_;
        RuntimeState state_ = RuntimeState::Uninitialized;
        std::shared_ptr<libMesh::LibMeshInit> process_libmesh_init_;
        int active_instances_ = 0;
        bool logging_initialized_ = false;
        bool banner_printed_ = false;
    };

#define APEC_M_STR_IMPL(x) #x
#define APEC_M_STR(x) APEC_M_STR_IMPL(x)

    const char *apec_m_version_string()
    {
#ifdef APEC_M_VERSION
        return APEC_M_STR(APEC_M_VERSION);
#else
        return "dev";
#endif
    }

    const char *build_mode_string()
    {
#ifdef APEC_M_BUILD_TYPE
        return APEC_M_BUILD_TYPE;
#elif defined(NDEBUG)
        return "Release";
#else
        return "Debug";
#endif
    }

    void print_solver_banner_once()
    {
        if (!RuntimeManager::instance().consume_banner_needed())
            return;

        apec::logging::table_root(apec::logging::format_solver_banner(
            apec_m_version_string(),
            build_mode_string(),
            std::string(__DATE__) + " " + __TIME__));
    }

} // namespace

void ApecMSolver::ensure_libmesh_initialized()
{
    libmesh_init = RuntimeManager::instance().acquire(argc, argv, num_threads_per_mpi_process_);
    print_solver_banner_once();
}

ApecMSolver::ApecMSolver(int num_threads_per_mpi_process)
{
    if (num_threads_per_mpi_process < 1)
        apec::logging::error_and_abort("num_threads_per_mpi_process must be >= 1.");

    num_threads_per_mpi_process_ = num_threads_per_mpi_process;
    ensure_libmesh_initialized();
    closed_ = false;
    apec::logging::info("Solver created successfully");
}

ApecMSolver::~ApecMSolver()
{
    close();
}

void ApecMSolver::close()
{
    if (closed_)
        return;

    try
    {
        if (model_)
        {
            model_->clear();
            model_.reset();
        }

        loaded = false;
        closed_ = true;

        libmesh_init.reset();
        RuntimeManager::instance().release_instance();
    }
    catch (const std::exception &ex)
    {
        log_and_rethrow("ApecMSolver::close failed", ex);
    }
    catch (...)
    {
        log_unknown_and_rethrow("ApecMSolver::close failed");
    }
}

void ApecMSolver::shutdown_runtime()
{
    try
    {
        RuntimeManager::instance().shutdown_runtime();
    }
    catch (const std::exception &ex)
    {
        log_and_rethrow("ApecMSolver::shutdown_runtime failed", ex);
    }
    catch (...)
    {
        log_unknown_and_rethrow("ApecMSolver::shutdown_runtime failed");
    }
}

void ApecMSolver::load_model(pybind11::object apec_object)
{
    try
    {
        if (closed_)
            apec::logging::error_and_abort("Solver instance has been closed and cannot be reused.");

        if (!model_)
            model_ = std::make_shared<ApecMModel>();

        model_->load_into_solver(libmesh_init, std::move(apec_object));
        loaded = true;
    }
    catch (const std::exception &ex)
    {
        log_and_rethrow("ApecMSolver::load_model failed", ex);
    }
    catch (...)
    {
        log_unknown_and_rethrow("ApecMSolver::load_model failed");
    }
}

void ApecMSolver::solve()
{
    try
    {
        if (!loaded || !model_ || !model_->is_loaded())
            apec::logging::error_and_abort("Solver has no loaded model. Call load_model() first.");

        model_->solve();
    }
    catch (const std::exception &ex)
    {
        log_and_rethrow("ApecMSolver::solve failed", ex);
    }
    catch (...)
    {
        log_unknown_and_rethrow("ApecMSolver::solve failed");
    }
}

void ApecMSolver::write_system_to_vtk_grid(const std::string &system_name)
{
    try
    {
        if (!loaded || !model_ || !model_->is_loaded())
            apec::logging::error_and_abort("Solver is not fully initialized. Call load_model() first.");

        model_->write_system_to_vtk_grid(system_name);
    }
    catch (const std::exception &ex)
    {
        log_and_rethrow("ApecMSolver::write_system_to_vtk_grid failed", ex);
    }
    catch (...)
    {
        log_unknown_and_rethrow("ApecMSolver::write_system_to_vtk_grid failed");
    }
}

bool ApecMSolver::is_loaded() const
{
    return loaded && model_ && model_->is_loaded();
}
