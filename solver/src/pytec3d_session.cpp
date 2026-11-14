/**
 * @file pytec3d_session.cpp
 * @brief pyTEC3DSession implementation: runtime lifecycle, model loading, and solve dispatch.
 */
#include "pytec3d_session.h"

#include "git_commit.h"
#include "logging.h"
#include "pytec3d_model.h"
#include "report_tables.h"

#include <mfem.hpp>
#include <mpi.h>
#include <omp.h>
#include <petscsys.h>

#ifdef _WIN32
#include <winsock2.h>
#else
#include <unistd.h>
#endif

#include <climits>
#include <mutex>
#include <string>

#define pyTEC3D_STR_IMPL(x) #x
#define pyTEC3D_STR(x) pyTEC3D_STR_IMPL(x)

namespace
{
    enum class RuntimeState
    {
        Uninitialized,
        Initialized,
        Closed
    };

    /** @brief The process-wide runtime: MPI, PETSc (the linear solvers), hypre (the parallel matrices of MFEM). */
    class Runtime
    {
    public:
        explicit Runtime(int n_threads)
        {
            int initialized = 0;
            MPI_Initialized(&initialized);
            if (!initialized)
            {
                int provided = 0;
                // NULL-terminated: MPI reads argv up to the null pointer when it starts without mpirun (singleton).
                static char name[] = "pytec3d_solver";
                static char *argv_storage[] = {name, nullptr};
                int argc = 1;
                char **argv = argv_storage;
                MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
                owns_mpi_ = true;
            }
            PetscInitializeNoArguments();
            mfem::Hypre::Init();
            omp_set_num_threads(n_threads);
        }
        ~Runtime()
        {
            mfem::Hypre::Finalize();
            PetscFinalize();
            if (owns_mpi_)
                MPI_Finalize();
        }

    private:
        bool owns_mpi_ = false;
    };

    /** @brief The process-wide runtime, shared by the sessions. */
    class RuntimeManager
    {
    public:
        static RuntimeManager &instance()
        {
            static RuntimeManager manager;
            return manager;
        }

        std::shared_ptr<Runtime> acquire(int n_threads)
        {
            std::lock_guard<std::mutex> guard(mutex_);

            if (state_ == RuntimeState::Closed)
                pytec3d::logging::error_and_abort("The pyTEC3D runtime has been shut down in this process.");

            if (state_ == RuntimeState::Uninitialized)
            {
                runtime_ = std::make_shared<Runtime>(n_threads);
                pytec3d::logging::initialize(MPI_COMM_WORLD, "pytec3d");
                pytec3d::logging::info("Runtime started: {} MPI rank(s), {} thread(s) per rank.",
                                       pytec3d::logging::size(), n_threads);
                n_threads_ = n_threads;
                state_ = RuntimeState::Initialized;
            }
            else if (n_threads != n_threads_)
                pytec3d::logging::warn("The runtime runs with {} threads per MPI process (set by the first session): "
                                       "the {} requested are ignored.",
                                       n_threads_, n_threads);

            ++active_sessions_;
            return runtime_;
        }

        void release()
        {
            std::lock_guard<std::mutex> guard(mutex_);
            if (active_sessions_ > 0)
                --active_sessions_;
        }

        void shutdown()
        {
            std::lock_guard<std::mutex> guard(mutex_);

            if (state_ != RuntimeState::Initialized)
            {
                state_ = RuntimeState::Closed;
                return;
            }

            if (active_sessions_ != 0)
                pytec3d::logging::error_and_abort("Cannot shut down the pyTEC3D runtime while sessions are open.");

            pytec3d::logging::shutdown();
            runtime_.reset();
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
        std::shared_ptr<Runtime> runtime_;
        int n_threads_ = 1;
        int active_sessions_ = 0;
        bool banner_printed_ = false;
    };

    const char *pytec3d_version_string()
    {
#ifdef pyTEC3D_VERSION
        return pyTEC3D_STR(pyTEC3D_VERSION);
#else
        return "dev";
#endif
    }

    const char *build_mode_string()
    {
#ifdef pyTEC3D_BUILD_TYPE
        return pyTEC3D_BUILD_TYPE;
#elif defined(NDEBUG)
        return "Release";
#else
        return "Debug";
#endif
    }

    const char *git_commit_string()
    {
#ifdef pyTEC3D_GIT_COMMIT
        return pyTEC3D_GIT_COMMIT;
#else
        return "unknown";
#endif
    }

    std::string host_name()
    {
        char name[256] = {};
        return gethostname(name, sizeof(name) - 1) == 0 ? std::string(name) : std::string("unknown");
    }

    void print_solver_banner_once()
    {
        if (!RuntimeManager::instance().consume_banner_needed())
            return;

        const std::string stamp = std::string(__DATE__) + " " + __TIME__;
        pytec3d::logging::table_root(pytec3d::report::format_solver_banner(
            {pytec3d_version_string(), build_mode_string(), stamp, git_commit_string(), host_name(), pytec3d::logging::size(),
             omp_get_max_threads()}));
    }
} // namespace

pyTEC3DSession::pyTEC3DSession(int num_threads_per_mpi_process)
{
    if (num_threads_per_mpi_process < 1)
        pytec3d::logging::error_and_abort("num_threads_per_mpi_process must be >= 1.");

    runtime_ = RuntimeManager::instance().acquire(num_threads_per_mpi_process);
    print_solver_banner_once();
}

pyTEC3DSession::~pyTEC3DSession()
{
    close();
}

void pyTEC3DSession::close()
{
    if (closed_)
        return;

    model_.reset();
    runtime_.reset();
    RuntimeManager::instance().release();
    closed_ = true;
}

void pyTEC3DSession::shutdown_runtime()
{
    RuntimeManager::instance().shutdown();
}

void pyTEC3DSession::load_model(const pytec3d::ModelInput &input)
{
    if (closed_)
        pytec3d::logging::error_and_abort("The solver session has been closed and cannot be reused.");

    model_.reset();
    auto model = std::make_unique<pyTEC3DModel>();
    model->load(MPI_COMM_WORLD, input);
    model_ = std::move(model);
}

void pyTEC3DSession::set_arrays(const pytec3d::ModelArrays &arrays)
{
    require_loaded();
    model_->set_arrays(arrays);
}

void pyTEC3DSession::solve(const TimeStepInput &time_step, const SolverSettings &settings)
{
    require_loaded();
    model_->solve(time_step, settings);
}

void pyTEC3DSession::write_system_output(const std::string &system_name)
{
    require_loaded();
    model_->write_system_output(system_name);
}

const StepReport &pyTEC3DSession::last_step_report() const
{
    require_loaded();
    return model_->last_step_report();
}

bool pyTEC3DSession::is_loaded() const
{
    return model_ != nullptr;
}

void pyTEC3DSession::require_loaded() const
{
    if (!is_loaded())
        pytec3d::logging::error_and_abort("The solver session has no model. Call load_model() first.");
}
