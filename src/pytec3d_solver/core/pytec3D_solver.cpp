/**
 * @file pytec3d_solver.cpp
 * @brief pyTEC3DSolver implementation: runtime lifecycle, model loading, and solve dispatch.
 */
#include "pytec3d_solver.h"

#include "logging.h"
#include "pytec3d_model.h"
#include "registries.h"

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
            pytec3d::logging::error_all_ranks("[rank {}/{}] {}: {}",
                                              pytec3d::logging::rank(),
                                              pytec3d::logging::size(),
                                              context,
                                              ex.what());
        }
        catch (...)
        {
        }

        std::cerr << "[pytec3d] " << context << ": " << ex.what() << std::endl;
        throw;
    }

    [[noreturn]] void log_unknown_and_rethrow(const char *context)
    {
        try
        {
            pytec3d::logging::error_all_ranks("[rank {}/{}] {}: unknown exception",
                                              pytec3d::logging::rank(),
                                              pytec3d::logging::size(),
                                              context);
        }
        catch (...)
        {
        }

        std::cerr << "[pytec3d] " << context << ": unknown exception" << std::endl;
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
                pytec3d::logging::error_and_abort("pyTEC3DSolver runtime is permanently closed in this process.");

            if (state_ == RuntimeState::Uninitialized)
            {
                try
                {
                    static const char *default_argv[] = {"pytec3d_solver"};
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
                    pytec3d::logging::initialize(MPI_COMM_WORLD, "pytec3d");
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
                pytec3d::logging::error_and_abort("Cannot shutdown pyTEC3DSolver runtime while solver instances are still active.");

            if (logging_initialized_)
            {
                pytec3d::logging::shutdown();
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

#define pyTEC3D_STR_IMPL(x) #x
#define pyTEC3D_STR(x) pyTEC3D_STR_IMPL(x)

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

    void print_solver_banner_once()
    {
        if (!RuntimeManager::instance().consume_banner_needed())
            return;

        pytec3d::logging::table_root(pytec3d::logging::format_solver_banner(
            pytec3d_version_string(),
            build_mode_string(),
            std::string(__DATE__) + " " + __TIME__));
    }

} // namespace

void pyTEC3DSolver::ensure_libmesh_initialized()
{
    libmesh_init = RuntimeManager::instance().acquire(argc, argv, num_threads_per_mpi_process_);
    print_solver_banner_once();
}

pyTEC3DSolver::pyTEC3DSolver(int num_threads_per_mpi_process)
{
    if (num_threads_per_mpi_process < 1)
        pytec3d::logging::error_and_abort("num_threads_per_mpi_process must be >= 1.");

    num_threads_per_mpi_process_ = num_threads_per_mpi_process;
    ensure_libmesh_initialized();
    closed_ = false;
    pytec3d::logging::info("Solver created successfully");
}

pyTEC3DSolver::~pyTEC3DSolver()
{
    close();
}

void pyTEC3DSolver::close()
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
        log_and_rethrow("pyTEC3DSolver::close failed", ex);
    }
    catch (...)
    {
        log_unknown_and_rethrow("pyTEC3DSolver::close failed");
    }
}

void pyTEC3DSolver::shutdown_runtime()
{
    try
    {
        RuntimeManager::instance().shutdown_runtime();
    }
    catch (const std::exception &ex)
    {
        log_and_rethrow("pyTEC3DSolver::shutdown_runtime failed", ex);
    }
    catch (...)
    {
        log_unknown_and_rethrow("pyTEC3DSolver::shutdown_runtime failed");
    }
}

void pyTEC3DSolver::load_model(pybind11::object pytec3d_object, pybind11::object solver_settings)
{
    try
    {
        if (closed_)
            pytec3d::logging::error_and_abort("Solver instance has been closed and cannot be reused.");

        if (!model_)
            model_ = std::make_shared<pyTEC3DModel>();

        model_->load_into_solver(libmesh_init, std::move(pytec3d_object), std::move(solver_settings));
        loaded = true;
    }
    catch (const std::exception &ex)
    {
        log_and_rethrow("pyTEC3DSolver::load_model failed", ex);
    }
    catch (...)
    {
        log_unknown_and_rethrow("pyTEC3DSolver::load_model failed");
    }
}

void pyTEC3DSolver::solve(const pybind11::dict &time_step)
{
    try
    {
        if (!loaded || !model_ || !model_->is_loaded())
            pytec3d::logging::error_and_abort("Solver has no loaded model. Call load_model() first.");

        model_->solve(time_step);
    }
    catch (const std::exception &ex)
    {
        log_and_rethrow("pyTEC3DSolver::solve failed", ex);
    }
    catch (...)
    {
        log_unknown_and_rethrow("pyTEC3DSolver::solve failed");
    }
}

void pyTEC3DSolver::write_system_to_vtk_grid(const std::string &system_name)
{
    try
    {
        if (!loaded || !model_ || !model_->is_loaded())
            pytec3d::logging::error_and_abort("Solver is not fully initialized. Call load_model() first.");

        model_->write_system_to_vtk_grid(system_name);
    }
    catch (const std::exception &ex)
    {
        log_and_rethrow("pyTEC3DSolver::write_system_to_vtk_grid failed", ex);
    }
    catch (...)
    {
        log_unknown_and_rethrow("pyTEC3DSolver::write_system_to_vtk_grid failed");
    }
}

pybind11::dict pyTEC3DSolver::registries()
{
    namespace reg = pytec3d::registry;
    pybind11::dict out;

    pybind11::list domains;
    for (const auto name : reg::kDomainNames)
        domains.append(std::string(name));
    out["domains"] = domains;

    pybind11::dict fields;
    for (std::size_t i = 0; i < reg::kNumFieldId; ++i)
        fields[pybind11::str(std::string(reg::kFieldIdNames[i]))] = reg::kFieldIsVector[i] ? "vector" : "scalar";
    out["fields"] = fields;

    pybind11::dict element_types;
    for (std::size_t i = 0; i < reg::kNumElementTypeId; ++i)
        element_types[pybind11::str(std::string(reg::kElementTypeIdNames[i]))] = reg::kElementTypeDimensions[i];
    out["element_types"] = element_types;

    const auto types_of = [](const auto &names, const auto &types)
    {
        pybind11::dict result;
        for (std::size_t i = 0; i < names.size(); ++i)
            result[pybind11::str(std::string(names[i]))] = pytec3d::type_name(types[i]);
        return result;
    };
    out["material_properties"] = types_of(reg::kMaterialPropertyNames, reg::kMaterialPropertyTypes);
    out["real_constants"] = types_of(reg::kRealConstantNames, reg::kRealConstantTypes);

    // law -> {"source": "materials" | "constants" | "element_types", "curve": [section, x, y] or None}
    pybind11::dict laws;
    for (std::size_t i = 0; i < reg::kNumLawId; ++i)
    {
        const reg::LawInfo &info = reg::kLaws[i];
        pybind11::dict law;
        law["source"] = info.source == reg::LawSource::Material   ? "materials"
                        : info.source == reg::LawSource::Constant ? "constants"
                                                                  : "element_types";
        if (info.curve == reg::CurveId::COUNT)
            law["curve"] = pybind11::none();
        else
        {
            const reg::CurveSource &curve = reg::kCurveSources[static_cast<std::size_t>(info.curve)];
            law["curve"] = pybind11::make_tuple(std::string(curve.section), std::string(curve.x), std::string(curve.y));
        }
        laws[pybind11::str(std::string(reg::kLawIdNames[i]))] = law;
    }
    out["laws"] = laws;
    return out;
}

bool pyTEC3DSolver::is_loaded() const
{
    return loaded && model_ && model_->is_loaded();
}
