/**
 * @file python_bindings.cpp
 * @brief The pybind11 module: converts the Python objects into the solver's C++ types and exposes the session.
 *
 * This is the only file of the solver that depends on Python. It turns
 *  - the Python model into its grid (model.mesh), its data as JSON text (model.data_json()) and the rotations of its
 *    coordinate systems (model.coordinate_system(id), CoordinateSystem.build_local_to_global_coefficients);
 *  - the time_step and solver_settings dicts into TimeStepInput and SolverSettings.
 *
 * Errors: in a serial run, C++ exceptions reach Python (pybind11 translates them). On several MPI ranks an error
 * on some ranks only would leave the others waiting in a collective call, so it is logged and MPI_Abort()ed.
 */
#include "logging.h"
#include "parameter_value.h"
#include "pytec3d_session.h"
#include "registries.h"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <vtkObjectBase.h>
#include <vtkPythonUtil.h>
#include <vtkUnstructuredGrid.h>

#include <mpi.h>

#include <string>
#include <utility>

namespace py = pybind11;

namespace
{
    /** @brief On several ranks, logs the error and aborts all ranks; in a serial run the exception goes on to Python. */
    void abort_if_parallel(const char *context, const char *what)
    {
        if (pytec3d::logging::size() <= 1)
            return;
        pytec3d::logging::error_all_ranks("[rank {}/{}] {}: {}", pytec3d::logging::rank(), pytec3d::logging::size(),
                                          context, what);
        spdlog::default_logger()->flush();
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /** @brief Runs @p f, applying abort_if_parallel() to any exception. */
    template <typename F>
    decltype(auto) guarded(const char *context, F &&f)
    {
        try
        {
            return std::forward<F>(f)();
        }
        catch (const std::exception &ex)
        {
            abort_if_parallel(context, ex.what());
            throw;
        }
        catch (...)
        {
            abort_if_parallel(context, "unknown exception");
            throw;
        }
    }

    using DoubleArray = py::array_t<double, py::array::c_style | py::array::forcecast>;

    vtkUnstructuredGrid *grid_of(const py::object &model)
    {
        const py::object mesh = model.attr("mesh");
        auto *base = vtkPythonUtil::GetPointerFromObject(mesh.ptr(), "vtkUnstructuredGrid");
        auto *grid = vtkUnstructuredGrid::SafeDownCast(base);
        if (!grid)
            throw std::invalid_argument("model.mesh is not a vtkUnstructuredGrid.");
        return grid;
    }

    libMesh::RealTensorValue tensor(const double *v)
    {
        return libMesh::RealTensorValue(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]);
    }

    /** @brief The rotations of the coordinate systems of the Python model; only used while the model is loaded. */
    CoordinateRotations rotations_of(const py::object &model)
    {
        CoordinateRotations rotations;
        rotations.at_points = [model](unsigned int cs_id, const std::vector<libMesh::Point> &points)
        {
            const auto n = static_cast<py::ssize_t>(points.size());
            DoubleArray positions({n, py::ssize_t(3)});
            auto view = positions.mutable_unchecked<2>();
            for (py::ssize_t i = 0; i < n; ++i)
                for (py::ssize_t j = 0; j < 3; ++j)
                    view(i, j) = points[static_cast<std::size_t>(i)](static_cast<unsigned int>(j));

            const py::object coordinate_system = model.attr("coordinate_system")(cs_id);
            const auto matrices = py::cast<DoubleArray>(coordinate_system.attr("local_to_global_matrices")(positions));
            if (matrices.ndim() != 3 || matrices.shape(0) != n || matrices.shape(1) != 3 || matrices.shape(2) != 3)
                throw std::runtime_error("local_to_global_matrices of coordinate system " + std::to_string(cs_id) +
                                         " did not return an (N, 3, 3) array.");
            CoordinateRotations::Rotations out;
            out.reserve(points.size());
            for (py::ssize_t i = 0; i < n; ++i)
                out.push_back(tensor(matrices.data(i, 0, 0)));
            return out;
        };
        rotations.from_angles = [](const std::vector<std::array<libMesh::Real, 3>> &angles)
        {
            const auto n = static_cast<py::ssize_t>(angles.size());
            DoubleArray values({n, py::ssize_t(3)});
            auto view = values.mutable_unchecked<2>();
            for (py::ssize_t i = 0; i < n; ++i)
                for (py::ssize_t j = 0; j < 3; ++j)
                    view(i, j) = angles[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];

            const py::object coordinate_system = py::module_::import("pytec3d.coordinate_system").attr("CoordinateSystem");
            const auto coefficients = py::cast<DoubleArray>(coordinate_system.attr("build_local_to_global_coefficients")(values));
            if (coefficients.ndim() != 2 || coefficients.shape(0) != n || coefficients.shape(1) != 9)
                throw std::runtime_error("CoordinateSystem.build_local_to_global_coefficients did not return an (N, 9) array.");
            CoordinateRotations::Rotations out;
            out.reserve(angles.size());
            for (py::ssize_t i = 0; i < n; ++i)
                out.push_back(tensor(coefficients.data(i, 0)));
            return out;
        };
        return rotations;
    }

    /** @brief The value of @p key in the dict @p group, or @p fallback when absent. */
    template <typename T>
    T setting(const py::dict &group, const char *key, const T &fallback)
    {
        return group.contains(key) ? group[key].cast<T>() : fallback;
    }

    /** @brief The nested solver_settings dict ("linear", "nonlinear", "petsc_options"); missing keys keep the defaults. */
    SolverSettings solver_settings_of(const py::object &settings)
    {
        SolverSettings out;
        if (settings.is_none())
            return out;
        const py::dict all = settings.cast<py::dict>();
        if (all.contains("linear"))
        {
            const py::dict group = all["linear"].cast<py::dict>();
            LinearSolverSettings &l = out.linear;
            l.ksp_type = setting(group, "ksp_type", l.ksp_type);
            l.pc_type = setting(group, "pc_type", l.pc_type);
            l.relative_tolerance = setting(group, "relative_tolerance", l.relative_tolerance);
            l.absolute_tolerance = setting(group, "absolute_tolerance", l.absolute_tolerance);
            if (group.contains("divergence_tolerance") && !group["divergence_tolerance"].is_none())
                l.divergence_tolerance = group["divergence_tolerance"].cast<double>();
            l.max_iterations = setting(group, "max_iterations", l.max_iterations);
            l.initial_guess_nonzero = setting(group, "initial_guess_nonzero", l.initial_guess_nonzero);
            if (group.contains("factor_solver_type"))
                l.factor_solver_type = group["factor_solver_type"].is_none() ? std::string()
                                                                             : group["factor_solver_type"].cast<std::string>();
        }
        if (all.contains("nonlinear"))
        {
            const py::dict group = all["nonlinear"].cast<py::dict>();
            NewtonSettings &n = out.nonlinear;
            n.max_iterations = setting(group, "max_iterations", n.max_iterations);
            n.absolute_residual_tolerance = setting(group, "absolute_residual_tolerance", n.absolute_residual_tolerance);
            n.relative_residual_tolerance = setting(group, "relative_residual_tolerance", n.relative_residual_tolerance);
            n.absolute_step_tolerance = setting(group, "absolute_step_tolerance", n.absolute_step_tolerance);
            n.relative_step_tolerance = setting(group, "relative_step_tolerance", n.relative_step_tolerance);
            n.require_residual_reduction = setting(group, "require_residual_reduction", n.require_residual_reduction);
            n.require_finite_residual = setting(group, "require_finite_residual", n.require_finite_residual);
            n.brent_line_search = setting(group, "brent_line_search", n.brent_line_search);
            n.minimum_step_length = setting(group, "minimum_step_length", n.minimum_step_length);
            n.linear_tolerance_multiplier = setting(group, "linear_tolerance_multiplier", n.linear_tolerance_multiplier);
            n.continue_after_max_iterations = setting(group, "continue_after_max_iterations", n.continue_after_max_iterations);
            n.continue_after_backtrack_failure =
                setting(group, "continue_after_backtrack_failure", n.continue_after_backtrack_failure);
            n.verbose = setting(group, "verbose", n.verbose);
        }
        if (all.contains("petsc_options"))
        {
            out.petsc_options.clear();
            for (const auto &[name, value] : all["petsc_options"].cast<py::dict>())
                if (!value.is_none())
                    out.petsc_options.emplace_back(name.cast<std::string>(), py::str(value).cast<std::string>());
        }
        return out;
    }

    /** @brief time, theta and transient_effects are required; beta and gamma default to average acceleration. */
    TimeStepInput time_step_of(const py::dict &time_step)
    {
        TimeStepInput out;
        out.time = time_step["time"].cast<libMesh::Real>();
        out.theta = time_step["theta"].cast<libMesh::Real>();
        out.transient_effects = time_step["transient_effects"].cast<bool>();
        if (time_step.contains("beta"))
            out.beta = time_step["beta"].cast<libMesh::Real>();
        if (time_step.contains("gamma"))
            out.gamma = time_step["gamma"].cast<libMesh::Real>();
        return out;
    }

    /** @brief The step report as a dict (the Newton history and law ranges as lists of dicts). */
    py::dict step_report_dict(const StepReport &r)
    {
        py::dict out;
        out["time"] = r.time;
        out["dt"] = r.dt;
        out["index"] = r.index;
        out["transient"] = r.transient;
        out["advanced"] = r.advanced;
        out["elements"] = r.elements;
        out["total_dofs"] = r.total_dofs;
        out["unconstrained_dofs"] = r.unconstrained_dofs;
        out["element_matrices_computed"] = r.element_matrices_computed;
        out["system_matrix_updated"] = r.system_matrix_updated;
        out["nonlinear"] = r.nonlinear;
        out["linear_iterations"] = r.linear_iterations;
        out["linear_residual"] = r.linear_residual;
        out["linear_reason"] = r.linear_reason;
        out["true_residual"] = r.true_residual < 0 ? py::object(py::none()) : py::object(py::float_(r.true_residual));
        out["newton_iterations"] = r.newton_iterations;
        out["newton_result"] = r.newton_result;

        const FactorizationReport &f = r.factorization;
        if (f.available)
        {
            py::dict factorization;
            factorization["status"] = f.status;
            factorization["status_detail"] = f.status_detail;
            factorization["flops"] = f.flops;
            factorization["memory_max_mb"] = f.memory_max_mb;
            factorization["memory_total_mb"] = f.memory_total_mb;
            factorization["matrix_entries"] = f.matrix_entries;
            factorization["factor_entries"] = f.factor_entries;
            factorization["null_pivots"] = f.null_pivots < 0 ? py::object(py::none()) : py::object(py::int_(f.null_pivots));
            factorization["delayed_pivots"] = f.delayed_pivots;
            out["factorization"] = factorization;
        }
        else
            out["factorization"] = py::none();

        py::list history;
        for (const NewtonIteration &it : r.newton_history)
        {
            py::dict row;
            row["residual"] = it.residual;
            row["step"] = it.step < 0 ? py::object(py::none()) : py::object(py::float_(it.step));
            row["relative_step"] = it.relative_step < 0 ? py::object(py::none()) : py::object(py::float_(it.relative_step));
            row["step_length"] = it.step_length < 0 ? py::object(py::none()) : py::object(py::float_(it.step_length));
            history.append(row);
        }
        out["newton_history"] = history;

        py::list ranges;
        for (const LawRange &range : r.law_ranges)
        {
            py::dict row;
            row["law"] = range.law;
            row["input"] = range.input;
            row["source"] = range.source;
            row["id"] = range.id;
            row["max_input"] = range.max_input;
            row["curve_max"] = range.curve_max;
            row["extrapolated"] = range.extrapolated;
            ranges.append(row);
        }
        out["law_ranges"] = ranges;

        py::dict seconds;
        seconds["element_matrices"] = r.seconds_element_matrices;
        seconds["solve"] = r.seconds_solve;
        seconds["nodal_fields"] = r.seconds_nodal_fields;
        seconds["element_fields"] = r.seconds_element_fields;
        seconds["output"] = r.seconds_output;
        out["seconds"] = seconds;
        return out;
    }

    py::dict registries()
    {
        namespace reg = pytec3d::registry;
        py::dict out;

        py::list domains;
        for (const auto name : reg::kDomainNames)
            domains.append(std::string(name));
        out["domains"] = domains;

        py::dict fields;
        for (std::size_t i = 0; i < reg::kNumFieldId; ++i)
            fields[py::str(std::string(reg::kFieldIdNames[i]))] = reg::kFieldIsVector[i] ? "vector" : "scalar";
        out["fields"] = fields;

        py::dict element_types;
        for (std::size_t i = 0; i < reg::kNumElementTypeId; ++i)
            element_types[py::str(std::string(reg::kElementTypeIdNames[i]))] = reg::kElementTypeDimensions[i];
        out["element_types"] = element_types;

        const auto types_of = [](const auto &names, const auto &types)
        {
            py::dict result;
            for (std::size_t i = 0; i < names.size(); ++i)
                result[py::str(std::string(names[i]))] = pytec3d::type_name(types[i]);
            return result;
        };
        out["material_properties"] = types_of(reg::kMaterialPropertyNames, reg::kMaterialPropertyTypes);
        out["real_constants"] = types_of(reg::kRealConstantNames, reg::kRealConstantTypes);

        // law -> {"source": "materials" | "constants" | "element_types" (what the ids of its "nonlinear_laws" entries are),
        //         "curve": (x, y) array names of its entries, or None}
        py::dict laws;
        for (std::size_t i = 0; i < reg::kNumLawId; ++i)
        {
            const reg::LawInfo &info = reg::kLaws[i];
            py::dict law;
            law["source"] = info.source == reg::LawSource::Material   ? "materials"
                            : info.source == reg::LawSource::Constant ? "constants"
                                                                      : "element_types";
            if (info.curve == reg::CurveId::COUNT)
                law["curve"] = py::none();
            else
            {
                const reg::CurveSource &curve = reg::kCurveSources[static_cast<std::size_t>(info.curve)];
                law["curve"] = py::make_tuple(std::string(curve.x), std::string(curve.y));
            }
            laws[py::str(std::string(reg::kLawIdNames[i]))] = law;
        }
        out["laws"] = laws;
        return out;
    }
} // namespace

PYBIND11_MODULE(pytec3d_solver, module)
{
    module.doc() = "pyTEC3D solver module (Python-owned model, C++ execution core).";
    module.def("shutdown_runtime", []
               { guarded("shutdown_runtime", [] { pyTEC3DSession::shutdown_runtime(); }); },
               "Shuts down the process-wide libMesh / MPI runtime; call once, after closing every session.");
    module.def("set_log_file", [](const std::string &path) { pytec3d::logging::set_log_file(path); },
               py::arg("path"),
               "Also write the solver log (rank 0) to the file, appending; can be called before the first session. "
               "An empty path stops writing to a file.");
    module.def("message_counts", [] { return pytec3d::logging::message_counts(); },
               "(warnings, errors) logged by the solver on this rank.");
    module.def("registries", &registries,
               "Domains, fields, element types, material properties and real constants compiled into the module.");

    py::class_<pyTEC3DSession>(module, "Session")
        .def(py::init([](int num_threads_per_mpi_process)
                      { return guarded("Session", [&]
                                       { return std::make_unique<pyTEC3DSession>(num_threads_per_mpi_process); }); }),
             py::arg("num_threads_per_mpi_process") = 1)
        .def("load_model",
             [](pyTEC3DSession &session, const py::object &model, const py::object &solver_settings)
             {
                 guarded("Session.load_model", [&]
                         { session.load_model(grid_of(model), model.attr("data_json")().cast<std::string>(),
                                              solver_settings_of(solver_settings), rotations_of(model)); });
             },
             py::arg("model"), py::arg("solver_settings") = py::none(),
             "Loads the model: keeps a reference to model.mesh, copies model.data and the solver settings.")
        .def("solve",
             [](pyTEC3DSession &session, const py::dict &time_step)
             { guarded("Session.solve", [&] { session.solve(time_step_of(time_step)); }); },
             py::arg("time_step"))
        .def("write_system_to_vtk_grid",
             [](pyTEC3DSession &session, const std::string &system_name)
             { guarded("Session.write_system_to_vtk_grid", [&] { session.write_system_to_vtk_grid(system_name); }); },
             py::arg("system_name"))
        .def("close", [](pyTEC3DSession &session)
             { guarded("Session.close", [&] { session.close(); }); })
        .def("last_step_report",
             [](const pyTEC3DSession &session)
             { return guarded("Session.last_step_report", [&] { return step_report_dict(session.last_step_report()); }); },
             "What the last solve did: iterations, residuals, factorization, Newton history, law ranges, timings.")
        .def("is_loaded", &pyTEC3DSession::is_loaded);
}
