/**
 * @file python_bindings.cpp
 * @brief The pybind11 module: converts NumPy arrays and dicts into the solver's C++ types and exposes the session.
 *
 * This is the only file of the solver that depends on Python. It turns
 *  - the model input dict (built by pytec3d.solver_input) into a pytec3d::ModelInput: every NumPy array becomes a
 *    typed ArrayView of its own memory, never a copy: it must have the dtype of its role (float64 for real numbers,
 *    int32 for ids, int64 for mesh indices, uint8 for cell types), be C-contiguous and in the native byte order,
 *    else it is refused; the model data dict becomes a nlohmann::json;
 *  - the time_step and solver_settings dicts into TimeStepInput and SolverSettings.
 * The Session keeps a reference to every model array it views, so that they stay alive while the solver uses them;
 * the mesh and the rotations are only read by load_model(): their references are dropped when it returns, so that
 * arrays computed for the load (the rotations) are freed at once.
 *
 * Errors: in a serial run, C++ exceptions reach Python (pybind11 translates them). On several MPI ranks an error
 * on some ranks only would leave the others waiting in a collective call, so it is logged and MPI_Abort()ed.
 */
#include "logging.h"
#include "model_input.h"
#include "parameter_value.h"
#include "pytec3d_session.h"
#include "registries.h"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <mpi.h>

#include <cstdint>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>

#define pyTEC3D_STR_IMPL(x) #x
#define pyTEC3D_STR(x) pyTEC3D_STR_IMPL(x)

namespace py = pybind11;

namespace
{
    /// Version of the interface between the Python package and this module (the model input dict, the Session
    /// methods); pytec3d checks it when it imports the module.
    constexpr int kInterfaceVersion = 3;

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

    std::string type_name(const py::handle &object)
    {
        return py::str(py::type::handle_of(object).attr("__name__")).cast<std::string>();
    }

    /** @brief NumPy kind and name of the element type @p T of a view. */
    template <typename T>
    constexpr std::pair<char, const char *> numpy_type()
    {
        if constexpr (std::is_same_v<T, double>)
            return {'f', "float64"};
        else if constexpr (std::is_same_v<T, std::int32_t>)
            return {'i', "int32"};
        else if constexpr (std::is_same_v<T, std::int64_t>)
            return {'i', "int64"};
        else
        {
            static_assert(std::is_same_v<T, std::uint8_t>, "unsupported view type");
            return {'u', "uint8"};
        }
    }

    /**
     * @brief A view of the memory of the NumPy array @p object (first axis: tuples, the others: components), which is
     * appended to @p refs to keep it alive. The array must have exactly the dtype of @p T (no conversion: a converted
     * array would be a copy), and be writable for a non-const @p T.
     */
    template <typename T>
    pytec3d::ArrayView<T> view_of(const py::handle &object, const std::string &what, py::list &refs)
    {
        using Value = std::remove_const_t<T>;
        constexpr auto expected = numpy_type<Value>();
        if (!py::isinstance<py::array>(object))
            throw py::type_error(what + " must be a NumPy array, got " + type_name(object) + ".");
        const auto array = py::reinterpret_borrow<py::array>(object);
        if (array.dtype().kind() != expected.first || array.dtype().itemsize() != static_cast<py::ssize_t>(sizeof(Value)))
            throw py::type_error(what + " must be " + expected.second + ", got " + py::str(array.dtype()).cast<std::string>() +
                                 ": the solver views the array without converting it.");
        if (array.ndim() < 1)
            throw py::value_error(what + " must have at least one dimension.");
        if (!(array.flags() & py::array::c_style))
            throw py::value_error(what + " is not C-contiguous: the solver works on the array itself, not on a copy, "
                                         "and needs it contiguous.");
        if (!array.dtype().attr("isnative").cast<bool>())
            throw py::value_error(what + " is not in the native byte order.");
        if constexpr (!std::is_const_v<T>)
            if (!array.writeable())
                throw py::value_error(what + " is read-only: the solver writes into it.");

        std::size_t n_components = 1;
        for (py::ssize_t axis = 1; axis < array.ndim(); ++axis)
            n_components *= static_cast<std::size_t>(array.shape(axis));
        refs.append(array);
        return pytec3d::ArrayView<T>(static_cast<T *>(const_cast<void *>(array.data())),
                                     static_cast<std::size_t>(array.shape(0)), n_components);
    }

    /** @brief A view of @p object, or an empty view when it is None. */
    template <typename T>
    pytec3d::ArrayView<T> optional_view_of(const py::handle &object, const std::string &what, py::list &refs)
    {
        return object.is_none() ? pytec3d::ArrayView<T>() : view_of<T>(object, what, refs);
    }

    py::dict dict_item(const py::dict &dict, const char *key, const std::string &what)
    {
        if (!dict.contains(key))
            throw py::key_error(what + " has no '" + key + "'.");
        const py::object value = dict[key];
        if (!py::isinstance<py::dict>(value))
            throw py::type_error(what + "['" + key + "'] must be a dict, got " + type_name(value) + ".");
        return py::reinterpret_borrow<py::dict>(value);
    }

    /** @brief The views of the arrays of the dict {name: array} @p arrays. */
    template <typename T>
    std::map<std::string, pytec3d::ArrayView<T>> views_of(const py::dict &arrays, const std::string &what, py::list &refs)
    {
        std::map<std::string, pytec3d::ArrayView<T>> out;
        for (const auto &[key, value] : arrays)
        {
            const std::string name = py::str(key).cast<std::string>();
            out.emplace(name, view_of<T>(value, what + "['" + name + "']", refs));
        }
        return out;
    }

    /**
     * @brief The ModelArrays of the dict {"ids": {"TYPE", "MATERIAL", "CONSTANT", "COORDINATE_SYSTEM"} (int32),
     * "couplings": {"CP_<k>": int32}, "point_fields": {name: float64}, "cell_fields": {name: float64}}.
     */
    pytec3d::ModelArrays model_arrays_of(const py::dict &arrays, py::list &refs)
    {
        pytec3d::ModelArrays out;
        const py::dict ids = dict_item(arrays, "ids", "the model arrays");
        const auto id = [&](const char *name)
        {
            if (!ids.contains(name))
                throw py::key_error(std::string("the model arrays have no cell id array '") + name + "'.");
            return view_of<const std::int32_t>(ids[name], std::string("ids['") + name + "']", refs);
        };
        out.type = id("TYPE");
        out.material = id("MATERIAL");
        out.constant = id("CONSTANT");
        out.coordinate_system = id("COORDINATE_SYSTEM");
        out.couplings = views_of<const std::int32_t>(dict_item(arrays, "couplings", "the model arrays"), "couplings", refs);
        out.point_fields = views_of<double>(dict_item(arrays, "point_fields", "the model arrays"), "point_fields", refs);
        out.cell_fields = views_of<double>(dict_item(arrays, "cell_fields", "the model arrays"), "cell_fields", refs);
        return out;
    }

    /** @brief JSON of plain Python data (dicts, lists, tuples, numbers, strings, None, NumPy scalars and arrays). */
    nlohmann::json json_of(const py::handle &object, const std::string &path)
    {
        if (object.is_none())
            return nullptr;
        if (py::isinstance<py::bool_>(object))
            return object.cast<bool>();
        if (py::isinstance<py::int_>(object))
            return object.cast<long long>();
        if (py::isinstance<py::float_>(object))
            return object.cast<double>();
        if (py::isinstance<py::str>(object))
            return object.cast<std::string>();
        if (py::isinstance<py::dict>(object))
        {
            nlohmann::json out = nlohmann::json::object();
            for (const auto &[key, value] : py::reinterpret_borrow<py::dict>(object))
            {
                const std::string name = py::str(key).cast<std::string>();
                out[name] = json_of(value, path + "[" + name + "]");
            }
            return out;
        }
        if (py::isinstance<py::array>(object) || py::hasattr(object, "item"))
            return json_of(object.attr("tolist")(), path); // NumPy arrays and scalars
        if (py::isinstance<py::list>(object) || py::isinstance<py::tuple>(object))
        {
            nlohmann::json out = nlohmann::json::array();
            std::size_t i = 0;
            for (const auto &item : object)
                out.push_back(json_of(item, path + "[" + std::to_string(i++) + "]"));
            return out;
        }
        throw py::type_error(path + ": a " + type_name(object) + " cannot be model data.");
    }

    /**
     * @brief The ModelInput of the model input dict: "mesh" {"points" (float64), "offsets", "connectivity" (int64),
     * "cell_types" (uint8)}, "arrays" (model_arrays_of()), "element_rotations" and "node_rotations" (float64 or None),
     * "model_data" (dict). The arrays viewed are appended to @p load_refs (mesh and rotations: only read by
     * load_model()) and to @p array_refs (the model arrays, used until the end of the session).
     */
    pytec3d::ModelInput model_input_of(const py::dict &input, py::list &load_refs, py::list &array_refs)
    {
        pytec3d::ModelInput out;
        const py::dict mesh = dict_item(input, "mesh", "the model input");
        const auto mesh_item = [&](const char *key) -> py::object
        {
            if (!mesh.contains(key))
                throw py::key_error(std::string("the model input 'mesh' has no '") + key + "'.");
            return mesh[key];
        };
        out.mesh.points = view_of<const double>(mesh_item("points"), "mesh['points']", load_refs);
        out.mesh.offsets = view_of<const std::int64_t>(mesh_item("offsets"), "mesh['offsets']", load_refs);
        out.mesh.connectivity = view_of<const std::int64_t>(mesh_item("connectivity"), "mesh['connectivity']", load_refs);
        out.mesh.cell_types = view_of<const std::uint8_t>(mesh_item("cell_types"), "mesh['cell_types']", load_refs);
        out.arrays = model_arrays_of(dict_item(input, "arrays", "the model input"), array_refs);
        const auto optional_item = [&input](const char *key) -> py::object
        { return input.contains(key) ? py::object(input[key]) : py::object(py::none()); };
        out.element_rotations = optional_view_of<const double>(optional_item("element_rotations"), "element_rotations", load_refs);
        out.node_rotations = optional_view_of<const double>(optional_item("node_rotations"), "node_rotations", load_refs);
        if (!input.contains("model_data"))
            throw py::key_error("the model input has no 'model_data'.");
        out.model_data = json_of(input["model_data"], "model_data");
        return out;
    }

    /** @brief The session with references to the arrays it views. */
    struct Session
    {
        explicit Session(int num_threads_per_mpi_process) : session(num_threads_per_mpi_process) {}

        void close()
        {
            session.close();
            array_refs = py::list();
        }

        pyTEC3DSession session;
        py::list array_refs; ///< The model arrays of the loaded model (the mesh and the rotations are not kept).
    };

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
    module.doc() = "pyTEC3D solver module: the C++ solver session, driven with NumPy views of the model arrays.";
#ifdef pyTEC3D_VERSION
    module.attr("__version__") = pyTEC3D_STR(pyTEC3D_VERSION);
#else
    module.attr("__version__") = "dev";
#endif
    module.attr("INTERFACE_VERSION") = kInterfaceVersion;
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

    py::class_<Session>(module, "Session")
        .def(py::init([](int num_threads_per_mpi_process)
                      { return guarded("Session", [&]
                                       { return std::make_unique<Session>(num_threads_per_mpi_process); }); }),
             py::arg("num_threads_per_mpi_process") = 1)
        .def("load_model",
             [](Session &self, const py::dict &model_input)
             {
                 guarded("Session.load_model", [&]
                         {
                             py::list load_refs; // mesh and rotations: copied by load_model(), released here
                             py::list array_refs;
                             const pytec3d::ModelInput input = model_input_of(model_input, load_refs, array_refs);
                             self.session.load_model(input);
                             self.array_refs = array_refs; });
             },
             py::arg("model_input"),
             "Loads the model: typed views (not copies) of the mesh and of the model arrays, the rotations and the model data "
             "(see pytec3d.solver_input). The mesh and the rotations are only read here (the session keeps no "
             "reference to them); the model arrays must keep their memory while the session uses them, and the "
             "session keeps a reference to each.")
        .def("set_arrays",
             [](Session &self, const py::dict &arrays)
             {
                 guarded("Session.set_arrays", [&]
                         {
                             py::list array_refs;
                             self.session.set_arrays(model_arrays_of(arrays, array_refs));
                             self.array_refs = array_refs; });
             },
             py::arg("arrays"),
             "Replaces the model arrays of the loaded model: {'ids', 'couplings', 'point_fields', 'cell_fields'} of "
             "views (see pytec3d.solver_input.model_arrays).")
        .def("solve",
             [](Session &self, const py::dict &time_step, const py::object &solver_settings)
             {
                 guarded("Session.solve", [&]
                         { self.session.solve(time_step_of(time_step), solver_settings_of(solver_settings)); });
             },
             py::arg("time_step"), py::arg("solver_settings") = py::none(),
             "Solves one time step with the solver settings and writes the solution into the output arrays.")
        .def("write_system_output",
             [](Session &self, const std::string &system_name)
             { guarded("Session.write_system_output", [&] { self.session.write_system_output(system_name); }); },
             py::arg("system_name"))
        .def("close", [](Session &self)
             { guarded("Session.close", [&] { self.close(); }); })
        .def("last_step_report",
             [](const Session &self)
             { return guarded("Session.last_step_report", [&] { return step_report_dict(self.session.last_step_report()); }); },
             "What the last solve did: iterations, residuals, factorization, Newton history, law ranges, timings.")
        .def("is_loaded", [](const Session &self) { return self.session.is_loaded(); });
}
