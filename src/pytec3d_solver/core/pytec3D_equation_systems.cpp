/**
 * @file pytec3d_equation_systems.cpp
 * @brief pyTEC3DEquationSystems: system registration, variable initialisation, and solve orchestration.
 */
#include "pytec3d_equation_systems.h"
#include "element_builder.h"
#include "pytec3d_system.h"
#include "pytec3d_model.h"
#include "element_types_factory.h"
#include "vtkio.h"
#include "logging.h"
#include <libmesh/petsc_linear_solver.h>
#include <libmesh/solver_configuration.h>
#include "solver_configuration.h"

#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cmath>
#include <limits>
#include <utility>
#include <array>
#include <unordered_map>

#include <libmesh/elem.h>
#include <libmesh/node.h>

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

namespace py = pybind11;

namespace
{
    BHCurveMaterialsMap parse_bh_curves_from_pytec3d(const py::object &pytec3d_object)
    {
        if (!py::hasattr(pytec3d_object, "bh_curves"))
            return {};

        py::object bh_object = pytec3d_object.attr("bh_curves");
        if (bh_object.is_none())
            return {};

        BHCurveMaterialsMap curves;
        const py::dict materials_dict = py::cast<py::dict>(bh_object);
        for (const auto &material_item : materials_dict)
        {
            const unsigned int material_id = py::cast<unsigned int>(material_item.first);
            const auto curve_array = py::cast<py::array_t<double, py::array::c_style | py::array::forcecast>>(material_item.second);

            if (curve_array.ndim() != 2 || curve_array.shape(1) != 2 || curve_array.shape(0) < 2)
            {
                pytec3d::logging::error_and_abort(
                    "bh_curves[" + std::to_string(material_id) +
                    "]: expected shape (N, 2) with N >= 2 and columns [H, B].");
            }

            const std::size_t n_values = static_cast<std::size_t>(curve_array.shape(0));
            std::vector<double> h_values(n_values);
            std::vector<double> b_values(n_values);

            const auto unchecked = curve_array.unchecked<2>();
            for (ssize_t i = 0; i < curve_array.shape(0); ++i)
            {
                h_values[static_cast<std::size_t>(i)] = unchecked(i, 0);
                b_values[static_cast<std::size_t>(i)] = unchecked(i, 1);
            }

            curves.emplace(material_id, BHCurve(std::move(h_values), std::move(b_values)));
        }

        return curves;
    }
}

void pyTEC3DEquationSystems::load_model_data(py::object pytec3d_object,
                                             const std::shared_ptr<const pytec3d::ElementTypeRegistry> &registry)
{
    element_types_factory_ = std::make_shared<ElementTypesFactory>(pytec3d_object.attr("element_types"), registry);

    bh_curves_ = parse_bh_curves_from_pytec3d(pytec3d_object);

    materials_rows_.clear();
    constants_rows_.clear();
    element_local_to_global_rotations_.clear();
    node_local_to_global_rotations_.clear();

    py::object materials_attribute = pytec3d_object.attr("materials");
    py::object constants_attribute = pytec3d_object.attr("constants");

    for (const py::handle row_obj : materials_attribute.attr("to_dict")("records"))
    {
        const py::dict row = py::reinterpret_borrow<py::dict>(row_obj);
        const libMesh::dof_id_type id = py::cast<libMesh::dof_id_type>(row["MATERIAL_ID"]);
        auto &values = materials_rows_[id];
        values.resize(MAT_NUM_COLS);
        values[MAT_MUR_X] = py::cast<libMesh::Real>(row["MUR_X"]);
        values[MAT_MUR_Y] = py::cast<libMesh::Real>(row["MUR_Y"]);
        values[MAT_MUR_Z] = py::cast<libMesh::Real>(row["MUR_Z"]);
        values[MAT_RHO_X] = py::cast<libMesh::Real>(row["RHO_X"]);
        values[MAT_RHO_Y] = py::cast<libMesh::Real>(row["RHO_Y"]);
        values[MAT_RHO_Z] = py::cast<libMesh::Real>(row["RHO_Z"]);
        values[MAT_M_X] = py::cast<libMesh::Real>(row["M_X"]);
        values[MAT_M_Y] = py::cast<libMesh::Real>(row["M_Y"]);
        values[MAT_M_Z] = py::cast<libMesh::Real>(row["M_Z"]);
    }

    for (const py::handle row_obj : constants_attribute.attr("to_dict")("records"))
    {
        const py::dict row = py::reinterpret_borrow<py::dict>(row_obj);
        const libMesh::dof_id_type id = py::cast<libMesh::dof_id_type>(row["CONSTANTS_ID"]);
        auto &values = constants_rows_[id];
        values.resize(CONST_NUM_COLS);
        values[CONST_THK] = py::cast<libMesh::Real>(row["THK"]);
        values[CONST_CS] = py::cast<libMesh::Real>(row["CS"]);
        values[CONST_TURNS] = py::cast<libMesh::Real>(row["TURNS"]);
        values[CONST_VOL] = py::cast<libMesh::Real>(row["VOL"]);
        values[CONST_RAD] = py::cast<libMesh::Real>(row["RAD"]);
        values[CONST_T_X] = py::cast<libMesh::Real>(row["T_X"]);
        values[CONST_T_Y] = py::cast<libMesh::Real>(row["T_Y"]);
        values[CONST_T_Z] = py::cast<libMesh::Real>(row["T_Z"]);
        values[CONST_R] = py::cast<libMesh::Real>(row["R"]);
        values[CONST_S] = py::cast<libMesh::Real>(row["S"]);
    }
}

void pyTEC3DEquationSystems::build_rotation_caches(pybind11::object pytec3dmodel_object)
{
    element_local_to_global_rotations_.clear();
    node_local_to_global_rotations_.clear();

    auto *model = this->get_parent_object();

    vtkUnstructuredGrid *vtk_grid = model->get_active_vtk_grid();

    vtkCellData *cell_data = vtk_grid->GetCellData();
    vtkPointData *point_data = vtk_grid->GetPointData();
    vtkDataArray *coordinate_system_cell_array =
        cell_data ? vtkDataArray::SafeDownCast(cell_data->GetArray("COORDINATE_SYSTEM")) : nullptr;
    vtkDataArray *angles_array =
        point_data ? vtkDataArray::SafeDownCast(point_data->GetArray("angles")) : nullptr;

    const libMesh::TensorValue<libMesh::Real> identity_rotation(
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0);

    std::unordered_map<libMesh::dof_id_type, std::vector<libMesh::dof_id_type>> elem_ids_by_cs;
    std::unordered_map<libMesh::dof_id_type, std::vector<std::array<libMesh::Real, 3>>> centers_by_cs;
    const libMesh::MeshBase &mesh = this->get_mesh();

    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        const auto elem_id = elem->id();
        libMesh::dof_id_type cs_id = 0;
        if (coordinate_system_cell_array)
            cs_id = static_cast<libMesh::dof_id_type>(coordinate_system_cell_array->GetComponent(static_cast<vtkIdType>(elem->id()), 0));

        if (cs_id == 0)
        {
            element_local_to_global_rotations_[elem_id] = identity_rotation;
            continue;
        }

        const libMesh::Point center = elem->vertex_average();
        elem_ids_by_cs[cs_id].push_back(elem_id);
        centers_by_cs[cs_id].push_back({static_cast<libMesh::Real>(center(0)),
                                        static_cast<libMesh::Real>(center(1)),
                                        static_cast<libMesh::Real>(center(2))});
    }

    py::gil_scoped_acquire gil;

    py::dict coordinate_systems_dict;
    if (py::hasattr(pytec3dmodel_object, "coordinate_systems"))
    {
        py::object coordinate_systems_obj = pytec3dmodel_object.attr("coordinate_systems");
        if (!coordinate_systems_obj.is_none())
            coordinate_systems_dict = py::cast<py::dict>(coordinate_systems_obj);
    }

    for (const auto &entry : centers_by_cs)
    {
        const auto cs_id = entry.first;
        const auto dict_key = py::int_(cs_id);
        if (!coordinate_systems_dict.contains(dict_key))
            pytec3d::logging::error_and_abort("Missing coordinate_systems entry for id=" + std::to_string(cs_id));

        py::object coordinate_system = coordinate_systems_dict[dict_key];
        const auto &centers = entry.second;
        const auto &elem_ids = elem_ids_by_cs[cs_id];
        const py::ssize_t n = static_cast<py::ssize_t>(centers.size());

        py::array_t<double> positions({n, static_cast<py::ssize_t>(3)});
        auto positions_view = positions.mutable_unchecked<2>();
        for (py::ssize_t i = 0; i < n; ++i)
        {
            positions_view(i, 0) = centers[static_cast<std::size_t>(i)][0];
            positions_view(i, 1) = centers[static_cast<std::size_t>(i)][1];
            positions_view(i, 2) = centers[static_cast<std::size_t>(i)][2];
        }

        py::object matrices_obj = coordinate_system.attr("local_to_global_matrices")(positions);
        const auto matrices = py::cast<py::array_t<double, py::array::c_style | py::array::forcecast>>(matrices_obj);
        if (matrices.ndim() != 3 || matrices.shape(0) != n || matrices.shape(1) != 3 || matrices.shape(2) != 3)
            pytec3d::logging::error_and_abort("Coordinate system rotation cache shape mismatch for id=" + std::to_string(cs_id));

        auto matrices_view = matrices.unchecked<3>();
        for (py::ssize_t i = 0; i < n; ++i)
        {
            libMesh::TensorValue<libMesh::Real> rotation(
                static_cast<libMesh::Real>(matrices_view(i, 0, 0)),
                static_cast<libMesh::Real>(matrices_view(i, 0, 1)),
                static_cast<libMesh::Real>(matrices_view(i, 0, 2)),
                static_cast<libMesh::Real>(matrices_view(i, 1, 0)),
                static_cast<libMesh::Real>(matrices_view(i, 1, 1)),
                static_cast<libMesh::Real>(matrices_view(i, 1, 2)),
                static_cast<libMesh::Real>(matrices_view(i, 2, 0)),
                static_cast<libMesh::Real>(matrices_view(i, 2, 1)),
                static_cast<libMesh::Real>(matrices_view(i, 2, 2)));
            element_local_to_global_rotations_[elem_ids[static_cast<std::size_t>(i)]] = rotation;
        }
    }

    if (angles_array)
    {
        const vtkIdType n_tuples = angles_array->GetNumberOfTuples();
        py::array_t<double> angles({static_cast<py::ssize_t>(n_tuples), static_cast<py::ssize_t>(3)});
        auto angles_view = angles.mutable_unchecked<2>();
        for (vtkIdType i = 0; i < n_tuples; ++i)
        {
            angles_view(i, 0) = angles_array->GetComponent(i, 0);
            angles_view(i, 1) = angles_array->GetComponent(i, 1);
            angles_view(i, 2) = angles_array->GetComponent(i, 2);
        }

        py::module_ coordinate_module = py::module_::import("pytec3d.coordinate_system");
        py::object coordinate_type = coordinate_module.attr("CoordinateSystem");
        py::object coeffs_obj = coordinate_type.attr("build_local_to_global_coefficients")(angles);
        const auto coeffs = py::cast<py::array_t<double, py::array::c_style | py::array::forcecast>>(coeffs_obj);
        if (coeffs.ndim() != 2 || coeffs.shape(0) != n_tuples || coeffs.shape(1) != 9)
            pytec3d::logging::error_and_abort("Node rotation cache shape mismatch from CoordinateSystem.build_local_to_global_coefficients.");

        auto coeffs_view = coeffs.unchecked<2>();
        for (const auto *node : mesh.node_ptr_range())
        {
            const auto node_id = static_cast<libMesh::dof_id_type>(node->id());

            libMesh::TensorValue<libMesh::Real> rotation(
                static_cast<libMesh::Real>(coeffs_view(static_cast<vtkIdType>(node->id()), 0)),
                static_cast<libMesh::Real>(coeffs_view(static_cast<vtkIdType>(node->id()), 1)),
                static_cast<libMesh::Real>(coeffs_view(static_cast<vtkIdType>(node->id()), 2)),
                static_cast<libMesh::Real>(coeffs_view(static_cast<vtkIdType>(node->id()), 3)),
                static_cast<libMesh::Real>(coeffs_view(static_cast<vtkIdType>(node->id()), 4)),
                static_cast<libMesh::Real>(coeffs_view(static_cast<vtkIdType>(node->id()), 5)),
                static_cast<libMesh::Real>(coeffs_view(static_cast<vtkIdType>(node->id()), 6)),
                static_cast<libMesh::Real>(coeffs_view(static_cast<vtkIdType>(node->id()), 7)),
                static_cast<libMesh::Real>(coeffs_view(static_cast<vtkIdType>(node->id()), 8)));
            node_local_to_global_rotations_[node_id] = rotation;
        }
    }
    else
    {
        for (const auto *node : mesh.node_ptr_range())
            node_local_to_global_rotations_[static_cast<libMesh::dof_id_type>(node->id())] = identity_rotation;
    }
}

void pyTEC3DEquationSystems::prepare_time_step_solve(const pybind11::dict &time_step_dict)
{
    const libMesh::Real requested_time = pybind11::cast<libMesh::Real>(time_step_dict["time"]);
    const libMesh::Real requested_theta = pybind11::cast<libMesh::Real>(time_step_dict["theta"]);
    const bool requested_transient_effect = pybind11::cast<bool>(time_step_dict["transient_effects"]);

    if (!time_step_)
    {
        time_step_ = std::make_unique<TimeStep>();
        time_step_->time = requested_time;
        time_step_->theta = requested_theta;
        time_step_->transient_effects = requested_transient_effect;
        time_step_->time_step_idx = 0;
        time_step_->dt = requested_time;
        time_step_->updated = true;
        return;
    }

    const libMesh::Real current_time = time_step_->time;
    const libMesh::Real time_scale = std::max(std::abs(requested_time), std::abs(current_time));
    const libMesh::Real time_tol = std::numeric_limits<libMesh::Real>::epsilon() * time_scale;
    const libMesh::Real time_delta = requested_time - current_time;

    if (time_delta > time_tol)
    {
        old_time_step_ = std::move(time_step_);
        time_step_ = std::make_unique<TimeStep>();
        time_step_->time = requested_time;
        time_step_->theta = requested_theta;
        time_step_->transient_effects = requested_transient_effect;
        time_step_->time_step_idx = old_time_step_->time_step_idx + 1;
        time_step_->dt = requested_time - current_time;
        time_step_->updated = true;
        return;
    }

    if (std::abs(time_delta) <= time_tol)
    {
        const libMesh::Real current_theta = time_step_->theta;
        const bool current_transient_effect = time_step_->transient_effects;
        const libMesh::Real theta_tol = std::numeric_limits<libMesh::Real>::epsilon() * std::abs(current_theta);

        if (std::abs(requested_theta - current_theta) > theta_tol ||
            requested_transient_effect != current_transient_effect)
        {
            pytec3d::logging::error_and_abort("pyTEC3DEquationSystems::prepare_time_step_solve: theta/transient_effects changed while time is unchanged.");
        }

        time_step_->time = requested_time;
        time_step_->theta = requested_theta;
        time_step_->transient_effects = requested_transient_effect;
        time_step_->updated = false;
        return;
    }

    pytec3d::logging::error_and_abort("pyTEC3DEquationSystems::prepare_time_step_solve: received time smaller than current time step.");
}

void pyTEC3DEquationSystems::solve_equation_systems(const pybind11::dict &time_step)
{
    pytec3d::logging::info("Preparing time-step solution");
    auto &system = get_dof_solution_system();

    system.init_time_step(time_step);
    system.write_element_matrices();

    this->solve();
}

void pyTEC3DEquationSystems::configure_systems()
{
    use_nonlinear_system_ = this->has_bh_curves();

    for (int i = 0; i < pyTEC3DSystemBase::NUM_SYSTEMS; ++i)
    {
        pyTEC3DSystemBase::SystemType type = static_cast<pyTEC3DSystemBase::SystemType>(i);
        const std::string &sys_name = pyTEC3DSystemBase::system_name_from_type(type);

        if (type == pyTEC3DSystemBase::DOF_SOLUTION_SYSTEM)
        {
            if (use_nonlinear_system_)
                this->add_system<pyTEC3DNonLinearSystem>(sys_name);
            else
                this->add_system<pyTEC3DLinearSystem>(sys_name);
        }
        else if (type == pyTEC3DSystemBase::NODAL_SOLUTION)
            this->add_system<pyTEC3DNodalSolutionSystem>(sys_name);
        else if (type == pyTEC3DSystemBase::ELEMENT_SOLUTION)
            this->add_system<pyTEC3DElementSolutionSystem>(sys_name);
    }
}

void pyTEC3DEquationSystems::post_initialize_systems()
{
    for (int i = 0; i < pyTEC3DSystemBase::NUM_SYSTEMS; ++i)
    {
        const auto type = static_cast<pyTEC3DSystemBase::SystemType>(i);
        const auto &info = pyTEC3DSystemBase::system_name_from_type(type);

        if (type == pyTEC3DSystemBase::DOF_SOLUTION_SYSTEM)
        {
            auto *system = dynamic_cast<pyTEC3DSystemBase *>(&get_dof_solution_system());
            if (!system)
                pytec3d::logging::error_and_abort("DOF solution system is not an pyTEC3DSystemBase.");
            system->post_init();
        }
        else if (type == pyTEC3DSystemBase::NODAL_SOLUTION)
        {
            auto &system = this->get_system<pyTEC3DNodalSolutionSystem>(info);
            system.post_init();
        }
        else if (type == pyTEC3DSystemBase::ELEMENT_SOLUTION)
        {
            auto &system = this->get_system<pyTEC3DElementSolutionSystem>(info);
            system.post_init();
        }
    }
}

void pyTEC3DEquationSystems::initialize_systems()
{
    pytec3d::logging::info("Setting systems variables");

    for (int i = 0; i < pyTEC3DSystemBase::NUM_SYSTEMS; ++i)
    {
        const auto type = static_cast<pyTEC3DSystemBase::SystemType>(i);
        const auto &info = pyTEC3DSystemBase::system_name_from_type(type);

        if (type == pyTEC3DSystemBase::DOF_SOLUTION_SYSTEM)
        {
            auto *system = dynamic_cast<pyTEC3DSystemBase *>(&get_dof_solution_system());
            if (!system)
                pytec3d::logging::error_and_abort("DOF solution system is not an pyTEC3DSystemBase.");
            system->init_variables();
        }
        else if (type == pyTEC3DSystemBase::NODAL_SOLUTION)
        {
            auto &system = this->get_system<pyTEC3DNodalSolutionSystem>(info);
            system.init_variables();
        }
        else if (type == pyTEC3DSystemBase::ELEMENT_SOLUTION)
        {
            auto &system = this->get_system<pyTEC3DElementSolutionSystem>(info);
            system.init_variables();
        }
    }
}

void pyTEC3DEquationSystems::initialize_physics()
{
    this->init();
    post_initialize_systems();
    pytec3d::logging::table_root(pytec3d::logging::format_equation_systems_summary_table(*this));
}

bool pyTEC3DEquationSystems::uses_nonlinear_system() const
{
    return use_nonlinear_system_;
}

pyTEC3DDOFSystemBase &pyTEC3DEquationSystems::get_dof_solution_system()
{
    const auto &name = pyTEC3DSystemBase::SYSTEM_INFO_MAP.at(pyTEC3DSystemBase::DOF_SOLUTION_SYSTEM);
    if (use_nonlinear_system_)
        return this->get_system<pyTEC3DNonLinearSystem>(name);
    return this->get_system<pyTEC3DLinearSystem>(name);
}

const pyTEC3DDOFSystemBase &pyTEC3DEquationSystems::get_dof_solution_system() const
{
    const auto &name = pyTEC3DSystemBase::SYSTEM_INFO_MAP.at(pyTEC3DSystemBase::DOF_SOLUTION_SYSTEM);
    if (use_nonlinear_system_)
        return this->get_system<pyTEC3DNonLinearSystem>(name);
    return this->get_system<pyTEC3DLinearSystem>(name);
}

const TimeStep *pyTEC3DEquationSystems::current_time_step() const
{
    return time_step_.get();
}

const TimeStep *pyTEC3DEquationSystems::previous_time_step() const
{
    return old_time_step_.get();
}

bool pyTEC3DEquationSystems::has_bh_curves() const
{
    return !bh_curves_.empty();
}

const BHCurveMaterialsMap &pyTEC3DEquationSystems::get_bh_curves() const
{
    return bh_curves_;
}

const std::shared_ptr<ElementTypesFactory> &pyTEC3DEquationSystems::get_element_types_factory() const
{
    return element_types_factory_;
}

const pyTEC3DEquationSystems::PropertyRowsMap &pyTEC3DEquationSystems::get_material_rows() const
{
    return materials_rows_;
}

const pyTEC3DEquationSystems::PropertyRowsMap &pyTEC3DEquationSystems::get_constant_rows() const
{
    return constants_rows_;
}

const pyTEC3DEquationSystems::RotationMap &pyTEC3DEquationSystems::get_element_local_to_global_rotations() const
{
    return element_local_to_global_rotations_;
}

const pyTEC3DEquationSystems::RotationMap &pyTEC3DEquationSystems::get_node_local_to_global_rotations() const
{
    return node_local_to_global_rotations_;
}
