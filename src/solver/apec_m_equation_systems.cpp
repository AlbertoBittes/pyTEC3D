/**
 * @file apec_m_equation_systems.cpp
 * @brief ApecMEquationSystems: system registration, variable initialisation, and solve orchestration.
 */
#include "apec_m_solver/apec_m_equation_systems.h"
#include "apec_m_solver/element_builder.h"
#include "apec_m_solver/apec_m_system.h"
#include "apec_m_solver/apec_m_model.h"
#include "apec_m_coordinate/coordinate_system.h"
#include "apec_m_solver/element_types_factory.h"
#include "apec_m_solver/vtkio_updated.h"
#include "apec_m_solver/apec_m_logging.h"
#include <libmesh/petsc_linear_solver.h>
#include <libmesh/solver_configuration.h>
#include "apec_m_solver/solver_configuration.h"

#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cmath>
#include <limits>
#include <utility>

namespace py = pybind11;

namespace
{
    BHCurveMaterialsMap parse_bh_curves_from_apec(const py::object &apec_object)
    {
        if (!py::hasattr(apec_object, "bh_curves"))
            return {};

        py::object bh_object = apec_object.attr("bh_curves");
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
                apec::logging::error_and_abort(
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

void ApecMEquationSystems::load_model_data(py::object apec_object)
{
    element_types_factory_ = std::make_shared<ElementTypesFactory>(apec_object.attr("element_types"));

    bh_curves_ = parse_bh_curves_from_apec(apec_object);

    materials_rows_.clear();
    constants_rows_.clear();
    coordinate_systems_.clear();

    py::object materials_attribute = apec_object.attr("materials");
    py::object constants_attribute = apec_object.attr("constants");

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

    const py::object coordinate_systems_object = apec_object.attr("coordinate_systems");
    if (!coordinate_systems_object.is_none())
    {
        const py::dict coordinate_systems_dict = py::cast<py::dict>(coordinate_systems_object);
        for (const auto &item : coordinate_systems_dict)
        {
            const libMesh::dof_id_type id = py::cast<libMesh::dof_id_type>(item.first);
            coordinate_systems_.emplace(
                id,
                py::cast<std::shared_ptr<apec_m::CoordinateSystem>>(item.second));
        }
    }
}

void ApecMEquationSystems::prepare_time_step_solve(const pybind11::dict &time_step_dict)
{
    const libMesh::Real requested_time = pybind11::cast<libMesh::Real>(time_step_dict["time"]);
    const libMesh::Real requested_theta = pybind11::cast<libMesh::Real>(time_step_dict["theta"]);
    const bool requested_transient_effect = pybind11::cast<bool>(time_step_dict["transient_effects"]);

    if (!time_step_)
    {
        time_step_ = std::make_unique<ApecMTimeStep>();
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
        time_step_ = std::make_unique<ApecMTimeStep>();
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
            apec::logging::error_and_abort("ApecMEquationSystems::prepare_time_step_solve: theta/transient_effects changed while time is unchanged.");
        }

        time_step_->time = requested_time;
        time_step_->theta = requested_theta;
        time_step_->transient_effects = requested_transient_effect;
        time_step_->updated = false;
        return;
    }

    apec::logging::error_and_abort("ApecMEquationSystems::prepare_time_step_solve: received time smaller than current time step.");
}

void ApecMEquationSystems::solve_equation_systems(const pybind11::dict &time_step)
{
    apec::logging::info("Preparing time-step solution");
    auto &system = get_dof_solution_system();

    system.init_time_step(time_step);
    system.write_element_data();

    this->solve();
}

void ApecMEquationSystems::configure_systems()
{
    use_nonlinear_system_ = this->has_bh_curves();

    for (int i = 0; i < ApecMSystemBase::NUM_SYSTEMS; ++i)
    {
        ApecMSystemBase::SystemType type = static_cast<ApecMSystemBase::SystemType>(i);
        const std::string &sys_name = ApecMSystemBase::system_name_from_type(type);

        if (type == ApecMSystemBase::DOF_SOLUTION_SYSTEM)
        {
            if (use_nonlinear_system_)
                this->add_system<ApecMNonLinearSystem>(sys_name);
            else
                this->add_system<ApecMLinearSystem>(sys_name);
        }
        else if (type == ApecMSystemBase::NODAL_SOLUTION)
            this->add_system<ApecMNodalSolutionSystem>(sys_name);
        else if (type == ApecMSystemBase::ELEMENT_SOLUTION)
            this->add_system<ApecMElementSolutionSystem>(sys_name);
    }
}

void ApecMEquationSystems::post_initialize_systems()
{
    for (int i = 0; i < ApecMSystemBase::NUM_SYSTEMS; ++i)
    {
        const auto type = static_cast<ApecMSystemBase::SystemType>(i);
        const auto &info = ApecMSystemBase::system_name_from_type(type);

        if (type == ApecMSystemBase::DOF_SOLUTION_SYSTEM)
        {
            auto *system = dynamic_cast<ApecMSystemBase *>(&get_dof_solution_system());
            if (!system)
                apec::logging::error_and_abort("DOF solution system is not an ApecMSystemBase.");
            system->post_init();
        }
        else if (type == ApecMSystemBase::NODAL_SOLUTION)
        {
            auto &system = this->get_system<ApecMNodalSolutionSystem>(info);
            system.post_init();
        }
        else if (type == ApecMSystemBase::ELEMENT_SOLUTION)
        {
            auto &system = this->get_system<ApecMElementSolutionSystem>(info);
            system.post_init();
        }
    }
}

void ApecMEquationSystems::initialize_systems()
{
    apec::logging::info("Setting systems variables");

    for (int i = 0; i < ApecMSystemBase::NUM_SYSTEMS; ++i)
    {
        const auto type = static_cast<ApecMSystemBase::SystemType>(i);
        const auto &info = ApecMSystemBase::system_name_from_type(type);

        if (type == ApecMSystemBase::DOF_SOLUTION_SYSTEM)
        {
            auto *system = dynamic_cast<ApecMSystemBase *>(&get_dof_solution_system());
            if (!system)
                apec::logging::error_and_abort("DOF solution system is not an ApecMSystemBase.");
            system->init_variables();
        }
        else if (type == ApecMSystemBase::NODAL_SOLUTION)
        {
            auto &system = this->get_system<ApecMNodalSolutionSystem>(info);
            system.init_variables();
        }
        else if (type == ApecMSystemBase::ELEMENT_SOLUTION)
        {
            auto &system = this->get_system<ApecMElementSolutionSystem>(info);
            system.init_variables();
        }
    }
}

void ApecMEquationSystems::initialize_physics()
{
    this->init();
    post_initialize_systems();
    apec::logging::table_root(apec::logging::format_equation_systems_summary_table(*this));
}

bool ApecMEquationSystems::uses_nonlinear_system() const
{
    return use_nonlinear_system_;
}

ApecMDOFSystemBase &ApecMEquationSystems::get_dof_solution_system()
{
    const auto &name = ApecMSystemBase::SYSTEM_INFO_MAP.at(ApecMSystemBase::DOF_SOLUTION_SYSTEM);
    if (use_nonlinear_system_)
        return this->get_system<ApecMNonLinearSystem>(name);
    return this->get_system<ApecMLinearSystem>(name);
}

const ApecMDOFSystemBase &ApecMEquationSystems::get_dof_solution_system() const
{
    const auto &name = ApecMSystemBase::SYSTEM_INFO_MAP.at(ApecMSystemBase::DOF_SOLUTION_SYSTEM);
    if (use_nonlinear_system_)
        return this->get_system<ApecMNonLinearSystem>(name);
    return this->get_system<ApecMLinearSystem>(name);
}

const ApecMTimeStep *ApecMEquationSystems::current_time_step() const
{
    return time_step_.get();
}

const ApecMTimeStep *ApecMEquationSystems::previous_time_step() const
{
    return old_time_step_.get();
}

bool ApecMEquationSystems::has_bh_curves() const
{
    return !bh_curves_.empty();
}

const BHCurveMaterialsMap &ApecMEquationSystems::get_bh_curves() const
{
    return bh_curves_;
}

const std::shared_ptr<ElementTypesFactory> &ApecMEquationSystems::get_element_types_factory() const
{
    return element_types_factory_;
}

const ApecMEquationSystems::PropertyRowsMap &ApecMEquationSystems::get_material_rows() const
{
    return materials_rows_;
}

const ApecMEquationSystems::PropertyRowsMap &ApecMEquationSystems::get_constant_rows() const
{
    return constants_rows_;
}

const ApecMEquationSystems::CoordinateSystemsMap &ApecMEquationSystems::get_coordinate_systems() const
{
    return coordinate_systems_;
}
