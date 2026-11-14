/**
 * @file pytec3d_model.cpp
 * @brief pyTEC3DModel: mesh import, equation systems setup, solve and output to the caller's arrays.
 */
#include "pytec3d_model.h"

#include "array_io.h"
#include "fe_mesh.h"
#include "logging.h"
#include "mesh_quality.h"
#include "model_element_types.h"
#include "pytec3d_equation_systems.h"
#include "report_tables.h"

namespace
{
    /** @brief The systems whose fields are written into the output arrays after each solve. */
    constexpr pyTEC3DEquationSystems::SystemType output_systems[] = {pyTEC3DEquationSystems::DOF_SOLUTION_SYSTEM,
                                                                     pyTEC3DEquationSystems::NODAL_SOLUTION,
                                                                     pyTEC3DEquationSystems::ELEMENT_SOLUTION};
} // namespace

pyTEC3DModel::pyTEC3DModel() = default;

pyTEC3DModel::~pyTEC3DModel()
{
    clear();
}

void pyTEC3DModel::load(MPI_Comm comm, const pytec3d::ModelInput &input)
{
    clear();

    input.validate();
    n_points_ = input.mesh.n_points();
    n_cells_ = input.mesh.n_cells();
    arrays_ = input.arrays;
    pytec3d::logging::info("Loading the model: {} points, {} cells.", n_points_, n_cells_);

    {
        pytec3d::logging::ScopedTimer timer("Mesh imported, partitioned and distributed");
        mesh_ = std::make_unique<pytec3d::fem::FEMesh>(comm, input.mesh, arrays_.type);
    }
    if (mesh_->flipped_elements() > 0)
        pytec3d::logging::info("{} element(s) with inverted node ordering reoriented.", mesh_->flipped_elements());

    equation_systems_ = std::make_unique<pyTEC3DEquationSystems>(*mesh_, *this);
    equation_systems_->load_model_data(input.model_data);
    required_fields_ = pytec3d::required_field_arrays(equation_systems_->element_types());
    pytec3d::check_field_arrays(arrays_, required_fields_, n_points_, n_cells_, true);
    pytec3d::logging::info("Model data: {}", equation_systems_->get_model_data().summary());
    equation_systems_->validate_model_data();

    pytec3d::logging::table_root(pytec3d::report::format_mesh_summary_table(*mesh_));
    {
        pytec3d::logging::ScopedTimer timer("Mesh quality checked");
        const auto quality = pytec3d::compute_mesh_quality(*mesh_);
        pytec3d::logging::table_root(pytec3d::report::format_mesh_quality_table(quality, mesh_->flipped_elements()));
        pytec3d::check_mesh_quality(quality);
    }

    equation_systems_->build_rotation_caches(input.element_rotations, input.node_rotations);
    equation_systems_->initialize_systems();
    pytec3d::logging::table_root(pytec3d::report::format_model_table(*equation_systems_));
    check_outputs();

    loaded_ = true;
}

void pyTEC3DModel::set_arrays(const pytec3d::ModelArrays &arrays)
{
    require_loaded();
    arrays.validate(n_points_, n_cells_);
    pytec3d::check_field_arrays(arrays, required_fields_, n_points_, n_cells_, false);
    arrays_ = arrays;
}

void pyTEC3DModel::check_outputs() const
{
    std::string missing;
    for (const auto type : output_systems)
        for (const std::string &name : equation_systems_->missing_outputs(type, arrays_))
            missing += (missing.empty() ? "" : ", ") + name;
    if (!missing.empty())
        pytec3d::logging::error_and_abort("The solution systems write fields the element types do not declare: " + missing + ".");
}

void pyTEC3DModel::solve(const TimeStepInput &time_step, const SolverSettings &settings)
{
    require_loaded();
    equation_systems_->solve_equation_systems(time_step, settings);
    {
        pytec3d::logging::ScopedTimer timer("Solution written to the output arrays",
                                            &equation_systems_->step_report().seconds_output);
        for (const auto type : output_systems)
            equation_systems_->write_output(type, arrays_);
    }
    const StepReport &report = equation_systems_->step_report();
    if (!report.newton_history.empty())
        pytec3d::logging::table_root(pytec3d::report::format_newton_history_table(report));
    pytec3d::logging::table_root(pytec3d::report::format_step_report_table(report));
}

const StepReport &pyTEC3DModel::last_step_report() const
{
    require_loaded();
    return equation_systems_->step_report();
}

void pyTEC3DModel::write_system_output(const std::string &system_name)
{
    require_loaded();
    for (int t = 0; t < pyTEC3DEquationSystems::NUM_SYSTEMS; ++t)
    {
        const auto type = static_cast<pyTEC3DEquationSystems::SystemType>(t);
        if (pyTEC3DEquationSystems::system_name(type) == system_name)
        {
            equation_systems_->write_output(type, arrays_);
            return;
        }
    }
    pytec3d::logging::error_and_abort("Solution output: system '" + system_name + "' not found.");
}

void pyTEC3DModel::clear()
{
    equation_systems_.reset();
    mesh_.reset();
    arrays_ = {};
    required_fields_.clear();
    n_points_ = 0;
    n_cells_ = 0;
    loaded_ = false;
}

bool pyTEC3DModel::is_loaded() const
{
    return loaded_;
}

void pyTEC3DModel::require_loaded() const
{
    if (!loaded_)
        pytec3d::logging::error_and_abort("Solver model is not fully initialized. Call load_model() first.");
}
