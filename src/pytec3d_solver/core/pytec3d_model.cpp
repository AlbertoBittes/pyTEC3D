/**
 * @file pytec3d_model.cpp
 * @brief pyTEC3DModel: mesh import, equation systems setup, solve and output to the grid.
 */
#include "pytec3d_model.h"

#include "augment_sparsity_on_nodes.h"
#include "logging.h"
#include "mesh_quality.h"
#include "model_element_types.h"
#include "pytec3d_equation_systems.h"
#include "pytec3d_system.h"
#include "report_tables.h"
#include "vtkio.h"

#include <libmesh/distributed_mesh.h>
#include <libmesh/elem.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_tools.h>

pyTEC3DModel::pyTEC3DModel() = default;

pyTEC3DModel::~pyTEC3DModel()
{
    clear();
}

void pyTEC3DModel::load(const libMesh::Parallel::Communicator &comm,
                        vtkUnstructuredGrid *grid,
                        const std::string &model_json,
                        const SolverSettings &solver_settings,
                        const CoordinateRotations &rotations)
{
    if (!grid)
        pytec3d::logging::error_and_abort("pyTEC3DModel::load: no grid.");

    clear();

    pytec3d::logging::info("Loading the model: {} points, {} cells.", grid->GetNumberOfPoints(), grid->GetNumberOfCells());
    grid_ = grid;

    mesh_ = std::make_unique<libMesh::DistributedMesh>(comm);
    {
        pytec3d::logging::ScopedTimer timer("Mesh imported from the VTK grid");
        VTKIO vtkio(*mesh_, grid_);
        vtkio.read(grid_);
        // Every rank holds the whole mesh until it is distributed: the count is the same on all of them.
        for (const libMesh::Elem *elem : mesh_->element_ptr_range())
            flipped_elements_ += elem->is_flipped() ? 1 : 0;
        libMesh::MeshTools::Modification::orient_elements(*mesh_);
    }
    if (flipped_elements_ > 0)
        pytec3d::logging::info("{} element(s) with inverted node ordering reoriented.", flipped_elements_);

    equation_systems_ = std::make_unique<pyTEC3DEquationSystems>(*mesh_, *this);
    equation_systems_->load_model_data(model_json, solver_settings);
    pytec3d::logging::info("Model data: {}", equation_systems_->get_model_data().summary());
    equation_systems_->validate_model_data();
    equation_systems_->configure_systems();
    pytec3d::logging::info("DOF system: {}.", equation_systems_->uses_nonlinear_system()
                                                  ? "nonlinear (cells with nonlinear laws), Newton"
                                                  : "linear");
    equation_systems_->initialize_systems();

    auto &main_system = equation_systems_->get_system<libMesh::System>(
        pyTEC3DSystemBase::system_name_from_type(pyTEC3DSystemBase::DOF_SOLUTION_SYSTEM));
    coupling_functor_ = std::make_unique<AugmentSparsityOnNodes>(*mesh_, grid_);
    main_system.get_dof_map().add_coupling_functor(*coupling_functor_);

    {
        pytec3d::logging::ScopedTimer timer("Mesh prepared and distributed");
        mesh_->prepare_for_use(true);
    }
    pytec3d::logging::table_root(pytec3d::report::format_mesh_summary_table(*mesh_));
    {
        pytec3d::logging::ScopedTimer timer("Mesh quality checked");
        const auto quality = pytec3d::compute_mesh_quality(*mesh_);
        pytec3d::logging::table_root(pytec3d::report::format_mesh_quality_table(quality, flipped_elements_));
        pytec3d::check_mesh_quality(quality);
    }

    equation_systems_->build_rotation_caches(rotations);

    {
        pytec3d::logging::ScopedTimer timer("Degrees of freedom distributed");
        equation_systems_->initialize_physics();
    }
    pytec3d::logging::table_root(pytec3d::report::format_model_table(*equation_systems_));

    loaded_ = true;
}

void pyTEC3DModel::solve(const TimeStepInput &time_step)
{
    require_loaded();

    equation_systems_->solve_equation_systems(time_step);

    {
        pytec3d::logging::ScopedTimer timer("Solution written to the VTK grid",
                                            &equation_systems_->step_report().seconds_output);
        for (const auto type : {pyTEC3DSystemBase::DOF_SOLUTION_SYSTEM,
                                pyTEC3DSystemBase::NODAL_SOLUTION,
                                pyTEC3DSystemBase::ELEMENT_SOLUTION})
            write_system_to_vtk_grid(pyTEC3DSystemBase::system_name_from_type(type));
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

void pyTEC3DModel::write_system_to_vtk_grid(const std::string &system_name)
{
    require_loaded();

    VTKIO vtkio(*mesh_, equation_systems_.get(), grid_);
    vtkio.write_system_grid(system_name);
}

void pyTEC3DModel::clear()
{
    if (equation_systems_)
        equation_systems_->clear();
    equation_systems_.reset();
    coupling_functor_.reset();

    if (mesh_)
        mesh_->clear();
    mesh_.reset();

    grid_ = nullptr;
    flipped_elements_ = 0;
    loaded_ = false;
}

bool pyTEC3DModel::is_loaded() const
{
    return loaded_;
}

vtkUnstructuredGrid *pyTEC3DModel::get_active_vtk_grid() const
{
    return grid_;
}

const AugmentSparsityOnNodes *pyTEC3DModel::get_coupling_functor() const
{
    if (!coupling_functor_)
        pytec3d::logging::error_and_abort("Solver coupling functor is not initialized.");
    return coupling_functor_.get();
}

void pyTEC3DModel::require_loaded() const
{
    if (!loaded_)
        pytec3d::logging::error_and_abort("Solver model is not fully initialized. Call load_model() first.");
}
