/**
 * @file pytec3d_equation_systems.h
 * @brief pyTEC3DEquationSystems: the systems of one model on its mesh, and the solve of a time step.
 *
 * Three systems, as before:
 *   - the DOF system (dof_system.h): the unknowns of the element types (fem::DofSpace), linear or nonlinear;
 *   - the nodal fields: L2 projections on the element's nodal shape functions, averaged at the points;
 *   - the element fields: one value per element (quadrature averages).
 * Their variables follow the element types of the model data: one per field label, active on the TYPE ids whose
 * element types have the field, in the order of the variable registry.
 */
#pragma once

#include <memory>
#include <string>
#include <vector>

#include "dof_space.h"
#include "fe_mesh.h"
#include "quadrature.h"
#include "la.h"
#include "model_data.h"
#include "model_input.h"
#include "solver_configuration.h"
#include "time_step.h"

class ModelElementTypes;
class pyTEC3DModel;
class DofSystem;

class pyTEC3DEquationSystems
{
public:
    /** @brief The systems, as the output names them. */
    enum SystemType
    {
        DOF_SOLUTION_SYSTEM = 0,
        NODAL_SOLUTION,
        ELEMENT_SOLUTION,
        NUM_SYSTEMS
    };
    static const std::string &system_name(SystemType type);

    pyTEC3DEquationSystems(pytec3d::fem::FEMesh &mesh, pyTEC3DModel &model);
    ~pyTEC3DEquationSystems();

    pyTEC3DModel &model() const { return model_; }
    const pytec3d::fem::FEMesh &mesh() const { return mesh_; }
    const pytec3d::fem::FECatalog &catalog() const { return *catalog_; }
    /** @brief The quadrature of the local elements (batched), and the reference shape tables of their elements. */
    const pytec3d::fem::QuadratureData &quadrature() const { return *quadrature_; }
    const pytec3d::fem::ReferenceTables &reference_tables() const { return reference_tables_; }

    /** @brief Reads the model data (model_data.h) and builds the element types it selects. */
    void load_model_data(const nlohmann::json &model_data);

    /**
     * @brief Dry run of the coefficients on one element per (TYPE, MATERIAL, CONSTANT) combination of the mesh.
     *
     * Aborts, naming the element type, domain and material or constants, when a term active in the domain
     * needs a material property or real constant the model data leaves undefined.
     */
    void validate_model_data() const;

    const pytec3d::ModelData &get_model_data() const { return model_data_; }
    const SolverSettings &get_solver_settings() const { return solver_settings_; }

    /** @brief The variables of the three systems, the DOFs and the DOF system. */
    void initialize_systems();

    /** @brief Updates time-step state and decides whether this solve advances time. */
    void prepare_time_step_solve(const TimeStepInput &time_step);

    /** @brief Solves one time step with @p settings: constraints, assembly, solve, nodal and element fields. */
    void solve_equation_systems(const TimeStepInput &time_step, const SolverSettings &settings);

    /** @brief Writes the fields of @p system into the caller's arrays (collective). */
    void write_output(SystemType system, const pytec3d::ModelArrays &arrays) const;
    /** @brief The arrays @p system writes that @p arrays lack ("name (point)" / "name (cell)"). */
    std::vector<std::string> missing_outputs(SystemType system, const pytec3d::ModelArrays &arrays) const;

    bool uses_nonlinear_system() const { return use_nonlinear_system_; }
    DofSystem &dof_system();
    const DofSystem &dof_system() const;
    const pytec3d::fem::DofSpace &dof_space() const { return *dof_space_; }
    const pytec3d::fem::SystemLayout &nodal_layout() const { return *nodal_layout_; }
    const pytec3d::fem::SystemLayout &element_layout() const { return *element_layout_; }
    const pytec3d::fem::SystemLayout &layout(SystemType type) const;

    StepReport &step_report() { return step_report_; }
    const StepReport &step_report() const { return step_report_; }
    const TimeStep *current_time_step() const { return time_step_.get(); }
    const TimeStep *previous_time_step() const { return old_time_step_.get(); }

    /** @brief True when a law is selected on an element of the grid: the DOF system is solved with Newton. */
    bool has_nonlinear_elements() const;
    /** @brief Highest time order of the element types of the model data: 0 static, 1, or 2. */
    unsigned int time_order() const;
    /** @brief Scheme of the current time step for the time order of the model (time_step.h). */
    TimeScheme time_scheme() const;
    /** @brief The element types of the model, by TYPE id; aborts before load_model_data(). */
    const ModelElementTypes &element_types() const;

    /**
     * @brief Takes the local-to-global rotations of the elements and of the points from the arrays of the caller
     * (model_input.h: (n_cells, 9) and (n_points, 9); empty for identity).
     */
    void build_rotation_caches(const pytec3d::ConstRealArray &element_rotations, const pytec3d::ConstRealArray &node_rotations);
    /** @brief Element axes to global axes of local element @p e. */
    const pytec3d::la::RealTensorValue &element_rotation(int e) const { return element_rotations_[static_cast<std::size_t>(e)]; }
    /** @brief Node frame to global axes of each point; empty: identity. */
    const std::vector<pytec3d::la::RealTensorValue> &node_rotations() const { return node_rotations_; }

private:
    std::vector<pytec3d::fem::Variable> system_variables(SystemType type) const;
    /** @brief The largest input of each nonlinear law per id, into the step report; warns on extrapolation. */
    void compute_law_ranges();

    pytec3d::fem::FEMesh &mesh_;
    pyTEC3DModel &model_;
    std::unique_ptr<pytec3d::fem::FECatalog> catalog_;
    std::unique_ptr<pytec3d::fem::QuadratureData> quadrature_;
    pytec3d::fem::ReferenceTables reference_tables_;
    bool use_nonlinear_system_ = false;
    std::unique_ptr<TimeStep> time_step_;
    std::unique_ptr<TimeStep> old_time_step_;
    StepReport step_report_;
    std::unique_ptr<ModelElementTypes> element_types_;
    pytec3d::ModelData model_data_;
    SolverSettings solver_settings_;
    std::vector<pytec3d::la::RealTensorValue> element_rotations_;
    std::vector<pytec3d::la::RealTensorValue> node_rotations_;

    std::unique_ptr<pytec3d::fem::DofSpace> dof_space_;
    std::unique_ptr<pytec3d::fem::SystemLayout> nodal_layout_;
    std::unique_ptr<pytec3d::fem::SystemLayout> element_layout_;
    std::unique_ptr<DofSystem> dof_system_;

};
