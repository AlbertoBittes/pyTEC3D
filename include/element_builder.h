/**
 * @file element_builder.h
 * @brief Per-element assembly: reinits the FEs, resizes the matrices and calls the kernels of the element type.
 *
 * ElementBuilder collects what an element needs (material and real constants, local-to-global rotation,
 * curves, VTK loads, local solution), finds where each field sits in the element vectors of the systems
 * (ElementType::element_indices) and has an ElementEvaluator run the evaluations of the element type:
 *   - Stiffness (Ke) and time-derivative (Ce) matrices of the DOF system.
 *   - Source (RHS) vector, and the tangent matrix / residual of the laws selected on the element.
 *   - Nodal and element post-processing fields, written to the solution systems.
 *
 * The formulations are those of element_types.json (formulations.h, generated).
 */
#pragma once

#include <string>
#include <map>
#include <unordered_map>
#include <vector>

#include <libmesh/dense_matrix.h>
#include <libmesh/dense_vector.h>
#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/equation_systems.h>
#include <libmesh/numeric_vector.h>

#include "element_evaluator.h"
#include "element_parameters.h"
#include "element_type.h"
#include "model_element_types.h"
#include "model_data.h"
#include "time_step.h"
#include "variables_map.h"

class pyTEC3DEquationSystems;
class vtkUnstructuredGrid;

namespace libMesh
{
    class ImplicitSystem;
}
class vtkDataArray;

/** @brief Per-element context bundle populated by set_element_data(). */
struct ElemData
{
    const libMesh::Elem *elem = nullptr;               ///< Active element.
    const ModelElementTypes::Entry *entry = nullptr; ///< Element type and domain of the element's subdomain.
    unsigned int material_num = 0;                     ///< Material ID (from VTK MATERIAL cell array).
    unsigned int real_constant_num = 0;                ///< Real-constant ID (from VTK CONSTANT cell array).
    unsigned int coordinate_system_num = 0;            ///< Coordinate-system ID (from VTK COORDINATE_SYSTEM cell array).
    pytec3d::registry::LawMask laws = 0;               ///< Laws selected by its TYPE, MATERIAL and CONSTANT rows.
};

/**
 * @brief Builds per-element matrices, vectors and post-processing fields from the element kernels.
 *
 * Constructed once per assembly pass. Call set_element_data() before any build_*() method.
 */
class ElementBuilder
{
public:
    enum MatricesSet
    {
        DOF_SOLUTION_MATRICES,      ///< {Ke, Ce} and Me for second-order element types
        TANGENT_STIFFNESS_MATRICES,
        NODAL_SOLUTION_MATRICES,    ///< Static nodal fields: {F u, F du/dt, F d2u/dt2} (terms without trial: at output)
        ELEMENT_SOLUTION_MATRICES,  ///< Static element fields: same layout
        NUM_MATRICES_SETS
    };

    enum VectorsSet
    {
        RHS_VECTORS,
        RESIDUAL_VECTORS,
        NUM_VECTORS_SETS
    };

    enum DOFMatrices
    {
        Ke,
        Ce,
        NUM_DOF_MATRICES,
    };

    /** @brief Constructs the builder and caches the model data of the equation systems. */
    ElementBuilder(pyTEC3DEquationSystems *equation_systems);

    /** @brief Populates the ElemData context for element @p elem. */
    void set_element_data(const libMesh::Elem *elem);

    /** @brief Returns the active element context. */
    const ElemData &get_elem_data() const;

    /** @brief True when laws are selected on the active element: it has a Jacobian and residual. */
    bool nonlinear_element() const;

    /** @brief Element matrix sets (indexed by MatricesSet). */
    std::vector<std::vector<libMesh::DenseMatrix<libMesh::Number>>> &get_matrices() { return matrices; }

    /** @brief Element vector sets (indexed by VectorsSet). */
    std::vector<std::vector<libMesh::DenseVector<libMesh::Number>>> &get_vectors() { return vectors; }

    /** @brief Builds Ke, Ce, (second-order element types) Me and the operators of the static nodal and element fields. */
    void build_element();

    /** @brief Builds the source vector (element loads, nodal loads, magnetization). */
    void build_element_rhs();

    /** @brief Builds the tangent matrix and residual of the law terms at the current solution; @p scheme gives the
     * velocity and acceleration of the law terms with dt() and the scaling of their Jacobian. */
    void build_element_tangent_stiffness_residual(const TimeScheme &scheme);

    /** @brief Nodal fields of the active element from the current solution (stored operators, read by the nodal system,
     * plus the dynamic fields), written to the nodal system. */
    void build_nodal_solution();

    /** @brief Element fields of the active element, as build_nodal_solution(), written to the element system. */
    void build_element_solution();

    /** @brief Number of DOFs of the active element. */
    unsigned int get_n_dofs() const { return n_dofs; }

    /** @brief Shares canonical nodal-load injector ownership data with this builder instance. */
    void set_nodal_load_map(const NodalLoadInjectorElemByVarMap *map)
    {
        nodal_load_map_ = map;
    }

private:
    /** @brief Material, constants, rotation and grid data of the active element. */
    pytec3d::ElementParameterValues::Inputs parameter_inputs() const;

    /** @brief Values of @p solution on the DOFs of the active element, in element DOF order (zeros when null). */
    libMesh::DenseVector<libMesh::Number> local_dof_vector(const libMesh::NumericVector<libMesh::Number> *solution) const;

    bool should_inject_nodal_load_dof(const std::string &dof_variable_name, libMesh::dof_id_type dof_id) const;

    /** @brief Evaluates the nodal or element fields on the DOF solution and writes the result into @p system. */
    void build_fields(pytec3d::FieldKind kind, libMesh::System &system);


    /** @brief Inputs of an ElementEvaluator on the active element. */
    pytec3d::ElementEvaluator::Inputs evaluator_inputs(pytec3d::Operation operation,
                                                       pytec3d::ElementParameterValues &values,
                                                       const pytec3d::ElementType::FieldIndices &dofs) const;

    std::vector<std::vector<libMesh::DenseMatrix<libMesh::Number>>> matrices; ///< By MatricesSet.
    std::vector<std::vector<libMesh::DenseVector<libMesh::Number>>> vectors;  ///< By VectorsSet.

    ElemData elem_data;
    unsigned int n_dofs = 0; ///< Number of DOFs on the active element.
    std::vector<libMesh::dof_id_type> dof_indices_; ///< DOF system indices of the active element.
    const NodalLoadInjectorElemByVarMap *nodal_load_map_ = nullptr;

    pyTEC3DEquationSystems *equation_systems_ = nullptr;
    libMesh::ImplicitSystem *dof_solution_system_ = nullptr;
    libMesh::System *nodal_solution_system_ = nullptr;
    libMesh::System *element_solution_system_ = nullptr;
    libMesh::DofMap *dof_solution_dof_map_ = nullptr;

    const ModelElementTypes *element_types_ = nullptr;
    const pytec3d::ModelData *model_data_ = nullptr;
    vtkUnstructuredGrid *vtk_grid_ = nullptr;
    vtkDataArray *material_cell_array_ = nullptr;
    vtkDataArray *real_constant_cell_array_ = nullptr;
    vtkDataArray *coordinate_system_cell_array_ = nullptr;
    const std::unordered_map<libMesh::dof_id_type, libMesh::TensorValue<libMesh::Real>> *element_local_to_global_rotations_ = nullptr;
};
