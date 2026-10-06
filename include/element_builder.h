/**
 * @file element_builder.h
 * @brief Per-element assembly: reinits the FEs, resizes the matrices and calls the kernels of the element type.
 *
 * ElementBuilder collects what the kernels need for one element (material and real constants,
 * local-to-global rotation, BH curve, VTK loads, local solution), finds where each field sits in the
 * element matrices (ElementType::element_indices) and calls the kernels of the element type. Results:
 *   - Stiffness (Ke) and time-derivative (Ce) matrices of the DOF system.
 *   - Source (RHS) vector, and the tangent matrix / residual of nonlinear (BH) elements.
 *   - Nodal and element post-processing fields, written to the solution systems.
 *
 * The kernels live in kernels.h, the element types in element_types.json (element_type.h).
 */
#ifndef ELEMENT_BUILDER_H
#define ELEMENT_BUILDER_H

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

#include "bh_curve.h"
#include "element_parameters.h"
#include "element_type.h"
#include "element_types_factory.h"
#include "fe_parent_object.h"
#include "model_data.h"
#include "variables_map.h"

class pyTEC3DEquationSystems;
class vtkUnstructuredGrid;
class vtkDataArray;

/** @brief Per-element context bundle populated by set_element_data(). */
struct ElemData
{
    const libMesh::Elem *elem = nullptr;               ///< Active element.
    const ElementTypesFactory::Entry *entry = nullptr; ///< Element type and domain of the element's subdomain.
    unsigned int material_num = 0;                     ///< Material ID (from VTK MATERIAL cell array).
    unsigned int real_constant_num = 0;                ///< Real-constant ID (from VTK CONSTANT cell array).
    unsigned int coordinate_system_num = 0;            ///< Coordinate-system ID (from VTK COORDINATE_SYSTEM cell array).
};

/**
 * @brief Builds per-element matrices, vectors and post-processing fields from the element kernels.
 *
 * Constructed once per assembly pass. Call set_element_data() before any build_*() method.
 */
class ElementBuilder
    : public FEParentObjectBase<pyTEC3DEquationSystems>
{
public:
    enum MatricesSet
    {
        DOF_SOLUTION_MATRICES, ///< {Ke, Ce}
        TANGENT_STIFFNESS_MATRICES,
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

    /** @brief Returns true when the active element material has an associated BH curve. */
    bool nonlinear_element() const;

    /** @brief Element matrix sets (indexed by MatricesSet). */
    std::vector<std::vector<libMesh::DenseMatrix<libMesh::Number>>> &get_matrices() { return matrices; }

    /** @brief Element vector sets (indexed by VectorsSet). */
    std::vector<std::vector<libMesh::DenseVector<libMesh::Number>>> &get_vectors() { return vectors; }

    /** @brief Builds Ke and Ce of the active element. */
    void build_element();

    /** @brief Builds the source vector (element loads, nodal loads, magnetization). */
    void build_element_rhs();

    /** @brief Builds the tangent matrix and residual of the BH kernels at the current solution. */
    void build_element_tangent_stiffness_residual();

    /** @brief Computes the nodal fields of the active element from the current solution and writes them to the nodal system. */
    void build_nodal_solution();

    /** @brief Computes the element fields of the active element from the current solution and writes them to the element system. */
    void build_element_solution();

    /** @brief Number of DOFs of the active element. */
    unsigned int get_n_dofs() const { return n_dofs; }

    /** @brief Shares canonical nodal-load injector ownership data with this builder instance. */
    void set_nodal_load_map(const NodalLoadInjectorElemByVarMap *map)
    {
        nodal_load_map_ = map;
    }

private:
    class Context;

    /** @brief Material, constants, rotation, grid data and BH curve of the active element. */
    pytec3d::ElementParameterValues::Inputs parameter_inputs() const;

    /** @brief Values of @p solution on the DOFs of the active element, in element DOF order (zeros when null). */
    libMesh::DenseVector<libMesh::Number> local_dof_vector(const libMesh::NumericVector<libMesh::Number> *solution) const;

    bool should_inject_nodal_load_dof(const std::string &dof_variable_name, libMesh::dof_id_type dof_id) const;

    /** @brief Evaluates the nodal or element field kernels on the DOF solution and writes the result into @p system. */
    void build_fields(pytec3d::FieldKind kind, libMesh::System &system);

    /** @brief Out-of-plane thickness of the active 2D element; 1 when its real constants do not define it. */
    libMesh::Real thickness() const;

    std::vector<std::vector<libMesh::DenseMatrix<libMesh::Number>>> matrices; ///< By MatricesSet.
    std::vector<std::vector<libMesh::DenseVector<libMesh::Number>>> vectors;  ///< By VectorsSet.

    ElemData elem_data;
    unsigned int n_dofs = 0; ///< Number of DOFs on the active element.
    const NodalLoadInjectorElemByVarMap *nodal_load_map_ = nullptr;

    pyTEC3DEquationSystems *equation_systems_ = nullptr;
    libMesh::ImplicitSystem *dof_solution_system_ = nullptr;
    libMesh::System *nodal_solution_system_ = nullptr;
    libMesh::System *element_solution_system_ = nullptr;
    libMesh::DofMap *dof_solution_dof_map_ = nullptr;

    ElementTypesFactory *element_types_factory_ = nullptr;
    const pytec3d::ModelData *model_data_ = nullptr;
    vtkUnstructuredGrid *vtk_grid_ = nullptr;
    vtkDataArray *material_cell_array_ = nullptr;
    vtkDataArray *real_constant_cell_array_ = nullptr;
    vtkDataArray *coordinate_system_cell_array_ = nullptr;
    const std::unordered_map<libMesh::dof_id_type, libMesh::TensorValue<libMesh::Real>> *element_local_to_global_rotations_ = nullptr;
    const BHCurveMaterialsMap *bh_curves_ = nullptr;
};

#endif // ELEMENT_BUILDER_H
