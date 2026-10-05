/**
 * @file element_builder.h
 * @brief Element-level data gathering on top of the kernel-driven ElementKernels.
 *
 * ElementBuilder collects what the kernels need for one element (material and real-constant rows,
 * local-to-global rotation, BH curve, VTK loads, local solution) and exposes the results as
 * per-element matrices and vectors:
 *   - Stiffness (Ke) and time-derivative (Ce) matrices of the DOF system.
 *   - Source (RHS) vector, and the tangent matrix / residual of nonlinear (BH) elements.
 *   - Nodal and element post-processing fields, written to the solution systems.
 *
 * The kernels themselves, and the per-formulation definitions, live in element_kernels.h and
 * element_types.json.
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
#include "element_kernels.h"
#include "element_types_factory.h"
#include "fe_parent_object.h"
#include "variables_map.h"

class pyTEC3DEquationSystems;
class vtkUnstructuredGrid;
class vtkDataArray;

/** @brief Per-element context bundle populated by set_element_data(). */
struct ElemData
{
    const libMesh::Elem *elem = nullptr;               ///< Active element.
    const ElementTypesFactory::Entry *entry = nullptr; ///< Element type and kernels of the element's subdomain.
    unsigned int material_num = 0;                     ///< Material ID (from VTK MATERIAL cell array).
    unsigned int real_constant_num = 0;                ///< Real-constant ID (from VTK CONSTANT cell array).
    unsigned int coordinate_system_num = 0;            ///< Coordinate-system ID (from VTK COORDINATE_SYSTEM cell array).
};

/**
 * @brief Column indices into the cached material property row.
 * @details MUR_{X,Y,Z} = relative permeability; RHO_{X,Y,Z} = resistivity; M_{X,Y,Z} = magnetization.
 */
enum MaterialColumnIndex : unsigned int
{
    MAT_MUR_X,
    MAT_MUR_Y,
    MAT_MUR_Z,
    MAT_RHO_X,
    MAT_RHO_Y,
    MAT_RHO_Z,
    MAT_M_X,
    MAT_M_Y,
    MAT_M_Z,
    MAT_NUM_COLS
};

/**
 * @brief Column indices into the cached real-constants row.
 * @details THK = thickness; CS = cross section; TURNS = number of turns; VOL = coil volume;
 * RAD = coil radius; T_{X,Y,Z} = current-direction unit vector; R = resistance; S = symmetry factor.
 */
enum ConstantColumnIndex : unsigned int
{
    CONST_THK = 0,
    CONST_CS = 1,
    CONST_TURNS = 2,
    CONST_VOL = 3,
    CONST_RAD = 4,
    CONST_T_X = 5,
    CONST_T_Y = 6,
    CONST_T_Z = 7,
    CONST_R = 8,
    CONST_S = 9,
    CONST_NUM_COLS = 10
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
    pytec3d::ElementParameters element_parameters() const;
    pytec3d::ElementLoads element_loads() const;

    /** @brief Element matrix index of every shape function of each DOF variable, from the DOF system's DofMap. */
    pytec3d::ElementDofs element_dofs() const;

    /** @brief Global DOF of every shape function of @p spec in @p system on the active element, in shape-function order. */
    std::vector<libMesh::dof_id_type> column_dofs(const libMesh::System &system,
                                                  pyTEC3DSystemBase::SystemType system_type,
                                                  const pytec3d::VariableSpec &spec) const;

    /** @brief libMesh variables the active element type has in @p system_type (cached per subdomain). */
    const std::vector<VariablesMap::SystemVariable> &system_variables(pyTEC3DSystemBase::SystemType system_type) const;

    /** @brief Values of @p solution on the DOFs of the active element, in element DOF order (zeros when null). */
    libMesh::DenseVector<libMesh::Number> local_dof_vector(const libMesh::NumericVector<libMesh::Number> *solution) const;

    bool should_inject_nodal_load_dof(const std::string &dof_variable_name, libMesh::dof_id_type dof_id) const;

    /** @brief Writes @p fields (one vector per field of @p kind, in element-type order) into the solution of @p system. */
    void write_fields(libMesh::System &system,
                      pyTEC3DSystemBase::SystemType system_type,
                      pytec3d::FieldKind kind,
                      const std::vector<libMesh::DenseVector<libMesh::Number>> &fields);

    const std::vector<libMesh::Real> &material_row(libMesh::dof_id_type id) const;
    const std::vector<libMesh::Real> &constant_row(libMesh::dof_id_type id) const;

    std::vector<std::vector<libMesh::DenseMatrix<libMesh::Number>>> matrices; ///< By MatricesSet.
    std::vector<std::vector<libMesh::DenseVector<libMesh::Number>>> vectors;  ///< By VectorsSet.

    ElemData elem_data;
    mutable std::map<std::pair<libMesh::subdomain_id_type, int>, std::vector<VariablesMap::SystemVariable>> system_variables_cache_;
    unsigned int n_dofs = 0; ///< Number of DOFs on the active element.
    const NodalLoadInjectorElemByVarMap *nodal_load_map_ = nullptr;

    pyTEC3DEquationSystems *equation_systems_ = nullptr;
    libMesh::ImplicitSystem *dof_solution_system_ = nullptr;
    libMesh::System *nodal_solution_system_ = nullptr;
    libMesh::System *element_solution_system_ = nullptr;
    libMesh::DofMap *dof_solution_dof_map_ = nullptr;

    ElementTypesFactory *element_types_factory_ = nullptr;
    const std::unordered_map<libMesh::dof_id_type, std::vector<libMesh::Real>> *materials_rows_ = nullptr;
    const std::unordered_map<libMesh::dof_id_type, std::vector<libMesh::Real>> *constants_rows_ = nullptr;
    vtkUnstructuredGrid *vtk_grid_ = nullptr;
    vtkDataArray *material_cell_array_ = nullptr;
    vtkDataArray *real_constant_cell_array_ = nullptr;
    vtkDataArray *coordinate_system_cell_array_ = nullptr;
    const std::unordered_map<libMesh::dof_id_type, libMesh::TensorValue<libMesh::Real>> *element_local_to_global_rotations_ = nullptr;
    const BHCurveMaterialsMap *bh_curves_ = nullptr;
};

#endif // ELEMENT_BUILDER_H
