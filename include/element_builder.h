/**
 * @file element_builder.h
 * @brief Core element-level FE matrix and vector builder.
 *
 * ElementBuilder is the heart of the APEC-M discretisation layer.  Given
 * material and constant data from the Python model, it evaluates shape
 * functions at Gauss points and assembles per-element:
 *   - Stiffness (Ke) and capacitance (Ce) matrices for the DOF system.
 *   - Base and load-modified RHS vectors.
 *   - Post-processing matrices for nodal (B, H) and element (Be, He, Jt) fields.
 *
 * Results are serialised to rank-local @c .emat files by BinaryReader.
 */
#ifndef ELEMENT_BUILDER_H
#define ELEMENT_BUILDER_H

#include <filesystem>
#include <memory>
#include <array>
#include <unordered_map>
#include <set>
#include <vector>

#include <libmesh/elem.h>
#include <libmesh/dense_matrix.h>
#include <libmesh/dense_vector.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/equation_systems.h>
#include <libmesh/dof_map.h>
#include <libmesh/tensor_value.h>
#include "fe_parent_object.h"
#include "apec_m_visibility.h"
#include "element_types_factory.h"
#include "variables_map.h"

class ApecMLinearSystem;
class ApecMEquationSystems;
class vtkUnstructuredGrid;
class vtkDataArray;

/** @brief Per-element context bundle populated by set_element_data(). */
struct ElemData
{
    libMesh::Elem *elem = nullptr;       ///< Pointer to the active element.
    ElementType *element_type = nullptr; ///< Pointer to the element-type descriptor.
    unsigned int material_num = 0;       ///< Material ID (from VTK MATERIAL cell array).
    unsigned int real_constant_num = 0;  ///< Real-constant ID (from VTK CONSTANT cell array).
};

/**
 * @brief Column indices into the cached material property row.
 * @details MUR_{X,Y,Z} = relative permeability components; RHO_{X,Y,Z} = resistivity components.
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
 * RAD = coil radius; T_{X,Y,Z} = current-direction unit vector; R = resistance; S = number of strands.
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
 * @brief Builds per-element FE matrices, RHS vectors, and post-processing matrices.
 *
 * Constructed once per assembly pass.  Call set_element_data() before any
 * build_*() method.  Owns the DenseMatrix / DenseVector storage organised
 * by MatricesSet and VectorsSet.  BinaryReader serialises the populated
 * containers to @c .emat files.
 */
class BinaryReader;

class APEC_M_HIDDEN ElementBuilder
    : public FEParentObjectBase<ApecMEquationSystems>
{
public:
    enum MatricesSet
    {
        DOF_SOLUTION_MATRICES,
        TANGENT_STIFFNESS_MATRICES,
        NODAL_SOLUTION_MATRICES,
        ELEMENT_SOLUTION_MATRICES,
        NUM_MATRICES_SETS
    };

    enum VectorsSet
    {
        RHS_VECTORS,
        NODAL_SOLUTION_VECTORS,
        ELEMENT_SOLUTION_VECTORS,
        NUM_VECTORS_SETS
    };

    enum DOFMatrices
    {
        Ke,
        Ce,
        NUM_DOF_MATRICES,
    };
    static constexpr unsigned int kNumNodalSolutionMatrices = 2;
    static constexpr unsigned int kNumElementSolutionMatrices = 3;
    static constexpr unsigned int kNodalMe = 0;
    static constexpr unsigned int kNodalAAe = 1;
    static constexpr unsigned int kElementJAe = 0;
    static constexpr unsigned int kElementJVe = 1;
    static constexpr unsigned int kElementBAe = 2;
    inline static const libMesh::Real kMu0 = 4.0 * libMesh::pi * 1e-7;

    /** @brief Resizes the vector storage for @p vectors_set according to @p sizes_map. */
    void build_vectors(std::vector<std::vector<unsigned int>> *, unsigned int vectors_set);

    /** @brief Resizes the matrix storage for @p matrices_set according to @p sizes_map. */
    void build_matrices(std::vector<std::vector<unsigned int>> *, unsigned int matrices_set);

    /** @brief Allocates all element matrix and vector containers for the active element. */
    void prepare_element_matrix_building();

    /** @brief Populates the ElemData context for element @p elem. */
    void set_element_data(libMesh::Elem *elem);

    /** @brief Returns a mutable reference to the active element context. */
    ElemData &get_elem_data();

    /** @brief Returns a const reference to the active element context. */
    const ElemData &get_elem_data() const;

    /** @brief Returns all element matrix sets (indexed by MatricesSet). */
    std::vector<std::vector<libMesh::DenseMatrix<libMesh::Number>>> &get_matrices() { return matrices; }

    /** @brief Returns all element vector sets (indexed by VectorsSet). */
    std::vector<std::vector<libMesh::DenseVector<libMesh::Number>>> &get_vectors() { return vectors; }

    /** @brief Dispatches element-matrix build to the correct formulation (AV3DNODAL / AV3DEDGE / AV2DNODAL). */
    void build_element();
    void build_av3dnodal(); ///< Builds element matrices for the 3D nodal A-V formulation.
    void build_av3dedge();  ///< Builds element matrices for the 3D Nedelec edge A-V formulation.
    void build_av2dnodal(); ///< Builds element matrices for the 2D nodal Az-V formulation.

    /** @brief Dispatches RHS-only build to the correct formulation (loads are applied at this stage). */
    void build_element_rhs();
    void build_av3dnodal_rhs(); ///< Applies source currents to the 3D nodal RHS.
    void build_av3dedge_rhs();  ///< Applies source currents to the 3D edge RHS.
    void build_av2dnodal_rhs(); ///< Applies source currents to the 2D nodal RHS.

    /** @brief Dispatches nodal-solution post-processing to the correct formulation. */
    void build_nodal_solution();

    /** @brief Dispatches element-solution post-processing to the correct formulation. */
    void build_element_solution();

    void build_av3dnodal_nodal_solution(); ///< Computes B and H at nodes for AV3DNODAL.
    void build_av3dedge_nodal_solution();  ///< Computes B and H at nodes for AV3DEDGE.
    void build_av2dnodal_nodal_solution(); ///< Computes B and H at nodes for AV2DNODAL.

    void build_av3dnodal_element_solution(); ///< Computes Be, He, Jt at element centroid for AV3DNODAL.
    void build_av3dedge_element_solution();  ///< Computes Be, He, Jt at element centroid for AV3DEDGE.
    void build_av2dnodal_element_solution(); ///< Computes Be, He, Jt at element centroid for AV2DNODAL.

    /** @brief Dispatches element-matrix build to the correct formulation (AV3DNODAL / AV3DEDGE / AV2DNODAL). */
    void build_element_tangent_stiffness();
    void build_av3dnodal_tangent_stiffness(); ///< Builds tangent stiffness matrix for the 3D nodal A-V formulation.
    void build_av3dedge_tangent_stiffness();  ///< Builds tangent stiffness matrix for the 3D Nedelec edge A-V formulation.
    void build_av2dnodal_tangent_stiffness(); ///< Builds tangent stiffness matrix for the 2D nodal Az-V formulation.
    /**
     * @brief Returns the per-variable local solution vectors for @p dof_map_index system.
     * @param solution       Global solution vector to extract local values from.
     * @param dof_map_index  ApecMSystemBase::SystemType index selecting which DOF map to use.
     */
    std::vector<libMesh::DenseVector<libMesh::Number>> get_variables_solution(
        const libMesh::NumericVector<libMesh::Number> *solution, int dof_map_index);

    /** @brief Sets the number of DOFs for the active element. */
    void set_n_dofs(unsigned int n) { n_dofs = n; }

    /** @brief Returns the number of DOFs for the active element. */
    unsigned int get_n_dofs() const { return n_dofs; }

    /** @brief Clears the cached DOF index set used to avoid duplicate nodal load applications. */
    void clear_nodal_loads_dof_cache() { nodal_loads_dof_cache.clear(); }

    /** @brief Constructs the builder and caches material/constant data from the Python model. */
    ElementBuilder(ApecMEquationSystems *equation_systems);

private:
    friend class BinaryReader;
    friend class AssemblerLinear;

    // Slot 0: stiffness/capacitance matrices (Ke, Ce).
    // Slot 1: nodal-solution matrices (Me, AAe used to interpolate B from A).
    // Slot 2: element-solution matrices (JAe, JVe, BAe for element-averaged fields).
    std::vector<std::vector<libMesh::DenseMatrix<libMesh::Number>>> matrices; ///< Per-element matrix containers (by MatricesSet).
    // Slot 0: element RHS vectors (Fe).
    std::vector<std::vector<libMesh::DenseVector<libMesh::Number>>> vectors; ///< Per-element vector containers (by VectorsSet).

    ElemData elem_data;                                   ///< Context for the element currently being processed.
    unsigned int n_dofs = 0;                              ///< Number of DOFs on the active element.
    std::set<libMesh::dof_id_type> nodal_loads_dof_cache; ///< Guards against double-application of nodal loads.

    ApecMEquationSystems *equation_systems_ = nullptr;   ///< Owning equation-systems container.
    ApecMLinearSystem *dof_solution_system_ = nullptr;   ///< DOF solution system reference.
    libMesh::System *nodal_solution_system_ = nullptr;   ///< Nodal post-processing system reference.
    libMesh::System *element_solution_system_ = nullptr; ///< Element post-processing system reference.
    libMesh::System *support_solution_system_ = nullptr; ///< Reserved for future auxiliary system.

    libMesh::DofMap *dof_solution_dof_map_ = nullptr;                                         ///< DOF map for the primary unknown system.
    std::unordered_map<libMesh::dof_id_type, std::vector<libMesh::Real>> materials_rows_;     ///< Cached material property rows keyed by material ID.
    std::unordered_map<libMesh::dof_id_type, std::vector<libMesh::Real>> constants_rows_;     ///< Cached real-constant rows keyed by constant ID.
    vtkUnstructuredGrid *vtk_grid_ = nullptr;                                                 ///< Non-owning VTK grid pointer.
    vtkDataArray *material_cell_array_ = nullptr;                                             ///< VTK MATERIAL cell-data array.
    vtkDataArray *real_constant_cell_array_ = nullptr;                                        ///< VTK CONSTANT cell-data array.
    vtkDataArray *system_cell_array_ = nullptr;                                               ///< VTK SYSTEM cell-data array (0=global).
    vtkDataArray *js_cell_array_ = nullptr;                                                   ///< VTK Js cell-data array (stored in element-local frame).
    std::unordered_map<libMesh::dof_id_type, std::array<libMesh::Real, 9>> system_rotations_; ///< Local-to-global rotation matrices keyed by SYSTEM id.
    std::vector<std::vector<libMesh::Real>> element_material_rows_;                           ///< Per-element material rows rotated to global Cartesian.
    std::vector<std::vector<libMesh::Real>> element_constant_rows_;                           ///< Per-element constants rows (T vector rotated to global Cartesian).
    std::vector<std::array<libMesh::Real, 3>> element_js_global_;                             ///< Per-element Js vectors rotated to global Cartesian.

    template <typename Resolver>
    void fill_vectors_by_api(const std::vector<APIVariableName> &api_vars,
                             std::vector<libMesh::DenseVector<libMesh::Number>> &output_vectors,
                             Resolver &&resolve_vector)
    {
        for (unsigned int k = 0; k < api_vars.size(); ++k)
        {
            const libMesh::DenseVector<libMesh::Number> *vec = resolve_vector(api_vars[k]);
            if (vec)
                output_vectors[k] = *vec;
        }
    }

    void append_and_write_solution(libMesh::System &system,
                                   const std::vector<libMesh::DenseVector<libMesh::Number>> &ordered_vectors);
    void build_coordinate_system_rotation_cache(pybind11::object apec_object);
    void build_per_element_property_cache();
    const std::array<libMesh::Real, 9> &element_rotation_matrix(vtkIdType elem_id) const;
    std::array<libMesh::Real, 3> current_js_global() const;
    static std::array<libMesh::Real, 3> rotate_vector_to_global(
        const std::array<libMesh::Real, 3> &local,
        const std::array<libMesh::Real, 9> &rotation);
    static std::array<libMesh::Real, 3> rotate_diagonal_tensor_to_global_diag(
        const std::array<libMesh::Real, 3> &local_diag,
        const std::array<libMesh::Real, 9> &rotation);
    const std::vector<libMesh::Real> &material_row(libMesh::dof_id_type id) const;
    const std::vector<libMesh::Real> &constant_row(libMesh::dof_id_type id) const;
    libMesh::Real material_value(libMesh::dof_id_type id, unsigned int col_idx) const;
    libMesh::Real constant_value(libMesh::dof_id_type id, unsigned int col_idx) const;
};

#endif // ELEMENT_BUILDER_H
