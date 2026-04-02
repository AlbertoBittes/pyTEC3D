#ifndef ELEMENT_BUILDER_H
#define ELEMENT_BUILDER_H

#include <filesystem>
#include <memory>
#include <unordered_map>
#include <set>
#include <vector>

#include <libmesh/elem.h>
#include <libmesh/dense_matrix.h>
#include <libmesh/dense_vector.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/equation_systems.h>
#include <libmesh/dof_map.h>
#include <pybind11/pybind11.h>

#include "fe_parent_object.h"
// #include "element_types.h"
#include "element_types_factory.h"
#include "variables_map.h"

using namespace libMesh;

class ApecMLinearSystem;
class ApecMEquationSystems;
class vtkUnstructuredGrid;
class vtkDataArray;

struct ElemData
{
    Elem *elem = nullptr;                // Pointer to the element
    ElementType *element_type = nullptr; // Pointer to the element type
    unsigned int material_num = 0;       // Material number
    unsigned int real_constant_num = 0;  // Real constant number
};

// This is probably the most important class of the program. Given the MAterials, the Loads, etc. it builds the
// Elements' matrices and RHS, or solution vectors (B, J, etc.), and writes the m down to files, which will then be
// Assembled by the assembler.

class BinaryReader;

class ElementBuilder
    : public FEParentObjectBase<EquationSystems>
{
public:
    enum MatricesSet
    {
        DOF_SOLUTION_MATRICES,
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
    inline static const Real kMu0 = 4.0 * libMesh::pi * 1e-7;

    void build_vectors(std::vector<std::vector<unsigned int>> *, unsigned int vectors_set);
    void build_matrices(std::vector<std::vector<unsigned int>> *, unsigned int matrices_set);
    void prepare_element_matrix_building();

    void set_element_data(Elem *elem);
    ElemData &get_elem_data();
    const ElemData &get_elem_data() const;

    std::vector<std::vector<DenseMatrix<Number>>> &get_matrices() { return matrices; }
    std::vector<std::vector<DenseVector<Number>>> &get_vectors() { return vectors; }

    // Methods to build elements
    void build_element();
    void build_av3dnodal();
    void build_av3dedge();
    void build_av2dnodal();

    // Methods to build elements
    void build_element_rhs();
    void build_av3dnodal_rhs();
    void build_av3dedge_rhs();
    void build_av2dnodal_rhs();

    void build_nodal_solution();
    void build_element_solution();

    void build_av3dnodal_nodal_solution();
    void build_av3dedge_nodal_solution();
    void build_av2dnodal_nodal_solution();

    void build_av3dnodal_element_solution();
    void build_av3dedge_element_solution();
    void build_av2dnodal_element_solution();

    std::vector<DenseVector<Number>> get_variables_solution(const NumericVector<Number> *solution, int dof_map_index);

    bool get_compute_nodal_solution_vectors() const { return compute_nodal_solution_vectors; }
    void set_compute_nodal_solution_vectors(bool compute) { compute_nodal_solution_vectors = compute; }

    void set_n_dofs(unsigned int n) { n_dofs = n; }
    unsigned int get_n_dofs() const { return n_dofs; }

    void clear_nodal_loads_dof_cache() { nodal_loads_dof_cache.clear(); }

    ElementBuilder(EquationSystems *equation_systems);

private:
    friend class BinaryReader;
    friend class AssemblerLinear;

    // Matrices' sets
    // Idx 0: element stiffness/capacitance matrices: Ke, Ce, etc.)
    // Idx 1: Nodal solution matrices (e.g. AAe, Me used to compute B from A, JAe, JVe used to compute J from dA/dt and dV/dt)
    std::vector<std::vector<DenseMatrix<Number>>> matrices; // Element matrices
    // Idx 0: elements' RHS (Fe)
    // Idx 1: elements' old RHS at previous LoadStep (used for ramping loads)
    std::vector<std::vector<DenseVector<Number>>> vectors; // Element vectors

    ElemData elem_data; // Element data structure
    bool compute_nodal_solution_vectors = true;
    unsigned int n_dofs = 0;
    std::set<dof_id_type> nodal_loads_dof_cache;

    ApecMEquationSystems *equation_systems_ = nullptr;
    ApecMLinearSystem *dof_solution_system_ = nullptr;
    System *nodal_solution_system_ = nullptr;
    System *element_solution_system_ = nullptr;
    System *support_solution_system_ = nullptr;

    DofMap *dof_solution_dof_map_ = nullptr;
    pybind11::object apec_object_;
    pybind11::object materials_attribute_;
    pybind11::object constants_attribute_;
    vtkUnstructuredGrid *vtk_grid_ = nullptr;
    vtkDataArray *material_cell_array_ = nullptr;
    vtkDataArray *real_constant_cell_array_ = nullptr;

    template <typename Resolver>
    void fill_vectors_by_api(const std::vector<APIVariableName> &api_vars,
                             std::vector<DenseVector<Number>> &output_vectors,
                             Resolver &&resolve_vector)
    {
        for (unsigned int k = 0; k < api_vars.size(); ++k)
        {
            const DenseVector<Number> *vec = resolve_vector(api_vars[k]);
            if (vec)
                output_vectors[k] = *vec;
        }
    }

    void append_and_write_solution(System &system,
                                   const std::vector<DenseVector<Number>> &ordered_vectors);
    pybind11::object get_row_by_id(const pybind11::object &table, const char *id_column, unsigned int id) const;
};

#endif // ELEMENT_BUILDER_H
