#ifndef ELEMENT_BUILDER_H
#define ELEMENT_BUILDER_H

#include <filesystem>
#include <memory>
#include <set>
#include <vector>

#include <libmesh/elem.h>
#include <libmesh/dense_matrix.h>
#include <libmesh/dense_vector.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/equation_systems.h>
#include <libmesh/dof_map.h>

#include "fe_parent_object.h"
#include "fem_service.h"
#include "fe_data_object.h"
#include "fem_data_repository.h"
// #include "element_types.h"
#include "element_types_factory.h"
#include "materials_object.h"
#include "real_constants_object.h"
#include "time_step.h"
#include "session.h"
#include "variables_map.h"

using namespace libMesh;

class ApecMLinearSystem;
class ApecMEquationSystems;

struct ElemData
{
    Elem *elem = nullptr;                // Pointer to the element
    ElementType *element_type = nullptr; // Pointer to the element type
    unsigned int material_num = 0;       // Material number
    unsigned int real_constant_num = 0;  // Real constant number
};

struct SolutionVectorSupport
{
    // std::vector<unsigned int> var_nums;
    // std::vector<std::vector<dof_id_type>> dof_indices_var;
    std::vector<DenseVector<Number>> sol_vectors;
    std::vector<std::vector<dof_id_type>> dof_indices_var;
    void clear(unsigned int idx = 0)
    {
        // var_nums.clear();
        // dof_indices_var.clear();
        sol_vectors.clear();
        dof_indices_var.clear();
    }
};

// This is probably the most important class of the program. Given the MAterials, the Loads, etc. it builds the
// Elements' matrices and RHS, or solution vectors (B, J, etc.), and writes the m down to files, which will then be
// Assembled by the assembler.

class BinaryReader;

class ElementBuilder
    : public FEMService<ElementTypesFactory, MaterialsObject, RealConstantsObject>,
      public FEParentObjectBase<EquationSystems>
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

    using ContextView = TypedContextView<ElementTypesFactory, MaterialsObject, RealConstantsObject>;

    void build_vectors(std::vector<std::vector<unsigned int>> *, unsigned int vectors_set);
    void build_matrices(std::vector<std::vector<unsigned int>> *, unsigned int matrices_set);
    void prepare_element_matrix_building();

    void set_element_data(Elem *elem);

    std::vector<std::vector<DenseMatrix<Number>>> &get_matrices() { return matrices; }
    std::vector<std::vector<DenseVector<Number>>> &get_vectors() { return vectors; }

    // Methods to build elements
    void build_element();
    void build_solid97();
    void build_solid236();
    void build_plane233();
    void build_infin110();
    void build_infin111();
    void build_infin112();

    // Methods to build elements
    void build_element_rhs();
    void build_solid97_rhs();
    void build_solid236_rhs();
    void build_plane233_rhs();
    void build_infin110_rhs();
    void build_infin111_rhs();
    void build_infin112_rhs();

    void build_solution_vectors();
    void build_nodal_solution();
    void build_element_solution();
    void build_support_solution();
    bool build_solution_vectors_from_emat_files(dof_id_type element_id,
                                                const std::filesystem::path &emat_folderpath);

    void build_solid97_solution_vectors();
    void build_solid236_solution_vectors();
    void build_plane233_solution_vectors();
    void build_infin110_solution_vectors();
    void build_infin111_solution_vectors();
    void build_infin112_solution_vectors();

    void build_solid97_nodal_solution();
    void build_solid236_nodal_solution();
    void build_plane233_nodal_solution();
    void build_infin110_nodal_solution();
    void build_infin111_nodal_solution();
    void build_infin112_nodal_solution();

    void build_solid97_element_solution();
    void build_solid236_element_solution();
    void build_plane233_element_solution();
    void build_infin110_element_solution();
    void build_infin111_element_solution();
    void build_infin112_element_solution();

    void build_solid97_support_solution();
    void build_solid236_support_solution();
    void build_plane233_support_solution();
    void build_infin110_support_solution();
    void build_infin111_support_solution();
    void build_infin112_support_solution();

    // WARNING: This fill solution vector support I can do it once, not twice for nodeal and element solution vectors
    void fill_solution_vector_support(const NumericVector<Number> *solution, int dof_map_index);

    bool get_compute_nodal_solution_vectors() const { return compute_nodal_solution_vectors; }
    void set_compute_nodal_solution_vectors(bool compute) { compute_nodal_solution_vectors = compute; }

    void set_n_dofs(unsigned int n) { n_dofs = n; }
    unsigned int get_n_dofs() const { return n_dofs; }

    void clear_nodal_loads_dof_cache() { nodal_loads_dof_cache.clear(); }

    static std::unique_ptr<ElementBuilder> make_element_builder(ApecMEquationSystems *eq_sys);

    ElementBuilder(EquationSystems *equation_systems, std::shared_ptr<ContextView> ctx_view);

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
    SolutionVectorSupport solution_vector_support;
    bool compute_nodal_solution_vectors = true;
    unsigned int n_dofs = 0;
    std::set<dof_id_type> nodal_loads_dof_cache;

    ApecMEquationSystems *equation_systems_ = nullptr;
    ApecMLinearSystem *dof_solution_system_ = nullptr;
    System *nodal_solution_system_ = nullptr;
    System *element_solution_system_ = nullptr;
    System *support_solution_system_ = nullptr;
    System *nodal_loads_system_ = nullptr;
    System *element_loads_system_ = nullptr;
    DofMap *dof_solution_dof_map_ = nullptr;
    MaterialsObject *materials_object_ = nullptr;
    RealConstantsObject *real_constants_object_ = nullptr;

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
};

#endif // ELEMENT_BUILDER_H
