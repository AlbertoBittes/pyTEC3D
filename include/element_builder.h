#ifndef ELEMENT_BUILDER_H
#define ELEMENT_BUILDER_H

#include <vector>

#include <libmesh/elem.h>
#include <libmesh/dense_matrix.h>
#include <libmesh/dense_vector.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/equation_systems.h>
#include <libmesh/dof_map.h>

#include "fe_parent_object.h"
#include "fem_service_abstract.h"
#include "fe_data_object.h"
#include "fem_data_repository.h"
// #include "element_types.h"
#include "element_types_factory.h"
#include "materials_object.h"
#include "real_constants_object.h"
#include "time_step.h"
#include "session.h"
#include "globals.h"

using namespace libMesh;

class ApecMLinearSystem;

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
    : public FEMServiceWithOldContext<ElementTypesFactory, MaterialsObject, RealConstantsObject>,
      public FEParentObjectBase<EquationSystems>
{
public:
    enum MatricesSet
    {
        DOF_SOLUTION_MATRICES,
        NODAL_ELEMENT_SOLUTION_MATRICES,
        NUM_MATRICES_SETS
    };

    enum VectorsSet
    {
        RHS_VECTORS,
        NUM_VECTORS_SETS
    };

    enum SOLID97DOFMatrices
    {
        Ke_SOLID97,
        Ce_SOLID97,
        NUM_SOLID97_DOF_MATRICES,
    };

    enum SOLID236DOFMatrices
    {
        Ke_SOLID236,
        Ce_SOLID236,
        NUM_SOLID236_DOF_MATRICES
    };

    enum PLANE13DOFMatrices
    {
        Ke_PLANE13,
        Ce_PLANE13,
        NUM_PLANE13_DOF_MATRICES
    };

    enum SOLID97SolutionMatrices
    {
        Me_SOLID97,
        AAe_SOLID97,
        JAe_SOLID97,
        JVe_SOLID97,
        BAe_SOLID97,
        NUM_SOLID97_SOLUTION_MATRICES
    };

    enum SOLID236SolutionMatrices
    {
        Me_SOLID236,
        AAe_SOLID236,
        JAe_SOLID236,
        JVe_SOLID236,
        BAe_SOLID236,
        NUM_SOLID236_SOLUTION_MATRICES
    };

    enum PLANE13SolutionMatrices
    {
        Me_PLANE13,
        AAe_PLANE13,
        JAe_PLANE13,
        BAe_PLANE13,
        NUM_PLANE13_SOLUTION_MATRICES
    };

    using ContextView = TypedContextView<Session, ElementTypesFactory, TimeStep, MaterialsObject, RealConstantsObject>;
    using FEParentObjectBase<EquationSystems>::FEParentObjectBase;

    ElementBuilder(EquationSystems *equation_systems = nullptr,
                   std::shared_ptr<ContextView> ctx_view = nullptr,
                   std::shared_ptr<ContextView> old_ctx_view = nullptr);

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

    void build_solid97_solution_vectors();
    void build_solid236_solution_vectors();
    void build_plane233_solution_vectors();
    void build_infin110_solution_vectors();
    void build_infin111_solution_vectors();
    void build_infin112_solution_vectors();

    // WARNING: This fill solution vector support I can do it once, not twice for nodeal and element solution vectors
    void fill_solution_vector_support(const NumericVector<Number> *solution, int dof_map_index);

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

    // std::shared_ptr<ElementLoadsObject> element_loads_object = nullptr; // Pointer to ElementLoadsObject
    std::vector<dof_id_type> dof_indices;
};

#endif // ELEMENT_BUILDER_H
