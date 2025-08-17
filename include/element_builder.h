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
#include "element_types.h"
#include "element_types_factory.h"
#include "materials_object.h"
#include "real_constants_object.h"
#include "element_loads_object.h"
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

struct MeshDataIndex
{
    unsigned int etype_num_index = 0;          // Element type number index
    unsigned int material_num_index = 0;       // Material number index
    unsigned int real_constants_num_index = 0; // Real constants number index
};

struct MiscData
{
    Real dt = 0.0;
    Substep *time_substep = nullptr; // Pointer to the current time substep
    bool ramped = true;
    bool transient_effect = false;                                          // Flag for transient effects
    std::shared_ptr<MaterialsObject> materials_object = nullptr;            // Pointer to MaterialsObject
    std::shared_ptr<ElementLoadsObject> element_loads_object = nullptr;     // Pointer to ElementLoadsObject
    std::shared_ptr<ElementLoadsObject> old_element_loads_object = nullptr; // Pointer to ElementLoadsObject
    ApecMLoadsSystem *loads_map = nullptr;                                  // Pointer to ApecMLoadsSystem

    ApecMLoadsSystem *nodal_loads_map = nullptr;
    DofMap *dof_map_nodal_loads = nullptr;
    BoundaryInfo *boundary_info = nullptr; // Pointer to BoundaryInfo
    // ExplicitSystem *old_loads_map = nullptr;                     // Shared pointer to ExplicitSystem
    unsigned int loads_sys_number = 0; // System number for loads
    unsigned int jsx_var_num = 0;
    unsigned int jsy_var_num = 0;
    unsigned int jsz_var_num = 0;
};

struct SolutionVectorSupport
{
    // std::vector<unsigned int> var_nums;
    // std::vector<std::vector<dof_id_type>> dof_indices_var;
    std::vector<DenseVector<Number>> sol_vectors[2];

    void clear(unsigned int idx = 0)
    {
        // var_nums.clear();
        // dof_indices_var.clear();
        sol_vectors[idx].clear();
    }
};

// This is probably the most important class of the program. Given the MAterials, the Loads, etc. it builds the
// Elements' matrices and RHS, or solution vectors (B, J, etc.), and writes the m down to files, which will then be
// Assembled by the assembler.

class ElementBuilder
    : public FEMServiceWithOldContext<Session, ElementTypesFactory, TimeStep, MaterialsObject, RealConstantsObject, ElementLoadsObject>,
      public FEParentObjectBase<EquationSystems>
{
public:
    using ContextView = TypedContextView<Session, ElementTypesFactory, TimeStep, MaterialsObject, RealConstantsObject, ElementLoadsObject>;
    using FEParentObjectBase<EquationSystems>::FEParentObjectBase;

    ElementBuilder(EquationSystems *equation_systems = nullptr,
                   std::shared_ptr<ContextView> ctx_view = nullptr,
                   std::shared_ptr<ContextView> old_ctx_view = nullptr);

    void init();

    void write_element_matrices(std::ofstream &outFile, unsigned int matrices_set);
    void read_element_matrices(std::ifstream &inFile, unsigned int matrices_set);
    void write_element_vectors(std::ofstream &outFile, unsigned int vectors_set);
    void read_element_vectors(std::ifstream &inFile, unsigned int vectors_set);

    void build_vectors(std::vector<std::vector<unsigned int>> *, unsigned int vectors_set);
    void build_matrices(std::vector<std::vector<unsigned int>> *, unsigned int matrices_set);
    void prepare_element_matrix_building();
    void prepare_solution_vectors();

    void ramp_vectors(unsigned int set_vector_1, unsigned int set_vector_2, unsigned int idx_vector_1, unsigned int idx_vector_2);

    void set_element_data(Elem *elem);
    void set_elem_builder_misc_data();

    std::vector<std::vector<DenseMatrix<Number>>> &get_matrices() { return matrices; }
    std::vector<std::vector<DenseVector<Number>>> &get_vectors() { return vectors; }

    void build_write_element_matrices_rhs(bool recompute_matrix = true);
    bool check_compute_matrix();

    // Methods to build elements
    void build_element();
    void build_solid97();
    void build_solid236();
    void build_plane13();
    void build_infin110();
    void build_infin111();
    void build_infin112();

    // Methods to build elements
    void build_element_rhs();
    void build_solid97_rhs();
    void build_solid236_rhs();
    void build_plane13_rhs();
    void build_infin110_rhs();
    void build_infin111_rhs();
    void build_infin112_rhs();

    void build_write_solution_vectors(const std::string sys_name);
    void build_solution_vector(Globals::SystemType sys_type);
    void build_nodal_solution_vector();
    void build_element_solution_vector();

    void build_solid97_nodal_solution_vector();
    void build_solid236_nodal_solution_vector();
    void build_plane13_nodal_solution_vector();
    void build_infin110_nodal_solution_vector();
    void build_infin111_nodal_solution_vector();
    void build_infin112_nodal_solution_vector();

    void build_solid97_element_solution_vector();
    void build_solid236_element_solution_vector();
    void build_plane13_element_solution_vector();
    void build_infin110_element_solution_vector();
    void build_infin111_element_solution_vector();
    void build_infin112_element_solution_vector();

    // WARNING: This fill solution vector support I can do it once, not twice for nodeal and element solution vectors
    void fill_solution_vector_support(const NumericVector<Number> *solution, unsigned int vectors_set);

private:
    // Matrices' sets
    // Idx 0: element stiffness/capacitance matrices: Ke, Ce, etc.)
    // Idx 1: Nodal solution matrices (e.g. AAe, Me used to compute B from A)
    // Idx 2: Element solution matrices (e.g. JAe, JVe used to compute J from dA/dt and dV/dt)
    std::vector<std::vector<DenseMatrix<Number>>> matrices; // Element matrices
    // Idx 0: elements' RHS (Fe)
    // Idx 1: elements' old RHS at previous LoadStep (used for ramping loads)
    // Idx 2: Nodal solution vectors (B, ...)
    // Idx 3: Element solution vectors (J, ...)
    // Idx 4: Other support vectors
    std::vector<std::vector<DenseVector<Number>>> vectors; // Element vectors
    std::vector<dof_id_type> dof_indices;
    ElemData elem_data;            // Element data structure
    MeshDataIndex mesh_data_index; // Mesh data structure
    SolutionVectorSupport solution_vector_support;
    const DofMap *main_dof_map = nullptr;
    ApecMLinearSystem *main_sys = nullptr; // Main system and dof map pointers
    MiscData elem_builder_misc_data;
};

#endif // ELEMENT_BUILDER_H