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
#include "materials_object.h"
#include "real_constants_object.h"
#include "loads_object.h"
#include "time_step.h"
#include "session.h"
#include "globals.h"
#include "system_type.h"

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
    bool transient_effect = false;                               // Flag for transient effects
    std::shared_ptr<MaterialsObject> materials_object = nullptr; // Pointer to MaterialsObject
    std::shared_ptr<LoadsObject> loads_object = nullptr;         // Pointer to LoadsObject
    std::shared_ptr<LoadsObject> old_loads_object = nullptr;     // Pointer to LoadsObject
    ExplicitSystem *loads_map = nullptr;                         // Pointer to ExplicitSystem
    ExplicitSystem *old_loads_map = nullptr;                     // Shared pointer to ExplicitSystem
    unsigned int loads_sys_number = 0;                           // System number for loads
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
class ElementBuilder
    : public FEMServiceAbstract<Session, ElementTypeFactory, TimeStep, MaterialsObject, RealConstantsObject, LoadsObject>,
      public FEParentObjectBase<EquationSystems>
{
public:
    using ContextView = TypedContextView<Session, ElementTypeFactory, TimeStep, MaterialsObject, RealConstantsObject, LoadsObject>;
    using FEParentObjectBase<EquationSystems>::FEParentObjectBase;

    ElementBuilder(std::shared_ptr<ContextView> ctx_view = nullptr,
                   EquationSystems *equation_systems = nullptr);

    std::vector<std::vector<DenseMatrix<Number>>> &get_matrices() { return matrices; }
    std::vector<std::vector<DenseVector<Number>>> &get_vectors() { return vectors; }

    void prepare_element_matrix_building()
    {
        // reinit to be done inside each element
        // as well as the reception of shape functions
        // this->elem_data.element_type->reinit_fe(this->elem);

        this->build_matrices(elem_data.element_type->get_dof_matrix_idx_size_map().get(), 0);

        this->build_matrices(elem_data.element_type->get_nodal_solution_matrix_idx_size_map().get(), 1);
        this->build_matrices(elem_data.element_type->get_element_solution_matrix_idx_size_map().get(), 2);

        this->build_vectors(elem_data.element_type->get_vector_idx_size_map().get(), 0);
        // this->build_vectors(vector_idx_size_map, 1);
    }

    void prepare_solution_vectors()
    {
        // reinit to be done inside each element
        // as well as the reception of shape functions
        // this->elem_data.element_type->reinit_fe(this->elem);
        this->build_vectors(elem_data.element_type->get_nodal_solution_vector_idx_size_map().get(), 2);
        this->build_vectors(elem_data.element_type->get_element_solution_vector_idx_size_map().get(), 3);
        this->build_vectors(elem_data.element_type->get_element_solution_vector_idx_size_map().get(), 4);
        // this->build_vectors(elem_data.element_type->get_vector_idx_size_map().get(), 0);
        // this->build_vectors(vector_idx_size_map, 1);
    }
    inline void ramp_vectors(Substep *time_substep, unsigned int set_vector_1, unsigned int set_vector_2, unsigned int idx_vector_1, unsigned int idx_vector_2)
    {
        Real y1;
        Real y2;
        Real dt = time_substep->get_parent_object()->get_delta_t();
        Real time = time_substep->get_time();
        Real time_beginning = time_substep->get_parent_object()->get_time_beginning();
        Real deltim_subst = time - time_beginning;
        unsigned int vec_size = vectors[set_vector_1][idx_vector_1].size();
        for (unsigned int i = 0; i < vec_size; i++)
        {
            y2 = vectors[set_vector_1][idx_vector_1](i);
            y1 = vectors[set_vector_2][idx_vector_2](i);

            vectors[set_vector_1][idx_vector_1](i) = y1 + (y2 - y1) / dt * deltim_subst;
        }
    }

    // Methods to build elements
    void build_element();
    void build_solid97();
    void build_solid236();
    void build_plane13();

    // Methods to build elements
    void build_element_rhs();
    void build_solid97_rhs();
    void build_solid236_rhs();
    void build_plane13_rhs();

    void build_vectors(std::vector<std::vector<unsigned int>> *, unsigned int vectors_set);
    void build_matrices(std::vector<std::vector<unsigned int>> *, unsigned int matrices_set);
    void write_element_matrices(std::ofstream &outFile, unsigned int matrices_set);
    void read_element_matrices(std::ifstream &inFile, unsigned int matrices_set);
    void write_element_vectors(std::ofstream &outFile, unsigned int vectors_set);
    void read_element_vectors(std::ifstream &inFile, unsigned int vectors_set);
    // void readElementData(std::ifstream &inFile, std::vector<DenseMatrix<Number>> &elem_matrices_vector, std::vector<DenseVector<Number>> &elem_vectors_vector);
    // void writeElementData(std::ofstream &outFile, const std::vector<DenseMatrix<Number>> &elem_matrices_vector, const std::vector<DenseVector<Number>> &elem_vectors_vector);
    ElemData &get_element_data()
    {
        return elem_data;
    }
    inline void set_element_data(Elem *elem)
    {
        auto element_types = this->context<ElementTypeFactory>();

        elem_data.element_type = element_types->get_element_type(elem->get_extra_integer(mesh_data_index.etype_num_index));
        elem_data.material_num = elem->get_extra_integer(mesh_data_index.material_num_index);
        elem_data.real_constant_num = elem->get_extra_integer(mesh_data_index.real_constants_num_index);

        elem_data.elem = elem;
    }

    inline void set_elem_builder_misc_data()
    {
        elem_builder_misc_data.time_substep = this->context<TimeStep>()->get_active_substep();
        elem_builder_misc_data.dt = elem_builder_misc_data.time_substep->get_dt();
        elem_builder_misc_data.ramped = elem_builder_misc_data.time_substep->get_ramped_load();
        elem_builder_misc_data.transient_effect = elem_builder_misc_data.time_substep->get_transient_effect();
        if (elem_builder_misc_data.time_substep->get_index() == 0)
        {
            elem_builder_misc_data.old_loads_object = elem_builder_misc_data.loads_object;
            elem_builder_misc_data.old_loads_map = elem_builder_misc_data.loads_map;
            elem_builder_misc_data.loads_object = this->context<LoadsObject>();
            elem_builder_misc_data.loads_map = elem_builder_misc_data.loads_object->get_data_map();

            elem_builder_misc_data.materials_object = this->context<MaterialsObject>();

            elem_builder_misc_data.loads_sys_number = elem_builder_misc_data.loads_map->number();
            if (elem_builder_misc_data.loads_map->has_variable("Js_x"))
            {
                elem_builder_misc_data.jsx_var_num = elem_builder_misc_data.loads_map->variable_number("Js_x");
            }
            if (elem_builder_misc_data.loads_map->has_variable("Js_y"))
            {
                elem_builder_misc_data.jsy_var_num = elem_builder_misc_data.loads_map->variable_number("Js_y");
            }
            if (elem_builder_misc_data.loads_map->has_variable("Js_z"))
            {
                elem_builder_misc_data.jsz_var_num = elem_builder_misc_data.loads_map->variable_number("Js_z");
            }
        }
    };

    void build_write_element_matrices_rhs(bool recompute_matrix = true);

    bool check_compute_matrix(std::shared_ptr<TypedContextView<Session, ElementTypeFactory, MaterialsObject, RealConstantsObject, LoadsObject>> old_view = nullptr);

    void build_write_solution_vectors(const std::string sys_name);

    void build_solution_vector(SystemType sys_type);

    void build_nodal_solution_vector();

    void build_element_solution_vector();

    void build_solid97_nodal_solution_vector();
    void build_solid236_nodal_solution_vector();
    void build_plane13_nodal_solution_vector();
    void build_solid97_element_solution_vector();
    void build_solid236_element_solution_vector();
    void build_plane13_element_solution_vector();

    // WARNING: This fill solution vector support I can do it once, not twice for nodeal and element solution vectors
    void fill_solution_vector_support(const NumericVector<Number> *solution, unsigned int vectors_set);

    void init();

private:
    std::vector<std::vector<DenseMatrix<Number>>> matrices; // Element matrices
    std::vector<std::vector<DenseVector<Number>>> vectors;  // Element vectors
    ElemData elem_data;                                     // Element data structure
    MeshDataIndex mesh_data_index;                          // Mesh data structure
    SolutionVectorSupport solution_vector_support;
    const DofMap *main_dof_map = nullptr;
    ApecMLinearSystem *main_sys = nullptr; // Main system and dof map pointers
    MiscData elem_builder_misc_data;
};

#endif // ELEMENT_BUILDER_H