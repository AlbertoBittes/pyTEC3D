#ifndef ELEMENT_BUILDER_H
#define ELEMENT_BUILDER_H

#include <vector>

#include <libmesh/elem.h>
#include <libmesh/dense_matrix.h>
#include <libmesh/dense_vector.h>
#include <libmesh/equation_systems.h>

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
#include "apec_m_system_objects.h"

using namespace libMesh;

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

struct SolutionVectorSupport
{
    std::vector<unsigned int> var_nums;
    std::vector<std::vector<dof_id_type>> dof_indices_var;
    std::vector<DenseVector<Number>> var_vectors;
    std::vector<DenseVector<Number>> sol_vectors;
    std::vector<dof_id_type> dof_indices;

    void clear()
    {
        var_nums.clear();
        dof_indices_var.clear();
        var_vectors.clear();
        sol_vectors.clear();
        dof_indices.clear();
    }
};
class ElementBuilder
    : public FEMServiceAbstract<Session, ElementTypeFactory, MaterialsObject, RealConstantsObject, LoadsObject>,
      public FEParentObjectBase<EquationSystems>
{
public:
    using ContextView = TypedContextView<Session, ElementTypeFactory, MaterialsObject, RealConstantsObject, LoadsObject>;
    using FEParentObjectBase<EquationSystems>::FEParentObjectBase;

    ElementBuilder(std::shared_ptr<ContextView> ctx_view = nullptr,
                   EquationSystems *equation_systems = nullptr)
        : FEMServiceAbstract<Session, ElementTypeFactory, MaterialsObject, RealConstantsObject, LoadsObject>(ctx_view),
          FEParentObjectBase<EquationSystems>(equation_systems),
          matrices(3),
          vectors(4)
    {
        MeshBase &mesh = equation_systems->get_mesh();
        mesh_data_index.etype_num_index = mesh.get_elem_integer_index("ansys_etype");
        mesh_data_index.material_num_index = mesh.get_elem_integer_index("ansys_material_type");
        mesh_data_index.real_constants_num_index = mesh.get_elem_integer_index("ansys_real_constant");
        main_sys = &equation_systems->get_system<ApecMLinearSystem>("ApecMLinearSystem");
        main_dof_map = &main_sys->get_dof_map();
    }

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
    void set_element_data(Elem *elem)
    {
        auto element_types = this->context<ElementTypeFactory>();

        elem_data.element_type = element_types->get_element_type(elem->get_extra_integer(mesh_data_index.etype_num_index));
        elem_data.material_num = elem->get_extra_integer(mesh_data_index.material_num_index);
        elem_data.real_constant_num = elem->get_extra_integer(mesh_data_index.real_constants_num_index);

        elem_data.elem = elem;
    }

    void build_write_element_matrices_rhs(bool recompute_matrix = true);

    bool check_compute_matrix(std::shared_ptr<TypedContextView<Session, ElementTypeFactory, MaterialsObject, RealConstantsObject, LoadsObject>> old_view = nullptr);

    void build_write_solution_vectors(SystemType sys_type);

    void build_solution_vector(SystemType sys_type);

    void build_nodal_solution_vector();

    void build_element_solution_vector();

    void build_solid97_nodal_solution_vector();
    void build_solid236_nodal_solution_vector();
    void build_plane13_nodal_solution_vector();
    void build_solid97_element_solution_vector();
    void build_solid236_element_solution_vector();
    void build_plane13_element_solution_vector();

    void fill_solution_vector_support(const NumericVector<Number> *solution)
    {
        solution_vector_support.clear();

        auto active_vars = elem_data.element_type->get_active_variables(Globals::VariableCategory::DOFVariables);
        for (const auto &var_enum : active_vars)
        {
            unsigned int system_var_num = elem_data.element_type->get_variable_information()[var_enum].system_var_num;
            solution_vector_support.var_nums.push_back(system_var_num);

            std::vector<dof_id_type> dof_indices;
            main_dof_map->dof_indices(elem_data.elem, dof_indices, system_var_num);
            solution_vector_support.dof_indices_var.push_back(dof_indices);

            DenseVector<Number> local_sol(dof_indices.size());
            for (unsigned int i = 0; i < dof_indices.size(); ++i)
            {
                local_sol(i) = solution->el(dof_indices[i]);
            }
            solution_vector_support.sol_vectors.push_back(local_sol);
        }
    }

private:
    std::vector<std::vector<DenseMatrix<Number>>> matrices; // Element matrices
    std::vector<std::vector<DenseVector<Number>>> vectors;  // Element vectors
    ElemData elem_data;                                     // Element data structure
    MeshDataIndex mesh_data_index;                          // Mesh data structure
    SolutionVectorSupport solution_vector_support;
    const DofMap *main_dof_map = nullptr;
    ApecMLinearSystem *main_sys = nullptr; // Main system and dof map pointers
};

#endif // ELEMENT_BUILDER_H