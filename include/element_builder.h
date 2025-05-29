#ifndef ELEMENT_BUILDER_H
#define ELEMENT_BUILDER_H

#include <libmesh/elem.h>
#include <libmesh/dense_matrix.h>
#include <libmesh/dense_vector.h>
#include <vector>
#include "element_types.h" // Include the header file for ElementTypesMap
#include "load_step.h"
#include "globals.h"
#include "fem_service_abstract.h"

using namespace libMesh;

struct ElemData
{
    Elem *elem = nullptr;                // Pointer to the element
    ElementType *element_type = nullptr; // Pointer to the element type
    unsigned int material_num = 0;       // Material number
    unsigned int real_constant_num = 0;  // Real constant number
};

class ElementBuilder : public FEMServiceAbstract // , public FEDataObjectBase<ApecMLinearSystem>
{
public:
    // Default constructor
    using FEMServiceAbstract::FEMServiceAbstract;

    // Getters
    std::vector<std::vector<DenseMatrix<Number>>> &get_matrices() { return matrices; }
    std::vector<std::vector<DenseVector<Number>>> &get_vectors() { return vectors; }

    void prepare_element_matrix_building()
    {
        this->scalar_shape_functions = element_type->get_scalar_shape_functions();
        this->element_type->reinit_fe(this->elem);
        if (element_type->get_fe_vec() != nullptr)
        {
            this->vector_shape_functions = element_type->get_vector_shape_functions();
            this->element_type->reinit_fe_vec(this->elem);
        }
        this->element_type->populate_all();
        this->build_matrices(element_type->get_matrix_idx_size_map(), matrices);
        this->build_vectors(element_type->get_vector_idx_size_map(), vectors);
        this->build_vectors(element_type->get_vector_idx_size_map(), prev_vectors);
        this->build_matrices(element_type->get_derived_matrix_idx_size_map(), matrices_solution);
    }
    void prepare_element_rhs_building()
    {
        this->scalar_shape_functions = element_type->get_scalar_shape_functions();
        this->element_type->reinit_fe(this->elem);
        this->element_type->populate_vector_idx_size_map();
        this->build_vectors(element_type->get_vector_idx_size_map(), 0);
        this->build_vectors(element_type->get_vector_idx_size_map(), 1);
    }
    inline void ramp_vectors(Substep *time_substep)
    {
        Real y1;
        Real y2;
        Real dt = time_substep->get_parent_object()->get_delta_t();
        Real time = time_substep->get_time();
        Real time_beginning = time_substep->get_parent_object()->get_time_beginning();
        Real deltim_subst = time - time_beginning;

        for (unsigned int i = 0; i < vectors[0].size(); i++)
        {
            y2 = vectors[0][0](i);
            y1 = vectors[0][1](i);

            vectors[0][0](i) = y1 + (y2 - y1) / dt * deltim_subst;
        }
    }
    // Methods to build elements
    void build_element_matrix();
    void build_solid97_matrix();
    void build_solid236_matrix();
    void build_plane13_matrix();

    // Methods to build elements
    void build_element_rhs();
    void build_solid97_rhs();
    void build_solid236_rhs();
    void build_plane13_rhs();

    void build_vectors(std::unique_ptr<std::vector<std::vector<unsigned int>>> sizes_map, unsigned int vectors_set);
    void build_matrices(std::unique_ptr<std::vector<std::vector<unsigned int>>> sizes_map, unsigned int matrices_set);
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
    void set_element_data(Elem *elem, ElementType *element_type, unsigned int material_num, unsigned int real_constant_num)
    {
        elem_data.elem = elem;
        elem_data.element_type = element_type;
        elem_data.material_num = material_num;
        elem_data.real_constant_num = real_constant_num;
    }

private:
    std::vector<std::vector<DenseMatrix<Number>>> matrices; // Element matrices
    std::vector<std::vector<DenseVector<Number>>> vectors;  // Element vectors
    ElemData elem_data;                                     // Element data structure
};

#endif // ELEMENT_BUILDER_H