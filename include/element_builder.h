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

class AssemblerLinear;

class ElementBuilder : public FEMServiceAbstract
{
public:
    // Default constructor
    using FEMServiceAbstract::FEMServiceAbstract;

    // Getters
    std::vector<std::vector<DenseMatrix<Number>>> &get_matrices() { return matrices; }
    std::vector<std::vector<DenseVector<Number>>> &get_vectors() { return vectors; }

    inline void set_element_data(Elem *elem, std::shared_ptr<ElementType> element_type, unsigned int material_num, unsigned int real_constant_num)
    {
        this->elem = elem;
        this->element_type = element_type;
        this->material_num = material_num;
        this->real_constant_num = real_constant_num;
        this->element_type->reinit_fe(this->elem);
        this->element_type->reinit_fe_vec(this->elem);
    }
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
        this->buildMatrices(element_type->get_matrix_idx_size_map(), matrices);
        this->buildVectors(element_type->get_vector_idx_size_map(), vectors);
        this->buildVectors(element_type->get_vector_idx_size_map(), prev_vectors);
        this->buildMatrices(element_type->get_derived_matrix_idx_size_map(), matrices_solution);
    }
    void prepare_element_rhs_building()
    {
        this->scalar_shape_functions = element_type->get_scalar_shape_functions();
        this->element_type->reinit_fe(this->elem);
        this->element_type->populate_vector_idx_size_map();
        this->buildVectors(element_type->get_vector_idx_size_map(), vectors);
        this->buildVectors(element_type->get_vector_idx_size_map(), prev_vectors);
    }
    // void ramp_vectors(Substep *time_substep)
    // {
    //     Real y1;
    //     Real y2;
    //     for (unsigned int i = 0; i < vectors[0].size(); i++)
    //     {
    //         y2 = vectors[0](i);
    //         y1 = prev_vectors[0](i);

    //         vectors[0](i) = y1 + ((y2 - y1) / (time_substep->get_parent_object()->get_delta_t())) * (time_substep->get_time() - time_substep->get_parent_object()->get_time_beginning());
    //     }
    // }
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

    void buildVectors(std::unique_ptr<std::vector<std::vector<unsigned int>>> sizes_map, std::vector<DenseVector<Number>> &elem_vectors_vector);
    void buildMatrices(std::unique_ptr<std::vector<std::vector<unsigned int>>> sizes_map, std::vector<DenseMatrix<Number>> &elem_matrices_vector);
    void writeElementMatrices(std::ofstream &outFile, unsigned int matrices_set);
    void readElementMatrices(std::ifstream &inFile, unsigned int matrices_set);
    void writeElementVectors(std::ofstream &outFile, unsigned int vectors_set);
    void readElementVectors(std::ifstream &inFile, unsigned int vectors_set);
    // void readElementData(std::ifstream &inFile, std::vector<DenseMatrix<Number>> &elem_matrices_vector, std::vector<DenseVector<Number>> &elem_vectors_vector);
    // void writeElementData(std::ofstream &outFile, const std::vector<DenseMatrix<Number>> &elem_matrices_vector, const std::vector<DenseVector<Number>> &elem_vectors_vector);

private:
    std::vector<std::vector<DenseMatrix<Number>>> matrices; // Element matrices
    std::vector<std::vector<DenseVector<Number>>> vectors;  // Element vectors
};

#endif // ELEMENT_BUILDER_H