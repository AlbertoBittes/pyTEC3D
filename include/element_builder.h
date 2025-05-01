#ifndef ELEMENT_BUILDER_H
#define ELEMENT_BUILDER_H

#include <libmesh/elem.h>
#include <libmesh/dense_matrix.h>
#include <libmesh/dense_vector.h>
#include <vector>
#include "../include/element_types.h" // Include the header file for ElementTypesMap
#include "../include/load_step.h"

using namespace libMesh;

class ElementBuilder
{
public:
    // Default constructor
    ElementBuilder()
        : load_step(nullptr) {}

    // Constructor with LoadStep
    ElementBuilder(std::shared_ptr<LoadStep> load_step)
        : load_step(load_step) {}

    // Constructor with LoadStep and old LoadStep
    ElementBuilder(std::shared_ptr<LoadStep> load_step, std::shared_ptr<LoadStep> old_load_step)
        : load_step(load_step), old_load_step(old_load_step) {}

    // Setters
    void set_load_step(std::shared_ptr<LoadStep> load_step) { this->load_step = load_step; }
    void set_old_load_step(std::shared_ptr<LoadStep> old_load_step) { this->old_load_step = old_load_step; }
    void set_matrices(const std::vector<DenseMatrix<Number>> &matrices) { this->matrices = matrices; }
    void set_matrices_derived(const std::vector<DenseMatrix<Number>> &matrices_derived) { this->matrices_derived = matrices_derived; }
    void set_vectors(const std::vector<DenseVector<Number>> &vectors) { this->vectors = vectors; }
    void set_prev_vectors(const std::vector<DenseVector<Number>> &prev_vectors) { this->prev_vectors = prev_vectors; }
    void set_elem(Elem *elem) { this->elem = elem; }
    void set_element_type(std::shared_ptr<ElementType> element_type) { this->element_type = element_type; }
    void set_material_num(unsigned int material_num) { this->material_num = material_num; }
    void set_real_constant_num(unsigned int real_constant_num) { this->real_constant_num = real_constant_num; }
    void set_scalar_shape_functions(std::shared_ptr<ScalarShapeFunctions> &scalar_shape_functions)
    {
        this->scalar_shape_functions = std::move(scalar_shape_functions);
    }
    void set_vector_shape_functions(std::shared_ptr<VectorShapeFunctions> &vector_shape_functions)
    {
        this->vector_shape_functions = std::move(vector_shape_functions);
    }
    // Getters
    std::shared_ptr<LoadStep> get_load_step() const { return load_step; }
    std::shared_ptr<LoadStep> get_old_load_step() const { return old_load_step; }
    std::vector<DenseMatrix<Number>> &get_matrices() { return matrices; }
    std::vector<DenseMatrix<Number>> &get_matrices_derived() { return matrices_derived; }
    std::vector<DenseVector<Number>> &get_vectors() { return vectors; }
    std::vector<DenseVector<Number>> &get_prev_vectors() { return prev_vectors; }
    Elem *get_elem() const { return elem; }
    std::shared_ptr<ElementType> get_element_type() const { return element_type; }
    unsigned int get_material_num() const { return material_num; }
    unsigned int get_real_constant_num() const { return real_constant_num; }
    std::shared_ptr<ScalarShapeFunctions> get_scalar_shape_functions() const { return scalar_shape_functions; }
    std::shared_ptr<VectorShapeFunctions> get_vector_shape_functions() const { return vector_shape_functions; }

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
        this->buildMatrices(element_type->get_derived_matrix_idx_size_map(), matrices_derived);
    }
    void prepare_element_rhs_building()
    {
        this->scalar_shape_functions = element_type->get_scalar_shape_functions();
        this->element_type->reinit_fe(this->elem);
        this->element_type->populate_vector_idx_size_map();
        this->buildVectors(element_type->get_vector_idx_size_map(), vectors);
    }
    void ramp_vectors(Substep *time_substep)
    {
        Real y1;
        Real y2;
        for (unsigned int i = 0; i < vectors[0].size(); i++)
        {
            y2 = vectors[0](i);
            y1 = prev_vectors[0](i);

            vectors[0](i) = y1 + ((y2 - y1) / (time_substep->get_parent().get_delta_t())) * (time_substep->get_time() - time_substep->get_parent().get_time_beginning());
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

    void buildVectors(const std::map<std::string, std::tuple<unsigned int, unsigned int>> &sizes_map, std::vector<DenseVector<Number>> &elem_vectors_vector);
    void buildMatrices(const std::map<std::string, std::tuple<unsigned int, unsigned int, unsigned int>> &sizes_map, std::vector<DenseMatrix<Number>> &elem_matrices_vector);
    void writeElementMatrices(std::ofstream &outFile, std::vector<DenseMatrix<Number>> &elem_matrices_vector);
    void readElementMatrices(std::ifstream &inFile, std::vector<DenseMatrix<Number>> &elem_matrices_vector);
    void writeElementVectors(std::ofstream &outFile, std::vector<DenseVector<Number>> &elem_matrices_vector);
    void readElementVectors(std::ifstream &inFile, std::vector<DenseVector<Number>> &elem_matrices_vector);
    // void readElementData(std::ifstream &inFile, std::vector<DenseMatrix<Number>> &elem_matrices_vector, std::vector<DenseVector<Number>> &elem_vectors_vector);
    // void writeElementData(std::ofstream &outFile, const std::vector<DenseMatrix<Number>> &elem_matrices_vector, const std::vector<DenseVector<Number>> &elem_vectors_vector);

private:
    std::shared_ptr<LoadStep> load_step;     // Pointer to the associated LoadStep
    std::shared_ptr<LoadStep> old_load_step; // Pointer to the previous LoadStep

    std::vector<DenseMatrix<Number>> matrices;         // Element matrices
    std::vector<DenseMatrix<Number>> matrices_derived; // Derived element matrices
    std::vector<DenseVector<Number>> vectors;          // Element vectors
    std::vector<DenseVector<Number>> prev_vectors;     // Derived element vectors

    Elem *elem;                                // Pointer to the current element being processed
    std::shared_ptr<ElementType> element_type; // Pointer to the current element type
    unsigned int material_num;                 // Current material number
    unsigned int real_constant_num;            // Current real constant number

    std::shared_ptr<ScalarShapeFunctions> scalar_shape_functions; // Pointer to the shape functions
    std::shared_ptr<VectorShapeFunctions> vector_shape_functions; // Pointer to the vector shape functions
};

#endif // ELEMENT_BUILDER_H