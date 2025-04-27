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
    explicit ElementBuilder(LoadStep *load_step)
        : load_step(load_step) {}

    // Setters
    void set_load_step(LoadStep *load_step) { this->load_step = load_step; }
    void set_matrices(const std::vector<DenseMatrix<Number>> &matrices) { this->matrices = matrices; }
    void set_matrices_derived(const std::vector<DenseMatrix<Number>> &matrices_derived) { this->matrices_derived = matrices_derived; }
    void set_vectors(const std::vector<DenseVector<Number>> &vectors) { this->vectors = vectors; }

    // Getters
    LoadStep *get_load_step() const { return load_step; }
    std::vector<DenseMatrix<Number>> &get_matrices() { return matrices; }
    std::vector<DenseMatrix<Number>> &get_matrices_derived() { return matrices_derived; }
    std::vector<DenseVector<Number>> &get_vectors() { return vectors; }

    // Methods to build elements
    void build_element(Elem *elem, std::shared_ptr<ElementType> element_type);
    void build_solid97(Elem *elem, std::shared_ptr<ElementType> element_type);
    void build_solid236(Elem *elem, std::shared_ptr<ElementType> element_type);
    void build_plane13(Elem *elem, std::shared_ptr<ElementType> element_type);

    // Methods to apply loads
    void apply_load(Elem *elem, std::shared_ptr<ElementType> element_type);
    void apply_load_solid97(Elem *elem, std::shared_ptr<ElementType> element_type);
    void apply_load_solid236(Elem *elem, std::shared_ptr<ElementType> element_type);
    void apply_load_plane13(Elem *elem, std::shared_ptr<ElementType> element_type);

private:
    LoadStep *load_step;                               // Pointer to the associated LoadStep
    std::vector<DenseMatrix<Number>> matrices;         // Element matrices
    std::vector<DenseMatrix<Number>> matrices_derived; // Derived element matrices
    std::vector<DenseVector<Number>> vectors;          // Element vectors
};

#endif // ELEMENT_BUILDER_H