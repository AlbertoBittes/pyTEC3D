#ifndef ELEMENT_BUILDER_H
#define ELEMENT_BUILDER_H

#include <libmesh/elem.h>
#include <libmesh/system.h>
#include <libmesh/dense_matrix.h>
#include <libmesh/dense_vector.h>
#include <libmesh/dense_submatrix.h>
#include <libmesh/dense_subvector.h>
#include <libmesh/dof_map.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/perf_log.h>
#include <libmesh/fe.h>
#include <vector>
#include <map>
#include "../include/time_stepper.h"
#include "../include/element_types.h"        // Include the header file for ElementTypesMap
#include "../include/materials_properties.h" // Add this line to include the definition of MaterialsProperties

using namespace libMesh;

class ElementDataManager;

class ElementBuilder
{
public:
    // Constructor
    ElementBuilder(std::shared_ptr<ElementDataManager> element_data_manager);

    // here for performance, delete all the arguments passing and add the objects as attributes
    void build_element(Elem *elem, std::shared_ptr<ElementType> element_type);

    void build_solid97(Elem *elem, std::shared_ptr<ElementType> element_type);
    void build_solid236(Elem *elem, std::shared_ptr<ElementType> element_type);
    void build_plane13(Elem *elem, std::shared_ptr<ElementType> element_type);

    void apply_load(Elem *elem, std::shared_ptr<ElementType> element_type);
    void apply_load_solid97(Elem *elem, std::shared_ptr<ElementType> element_type);
    void apply_load_solid236(Elem *elem, std::shared_ptr<ElementType> element_type);
    void apply_load_plane13(Elem *elem, std::shared_ptr<ElementType> element_type);

    std::vector<DenseMatrix<Number>> matrices;
    std::vector<DenseMatrix<Number>> matrices_derived;
    std::vector<DenseVector<Number>> vectors;

    std::shared_ptr<ElementDataManager> element_data_manager;
};

#endif // ELEMENT_BUILDER_H