/**
 * @file element_parameters.h
 * @brief Values of the parameters of an element type on one element.
 *
 * The kernels take plain values; this is where they come from. A parameter (element_types.json) is read
 * from the material (mat.), the real constants (rc.) or a cell-data array of the grid (me.), or is an
 * expression of other parameters (compiled C++, registry::evaluate_expression). Vectors and tensors read
 * from the model are given in the element's coordinate system and always put in the global axes with the
 * element rotation, so kernels and expressions only see global quantities. Values are computed when first
 * used and kept for the element, so a property that no kernel of the element's domain needs may stay
 * undefined; one that is needed and undefined aborts with the element type, domain and material or
 * constants it is missing from. A parameter only exists in the domains its element type declares for it.
 *
 * It is also the source of the parameters the nonlinear laws read (formulations.h, evaluate_law()).
 */
#pragma once

#include <vector>

#include <libmesh/libmesh_common.h>
#include <libmesh/tensor_value.h>
#include <libmesh/vector_value.h>
#include <vtkType.h>

#include "element_type.h"
#include "model_data.h"
#include "parameter_value.h"

class vtkCellData;

namespace pytec3d
{
    class ElementParameterValues
    {
    public:
        struct Inputs
        {
            const MaterialRow *material = nullptr;  ///< Null when the model does not define the material.
            const ConstantRow *constants = nullptr; ///< Null when the model does not define the constants.
            unsigned int material_id = 0;
            unsigned int constants_id = 0;
            registry::Domain domain = registry::Domain::COUNT;
            libMesh::TensorValue<libMesh::Real> rotation{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}; ///< Element axes to global axes.
            vtkCellData *cell_data = nullptr;
            vtkIdType element_id = 0;
        };

        ElementParameterValues(const ElementType &type, const Inputs &inputs);

        /** @brief Value of @p id on the element; aborts if the parameter does not exist in its domain or is undefined. */
        const ParameterValue &symbol(registry::ParameterId id);

        // What registry::evaluate_expression() reads.
        const libMesh::Real &get_scalar(registry::ParameterId id) { return symbol(id).scalar(); }
        const libMesh::RealVectorValue &get_vector(registry::ParameterId id) { return symbol(id).vector(); }
        const libMesh::RealTensorValue &get_tensor(registry::ParameterId id) { return symbol(id).tensor(); }
        unsigned int dimension() const { return type_.dimension(); }
        bool axisymmetric() const { return type_.options().axisymmetric; }

    private:
        ParameterValue compute(registry::ParameterId id);
        ParameterValue read_cell(const registry::ParameterInfo &info) const;
        [[noreturn]] void undefined(const registry::ParameterInfo &info, const char *what, unsigned int id, bool row_exists) const;

        const ElementType &type_;
        Inputs inputs_;
        std::vector<ParameterValue> values_;
        std::vector<char> computed_;
    };
} // namespace pytec3d
