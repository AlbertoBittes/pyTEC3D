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
 * A nonlinear law (BH) stands in a kernel where a parameter would. On an element without a BH curve it is
 * its linear parameter (a tensor); on an element with one it is applied to the trial value by the
 * nonlinear kernels, through apply_law().
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

class BHCurve;
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
            const BHCurve *bh_curve = nullptr; ///< Set for nonlinear elements.
        };

        ElementParameterValues(const ElementType &type, const Inputs &inputs);

        /** @brief Value of @p reference with its sign; a law gives its linear parameter. */
        ParameterValue value(const ParameterRef &reference);

        /** @brief True when @p reference is a law and the element has a curve for it. */
        bool nonlinear(const ParameterRef &reference) const;

        /** @brief H and dH/dB of the law of @p reference at B, with its sign. Only for nonlinear(). */
        void apply_law(const ParameterRef &reference, const libMesh::Real B[3], libMesh::Real H[3], libMesh::Real J[3][3]) const;

        const ParameterValue &symbol(registry::ParameterId id);

        // What registry::evaluate_expression() reads.
        const libMesh::Real &get_scalar(registry::ParameterId id) { return symbol(id).scalar; }
        const libMesh::VectorValue<libMesh::Real> &get_vector(registry::ParameterId id) { return symbol(id).vector; }
        const libMesh::TensorValue<libMesh::Real> &get_tensor(registry::ParameterId id) { return symbol(id).tensor; }
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
