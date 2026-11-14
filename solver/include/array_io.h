/**
 * @file array_io.h
 * @brief The field arrays the model needs (model_input.h): the fields of its element types and their parameters.
 */
#pragma once

#include <string>
#include <vector>

#include "model_input.h"

class ModelElementTypes;

namespace pytec3d
{
    /** @brief A field array the element types of a model read or write. */
    struct FieldArrayRequirement
    {
        std::string name;
        bool point = true;            ///< Point field, else cell field.
        std::size_t n_components = 1; ///< 3 for a vector field, 1 for a scalar.
    };

    /**
     * @brief The field arrays of the element types of the model, each once (the same as the Python
     * pytec3d.element_library.grid_arrays()): for each DOF field F the point fields F and _F and its nodal load; the
     * nodal fields (point); the element fields (cell); the point and cell parameters (mn., me.) of the element types.
     */
    std::vector<FieldArrayRequirement> required_field_arrays(const ModelElementTypes &element_types);

    /**
     * @brief Checks that @p arrays has every field of @p required, with n_points / n_cells rows and the components of
     * the field; aborts listing the problems. With @p warn_unused, warns about the fields @p arrays has and the model
     * does not use.
     */
    void check_field_arrays(const ModelArrays &arrays, const std::vector<FieldArrayRequirement> &required,
                            std::size_t n_points, std::size_t n_cells, bool warn_unused);

} // namespace pytec3d
