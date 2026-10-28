/**
 * @file array_io.h
 * @brief The libMesh mesh from the mesh view of the caller, the field arrays the model needs, and the solution
 * written into them.
 *
 * The mesh, the ids and the fields are the arrays of model_input.h: libMesh node i is point i, element i is cell i
 * (subdomain: its TYPE), and the field F of the solution systems is written into the point or cell field F.
 */
#pragma once

#include <string>
#include <vector>

#include "model_input.h"

class ModelElementTypes;

namespace libMesh
{
    class EquationSystems;
    class MeshBase;
}

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

    /**
     * @brief Fills @p mesh (cleared first) with the points and cells of @p view; each element's subdomain is its
     * value in @p type_ids (the TYPE cell array). Every rank reads the whole mesh.
     */
    void import_mesh(libMesh::MeshBase &mesh, const MeshView &view, const IdArray &type_ids);

    /**
     * @brief Output fields the solution system @p system_name writes and @p arrays does not have (or with too few
     * components), as "name (point)" / "name (cell)".
     */
    std::vector<std::string> missing_output_arrays(const libMesh::EquationSystems &systems, const std::string &system_name,
                                                   const ModelArrays &arrays);

    /**
     * @brief Writes the fields of the solution system @p system_name into the fields of @p arrays: the nodal systems
     * into the point fields, the element system into the cell fields. Collective; every rank writes every value.
     */
    void write_system_to_arrays(libMesh::EquationSystems &systems, const std::string &system_name, const ModelArrays &arrays);
} // namespace pytec3d
