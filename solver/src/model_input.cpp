/**
 * @file model_input.cpp
 * @brief Checks of the views handed to the solver (model_input.h).
 */
#include "model_input.h"

#include "logging.h"

#include <string>

namespace pytec3d
{
    namespace
    {
        [[noreturn]] void fail(const std::string &message)
        {
            logging::error_and_abort("Model input: " + message);
        }

        template <typename T>
        std::string shape(const ArrayView<T> &view)
        {
            return "(" + std::to_string(view.n_tuples()) + ", " + std::to_string(view.n_components()) + ")";
        }

        /** @brief @p view is given, with @p n_tuples rows of @p n_components values. */
        template <typename T>
        void check_shape(const ArrayView<T> &view, std::size_t n_tuples, std::size_t n_components, const std::string &what)
        {
            if (view.empty())
                fail(what + " is missing.");
            if (view.n_tuples() != n_tuples || view.n_components() != n_components)
                fail(what + " has shape " + shape(view) + ", expected (" + std::to_string(n_tuples) + ", " +
                     std::to_string(n_components) + ").");
        }

        void check_rotations(const ConstRealArray &view, std::size_t n_tuples, const std::string &what)
        {
            if (view.empty())
                return;
            if (view.n_tuples() != n_tuples || view.n_components() != 9)
                fail(what + " has shape " + shape(view) + ", expected (" + std::to_string(n_tuples) + ", 9), or (" +
                     std::to_string(n_tuples) + ", 3, 3).");
        }

        const RealArray *find(const std::map<std::string, RealArray> &fields, const std::string &name)
        {
            const auto it = fields.find(name);
            return it == fields.end() || it->second.empty() ? nullptr : &it->second;
        }
    } // namespace

    const RealArray *ModelArrays::point_field(const std::string &name) const
    {
        return find(point_fields, name);
    }

    const RealArray *ModelArrays::cell_field(const std::string &name) const
    {
        return find(cell_fields, name);
    }

    void ModelArrays::validate(std::size_t n_points, std::size_t n_cells) const
    {
        check_shape(type, n_cells, 1, "the cell array TYPE");
        check_shape(material, n_cells, 1, "the cell array MATERIAL");
        check_shape(constant, n_cells, 1, "the cell array CONSTANT");
        check_shape(coordinate_system, n_cells, 1, "the cell array COORDINATE_SYSTEM");
        for (const auto &[name, view] : couplings)
            check_shape(view, n_points, 1, "the coupling array '" + name + "'");
    }

    void ModelInput::validate() const
    {
        const std::size_t n_points = mesh.n_points();
        const std::size_t n_cells = mesh.n_cells();
        check_shape(mesh.points, n_points, 3, "the points");
        check_shape(mesh.offsets, n_cells + 1, 1, "the cell offsets");
        const std::int64_t n_connectivity = mesh.offsets(n_cells);
        if (n_connectivity < 0 || mesh.connectivity.n_tuples() * mesh.connectivity.n_components() != static_cast<std::size_t>(n_connectivity))
            fail("the connectivity has " + std::to_string(mesh.connectivity.n_tuples()) + " entries, the last offset is " +
                 std::to_string(n_connectivity) + ".");
        arrays.validate(n_points, n_cells);
        check_rotations(element_rotations, n_cells, "element_rotations");
        check_rotations(node_rotations, n_points, "node_rotations");
        if (!model_data.is_object())
            fail("the model data must be a JSON object.");
    }
} // namespace pytec3d
