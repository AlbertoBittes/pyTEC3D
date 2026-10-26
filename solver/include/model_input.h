/**
 * @file model_input.h
 * @brief What the solver reads and writes: typed views of the caller's arrays and the model data.
 *
 * The solver owns no model array. The caller (the Python API, or any C++ program) keeps the mesh, the model arrays
 * and the output fields in its own memory and hands the solver views of them (ArrayView): the solver reads the
 * mesh, the ids, the loads and the Dirichlet values from them, and writes the solution into the output fields, in
 * place, after each solve. Nothing is copied, so a load changed by the caller between two solves is the load of the
 * next solve.
 *
 * Every array has a fixed type: real numbers are double, ids int32, mesh indices int64, cell types uint8.
 *
 * Mesh (MeshView), with the conventions of the VTK unstructured grid (no VTK code is involved), read by load_model()
 * only (the solver copies what it needs, so it can be freed when load_model() returns):
 *  - points (n_points, 3), double;
 *  - cells: cell i has the points connectivity[offsets[i] : offsets[i + 1]] (int64), in the node order of the VTK cell
 *    of type cell_types[i] (uint8, the VTK cell type numbers, cell_type below).
 *
 * Model arrays (ModelArrays), which must stay alive and keep their address from load_model() to the end of the
 * session, or until they are replaced with set_arrays():
 *  - the cell ids TYPE, MATERIAL, CONSTANT, COORDINATE_SYSTEM (n_cells, int32): the rows of the model data;
 *  - the DOF couplings CP_<k> (n_points, int32), read when the model is loaded;
 *  - the fields (double), exactly those of the element types of the model data (required_field_arrays(), array_io.h):
 *      point: each DOF field F and its Dirichlet values _F (NaN where free), the nodal fields, the nodal loads;
 *      cell:  the element fields, the element loads;
 *    3 components for a vector field, 1 for a scalar. The loads and the Dirichlet values are read at every solve; the
 *    DOF, nodal and element fields are written at its end.
 *
 * Rotations, computed by the caller from its coordinate systems (row-major 3 x 3 matrices, local to global axes,
 * double), read by load_model() only:
 *  - element_rotations (n_cells, 9): the element coordinate system (COORDINATE_SYSTEM) of each cell at its vertex
 *    average (the mean of its corner points);
 *  - node_rotations (n_points, 9): the axes of the DOFs of each point.
 * Empty views mean identity everywhere.
 */
#pragma once

#include <cstddef>
#include <cstdint>
#include <map>
#include <string>

#include <nlohmann/json.hpp>

namespace pytec3d
{
    /**
     * @brief A view of an array owned by the caller: n_tuples rows of n_components values of type @p T, contiguous,
     * row-major. ArrayView<const T> is read-only. A default-constructed view is empty.
     */
    template <typename T>
    class ArrayView
    {
    public:
        ArrayView() = default;
        ArrayView(T *data, std::size_t n_tuples, std::size_t n_components = 1)
            : data_(data), n_tuples_(n_tuples), n_components_(n_components)
        {
        }

        bool empty() const { return data_ == nullptr; }
        std::size_t n_tuples() const { return n_tuples_; }
        std::size_t n_components() const { return n_components_; }
        T *data() const { return data_; }

        /** @brief Value @p component of tuple @p tuple. */
        T &operator()(std::size_t tuple, std::size_t component = 0) const { return data_[tuple * n_components_ + component]; }

    private:
        T *data_ = nullptr;
        std::size_t n_tuples_ = 0;
        std::size_t n_components_ = 1;
    };

    using RealArray = ArrayView<double>;              ///< A field: loads, Dirichlet values, outputs.
    using ConstRealArray = ArrayView<const double>;   ///< Points, rotations.
    using IdArray = ArrayView<const std::int32_t>;    ///< Cell ids, CP_<k> couplings.
    using IndexArray = ArrayView<const std::int64_t>; ///< Cell offsets and connectivity.
    using CellTypeArray = ArrayView<const std::uint8_t>;

    /** @brief VTK cell type numbers of the cells the solver reads (vtkCellType.h). */
    namespace cell_type
    {
        constexpr int LINE = 3;
        constexpr int TRIANGLE = 5;
        constexpr int QUAD = 9;
        constexpr int TETRA = 10;
        constexpr int HEXAHEDRON = 12;
        constexpr int WEDGE = 13;
        constexpr int PYRAMID = 14;
        constexpr int QUADRATIC_EDGE = 21;
        constexpr int QUADRATIC_TRIANGLE = 22;
        constexpr int QUADRATIC_QUAD = 23;
        constexpr int QUADRATIC_TETRA = 24;
        constexpr int QUADRATIC_HEXAHEDRON = 25;
        constexpr int QUADRATIC_WEDGE = 26;
        constexpr int QUADRATIC_PYRAMID = 27;
        constexpr int BIQUADRATIC_QUAD = 28;
        constexpr int TRIQUADRATIC_HEXAHEDRON = 29;
        constexpr int BIQUADRATIC_QUADRATIC_WEDGE = 32;
    } // namespace cell_type

    /** @brief The mesh: points and cells (see the file comment). */
    struct MeshView
    {
        ConstRealArray points;    ///< (n_points, 3).
        IndexArray offsets;       ///< (n_cells + 1): the points of cell i are connectivity[offsets[i] : offsets[i + 1]].
        IndexArray connectivity;  ///< (offsets[n_cells]) point ids.
        CellTypeArray cell_types; ///< (n_cells) VTK cell type numbers.

        std::size_t n_points() const { return points.n_tuples(); }
        std::size_t n_cells() const { return cell_types.n_tuples(); }
    };

    /** @brief The model arrays (see the file comment). */
    struct ModelArrays
    {
        IdArray type;
        IdArray material;
        IdArray constant;
        IdArray coordinate_system;
        std::map<std::string, IdArray> couplings;      ///< CP_<k> point arrays.
        std::map<std::string, RealArray> point_fields; ///< By name: F, _F, nodal fields, nodal loads.
        std::map<std::string, RealArray> cell_fields;  ///< By name: element fields, element loads.

        /** @brief The point field @p name, or nullptr. */
        const RealArray *point_field(const std::string &name) const;
        /** @brief The cell field @p name, or nullptr. */
        const RealArray *cell_field(const std::string &name) const;

        /** @brief Checks the ids and the couplings (given, one component, n_cells / n_points rows); aborts naming the array. */
        void validate(std::size_t n_points, std::size_t n_cells) const;
    };

    /** @brief Everything load_model() needs (see the file comment). */
    struct ModelInput
    {
        MeshView mesh;
        ModelArrays arrays;
        ConstRealArray element_rotations; ///< (n_cells, 9), or empty for identity.
        ConstRealArray node_rotations;    ///< (n_points, 9), or empty for identity.
        nlohmann::json model_data;        ///< The model data (model_data.h): element_types, materials, constants...

        /** @brief Checks the shapes of the mesh, the ids and the rotations; aborts naming what is wrong. */
        void validate() const;
    };
} // namespace pytec3d
