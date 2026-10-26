/**
 * @file square_conductor_2d.cpp
 * @brief The pyTEC3D solver driven from C++, without Python and without VTK.
 *
 * A planar 2D magnetostatic model built in memory: a 2 m x 2 m air box meshed with N x N QUAD4 cells, with a
 * 0.4 m x 0.4 m stranded conductor in the middle carrying Js = 1 MA/m^2 along z, and A_z = 0 on the outer boundary.
 * The program owns every array (std::vector); the solver receives views of them (pytec3d::ArrayView, model_input.h),
 * reads the loads and the Dirichlet values from them and writes the solution into A, B, H and Jt. Doubling Js in
 * place and solving again doubles A: the solver reads the arrays at every solve.
 *
 * Build: cmake -DPYTEC3D_BUILD_CPP_EXAMPLES=ON (build.sh: PYTEC3D_CMAKE_ARGS=-DPYTEC3D_BUILD_CPP_EXAMPLES=ON).
 * Run:   ./square_conductor_2d [N]        or  mpirun -n 2 ./square_conductor_2d [N]
 * The element matrices are stored in ./emat, as for every pyTEC3D run.
 */
#include "pytec3d_session.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace
{
    constexpr double kHalfBox = 1.0;       ///< The box is [-1, 1] x [-1, 1] m.
    constexpr double kHalfConductor = 0.2; ///< The conductor is [-0.2, 0.2] x [-0.2, 0.2] m.
    constexpr double kCurrentDensity = 1e6; ///< A/m^2, along z.

    /** @brief The model data: two element types (air, conductor) of AV2D_nodal, one material, one constants row. */
    nlohmann::json model_data()
    {
        const nlohmann::json identity = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
        return {
            {"element_types",
             {{"1", {{"element_type", "AV2D_nodal"}, {"option1", "NON_CONDUCTING"}, {"option2", 1}, {"option3", 0}}},
              {"2", {{"element_type", "AV2D_nodal"}, {"option1", "STRANDED_CONDUCTOR"}, {"option2", 1}, {"option3", 0}}}}},
            {"materials", {{"1", {{"mur", identity}}}}},
            {"constants", {{"1", {{"thickness", 1.0}}}}},
            {"dof_couplings", nlohmann::json::object()},
            {"nonlinear_laws", nlohmann::json::object()},
            {"coordinate_systems", nlohmann::json::object()},
        };
    }
} // namespace

int main(int argc, char **argv)
{
    const int n = argc > 1 ? std::max(2, std::atoi(argv[1])) : 40; // cells per side
    const std::size_t n_points = static_cast<std::size_t>((n + 1) * (n + 1));
    const std::size_t n_cells = static_cast<std::size_t>(n * n);
    const double h = 2.0 * kHalfBox / n;
    const double nan = std::numeric_limits<double>::quiet_NaN();

    // --- The mesh: points and QUAD4 cells (VTK node order: counterclockwise) --------------------------------------
    std::vector<double> points(3 * n_points);
    for (int j = 0; j <= n; ++j)
        for (int i = 0; i <= n; ++i)
        {
            const std::size_t p = static_cast<std::size_t>(j * (n + 1) + i);
            points[3 * p] = -kHalfBox + i * h;
            points[3 * p + 1] = -kHalfBox + j * h;
            points[3 * p + 2] = 0.0;
        }
    std::vector<std::int64_t> offsets(n_cells + 1);
    std::vector<std::int64_t> connectivity;
    std::vector<std::uint8_t> cell_types(n_cells, pytec3d::cell_type::QUAD);
    connectivity.reserve(4 * n_cells);
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
        {
            const std::int64_t p = j * (n + 1) + i;
            offsets[static_cast<std::size_t>(j * n + i)] = static_cast<std::int64_t>(connectivity.size());
            connectivity.insert(connectivity.end(), {p, p + 1, p + n + 2, p + n + 1});
        }
    offsets[n_cells] = static_cast<std::int64_t>(connectivity.size());

    // --- The cell arrays: ids, and the load Js (3 components: a vector field) -------------------------------------
    std::vector<std::int32_t> type(n_cells, 1), material(n_cells, 1), constants(n_cells, 1), coordinate_system(n_cells, 0);
    std::vector<double> js(3 * n_cells, 0.0);
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
        {
            const double x = -kHalfBox + (i + 0.5) * h;
            const double y = -kHalfBox + (j + 0.5) * h;
            if (std::abs(x) < kHalfConductor && std::abs(y) < kHalfConductor)
            {
                const std::size_t c = static_cast<std::size_t>(j * n + i);
                type[c] = 2;
                js[3 * c + 2] = kCurrentDensity;
            }
        }

    // --- The fields of the element type AV2D_nodal (pytec3d::required_field_arrays(); the solver checks they are all
    //     given): DOF fields F with their Dirichlet values _F (NaN: free), nodal loads, nodal and element fields. Only
    //     A, its boundary values and Js matter in this model; V, VOLT, EMF belong to the conductor domains not used.
    std::map<std::string, std::vector<double>> point_storage, cell_storage;
    for (const char *name : {"A", "V", "VOLT", "EMF"})
    {
        const std::size_t components = std::string(name) == "A" ? 3 : 1;
        point_storage[name].assign(components * n_points, 0.0);
        point_storage[std::string("_") + name].assign(components * n_points, nan);
    }
    for (const char *name : {"AMPS", "CURR"})
        point_storage[name].assign(n_points, 0.0);
    for (const char *name : {"B", "H"})
        point_storage[name].assign(3 * n_points, 0.0);
    for (const char *name : {"B", "H", "Jt"})
        cell_storage[name].assign(3 * n_cells, 0.0);
    cell_storage["Js"] = js;

    // A_z = 0 on the outer boundary (component 2 of _A).
    std::vector<double> &a_fixed = point_storage["_A"];
    for (int j = 0; j <= n; ++j)
        for (int i = 0; i <= n; ++i)
            if (i == 0 || j == 0 || i == n || j == n)
                a_fixed[3 * static_cast<std::size_t>(j * (n + 1) + i) + 2] = 0.0;

    // --- Views of all of them (typed: double fields and points, int32 ids, int64 mesh indices) ----------------------
    pytec3d::ModelInput input;
    input.mesh.points = pytec3d::ConstRealArray(points.data(), n_points, 3);
    input.mesh.offsets = pytec3d::IndexArray(offsets.data(), offsets.size());
    input.mesh.connectivity = pytec3d::IndexArray(connectivity.data(), connectivity.size());
    input.mesh.cell_types = pytec3d::CellTypeArray(cell_types.data(), n_cells);
    input.arrays.type = pytec3d::IdArray(type.data(), n_cells);
    input.arrays.material = pytec3d::IdArray(material.data(), n_cells);
    input.arrays.constant = pytec3d::IdArray(constants.data(), n_cells);
    input.arrays.coordinate_system = pytec3d::IdArray(coordinate_system.data(), n_cells);
    for (auto &[name, values] : point_storage)
        input.arrays.point_fields.emplace(name, pytec3d::RealArray(values.data(), n_points, values.size() / n_points));
    for (auto &[name, values] : cell_storage)
        input.arrays.cell_fields.emplace(name, pytec3d::RealArray(values.data(), n_cells, values.size() / n_cells));
    input.model_data = model_data();

    TimeStepInput time_step;
    time_step.time = 1.0;
    time_step.transient_effects = false; // magnetostatic
    SolverSettings settings;
    settings.petsc_options.clear(); // no PETSc monitors

    int rank = 0;
    {
        pyTEC3DSession session(1);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        session.load_model(input);

        const std::vector<double> &a = point_storage["A"];
        const std::vector<double> &b_nodes = point_storage["B"];
        std::vector<double> &js_view = cell_storage["Js"]; // the memory the solver views

        session.solve(time_step, settings);
        const std::size_t centre = static_cast<std::size_t>((n / 2) * (n + 1) + n / 2);
        const double a_centre = a[3 * centre + 2];
        double b_max = 0.0;
        for (std::size_t p = 0; p < n_points; ++p)
            b_max = std::max(b_max, std::hypot(b_nodes[3 * p], b_nodes[3 * p + 1]));

        // The loads are read at every solve: doubling Js in place doubles the solution.
        for (std::size_t c = 0; c < n_cells; ++c)
            js_view[3 * c + 2] *= 2.0;
        session.solve(time_step, settings);
        const double ratio = a[3 * centre + 2] / a_centre;

        if (rank == 0)
            std::printf("square_conductor_2d: %d x %d cells, I = %.3e A\n"
                        "  A_z at the centre      %.6e Wb/m\n"
                        "  max |B| at the nodes   %.6e T\n"
                        "  A_z ratio with 2 Js    %.12f\n",
                        n, n, kCurrentDensity * 4 * kHalfConductor * kHalfConductor, a_centre, b_max, ratio);
        session.close();
    }
    pyTEC3DSession::shutdown_runtime();
    return 0;
}
