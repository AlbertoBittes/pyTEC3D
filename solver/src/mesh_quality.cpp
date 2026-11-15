/**
 * @file mesh_quality.cpp
 * @brief Edge ratio, angles and scaled Jacobian of the mesh elements.
 */
#include "mesh_quality.h"

#include <array>
#include <cmath>

#include "fe_mesh.h"
#include "quadrature.h"

#include "logging.h"


#include <algorithm>
#include <limits>
#include <map>
#include <memory>
#include <utility>

namespace pytec3d
{
    namespace
    {
        constexpr std::size_t kTypes = 80; ///< VTK cell type numbers are below this (Lagrange cells: 69 - 73).

        using P = std::array<double, 3>;
        P minus(const double *a, const double *b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
        double length(const P &v) { return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); }
        double angle(const P &a, const P &b)
        {
            const double c = (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) / (length(a) * length(b));
            return std::acos(std::max(-1.0, std::min(1.0, c))) * 180.0 / 3.14159265358979323846;
        }
    } // namespace

    std::vector<ShapeQuality> compute_mesh_quality(const fem::FEMesh &mesh)
    {
        // Indexed by VTK cell type, so that the reductions line up on every rank.
        const double huge = std::numeric_limits<double>::max();
        std::vector<double> max_edge_ratio(kTypes, 0.0), min_angle(kTypes, huge), max_angle(kTypes, 0.0), min_scaled(kTypes, huge);
        std::vector<unsigned long long> count(kTypes, 0), distorted(kTypes, 0), degenerate(kTypes, 0);
        std::vector<long long> worst(kTypes, 0);

        // det(J) at the points of a rule of twice the geometric order, for all the elements at once.
        const mfem::ParMesh &pmesh = mesh.mesh();
        const std::vector<int> orders(static_cast<std::size_t>(mesh.n_local_elements()), 2 * mesh.geometric_order());
        const fem::QuadratureData quadrature(mesh, orders);
        mfem::Array<int> vertices, edges, orientations;
        for (int e = 0; e < mesh.n_local_elements(); ++e)
        {
            const auto t = static_cast<std::size_t>(mesh.cell_type(e));
            ++count[t];
            // Edge length ratio and corner angles, from the corners.
            pmesh.GetElementVertices(e, vertices);
            const mfem::Element *element = pmesh.GetElement(e);
            const int n_edges = element->GetNEdges();
            double shortest = huge, longest = 0.0;
            for (int k = 0; k < n_edges; ++k)
            {
                const int *ev = element->GetEdgeVertices(k);
                const double l = length(minus(pmesh.GetVertex(vertices[ev[1]]), pmesh.GetVertex(vertices[ev[0]])));
                shortest = std::min(shortest, l);
                longest = std::max(longest, l);
            }
            max_edge_ratio[t] = std::max(max_edge_ratio[t], shortest > 0.0 ? longest / shortest : huge);
            // The angles at the corners of the faces (2D: of the element).
            const auto corner_angles = [&](const int *face, int n)
            {
                for (int k = 0; k < n; ++k)
                {
                    const double *v = pmesh.GetVertex(vertices[face[k]]);
                    const double a = angle(minus(pmesh.GetVertex(vertices[face[(k + 1) % n]]), v),
                                           minus(pmesh.GetVertex(vertices[face[(k + n - 1) % n]]), v));
                    min_angle[t] = std::min(min_angle[t], a);
                    max_angle[t] = std::max(max_angle[t], a);
                }
            };
            if (mesh.dimension() == 2)
            {
                std::vector<int> all(static_cast<std::size_t>(vertices.Size()));
                for (int k = 0; k < vertices.Size(); ++k)
                    all[static_cast<std::size_t>(k)] = k;
                corner_angles(all.data(), vertices.Size());
            }
            else
                for (int f = 0; f < element->GetNFaces(); ++f)
                {
                    const int *face = element->GetFaceVertices(f);
                    const int n = element->GetNFaceVertices(f);
                    corner_angles(face, n);
                }

            const int nq = quadrature.n_points(e);
            const double *det = quadrature.det(e);
            double lowest = huge, highest = -huge;
            for (int q = 0; q < nq; ++q)
            {
                lowest = std::min(lowest, det[q]);
                highest = std::max(highest, det[q]);
            }
            const double scaled = highest > 0.0 ? lowest / highest : 0.0;
            if (scaled <= 0.0)
                ++degenerate[t];
            else if (scaled < kDistortedScaledJacobian)
                ++distorted[t];
            if (scaled < min_scaled[t])
            {
                min_scaled[t] = scaled;
                worst[t] = mesh.cell(e);
            }
        }

        const MPI_Comm comm = mesh.comm();
        MPI_Allreduce(MPI_IN_PLACE, count.data(), static_cast<int>(kTypes), MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
        MPI_Allreduce(MPI_IN_PLACE, distorted.data(), static_cast<int>(kTypes), MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
        MPI_Allreduce(MPI_IN_PLACE, degenerate.data(), static_cast<int>(kTypes), MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
        MPI_Allreduce(MPI_IN_PLACE, max_edge_ratio.data(), static_cast<int>(kTypes), MPI_DOUBLE, MPI_MAX, comm);
        MPI_Allreduce(MPI_IN_PLACE, min_angle.data(), static_cast<int>(kTypes), MPI_DOUBLE, MPI_MIN, comm);
        MPI_Allreduce(MPI_IN_PLACE, max_angle.data(), static_cast<int>(kTypes), MPI_DOUBLE, MPI_MAX, comm);

        std::vector<ShapeQuality> out;
        for (std::size_t t = 0; t < kTypes; ++t)
        {
            if (count[t] == 0)
                continue;
            // The worst element of all the ranks (the smallest scaled Jacobian, the lowest rank on ties).
            struct
            {
                double value;
                int rank;
            } local{min_scaled[t], mesh.rank()}, global{};
            MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
            long long worst_id = worst[t];
            MPI_Bcast(&worst_id, 1, MPI_LONG_LONG, global.rank, comm);

            ShapeQuality quality;
            quality.shape = fem::cell_type_name(static_cast<int>(t));
            quality.elements = static_cast<std::size_t>(count[t]);
            quality.max_edge_ratio = max_edge_ratio[t];
            quality.min_angle = min_angle[t] == huge ? 0.0 : min_angle[t];
            quality.max_angle = max_angle[t];
            quality.min_scaled_jacobian = global.value;
            quality.worst_element = worst_id;
            quality.distorted = static_cast<std::size_t>(distorted[t]);
            quality.degenerate = static_cast<std::size_t>(degenerate[t]);
            out.push_back(std::move(quality));
        }
        return out;
    }

    void check_mesh_quality(const std::vector<ShapeQuality> &quality)
    {
        for (const ShapeQuality &shape : quality)
        {
            if (shape.degenerate > 0)
                logging::warn("{} {} element(s) are degenerate or inverted (det(J) <= 0 at a quadrature point; e.g. cell "
                              "{}): the solution there is meaningless.",
                              shape.degenerate, shape.shape, shape.worst_element);
            else if (shape.distorted > 0)
                logging::warn("{} {} element(s) are distorted (scaled Jacobian < {}; the worst, {:.3f}, at cell {}).",
                              shape.distorted, shape.shape, kDistortedScaledJacobian, shape.min_scaled_jacobian,
                              shape.worst_element);
        }
    }
} // namespace pytec3d
