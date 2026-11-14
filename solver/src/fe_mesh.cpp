/**
 * @file fe_mesh.cpp
 * @brief FEMesh: import of the caller's cells into an MFEM ParMesh (fe_mesh.h).
 */
#include "fe_mesh.h"

#include <algorithm>
#include <array>
#include <map>
#include <set>
#include <string>
#include <unordered_map>

#include "fe_collections.h"
#include "logging.h"

namespace pytec3d::fem
{
    namespace
    {
        [[noreturn]] void fail(const std::string &message)
        {
            logging::error_and_abort("Mesh import: " + message);
        }

        /** @brief A non-corner point of a VTK cell: the corners it is the midpoint (edge), centre (face, cell) of. */
        struct HighOrderNode
        {
            int node;
            std::vector<int> corners;
        };

        /** @brief What the solver needs of a VTK cell type. */
        struct CellSpec
        {
            int vtk = 0;
            mfem::Geometry::Type geometry = mfem::Geometry::INVALID;
            int dimension = 0;
            int n_corners = 0;
            int n_nodes = 0;
            int order = 1;
            int quad_faces = 0;                   ///< 0: none or not quadratic, 1: serendipity (no centre), 2: with centre.
            std::vector<HighOrderNode> high;      ///< Edge midpoints, face and cell centres.
            std::vector<int> flip;                ///< Corner permutation that reverses the orientation.
        };

        const CellSpec &cell_spec(int vtk, std::size_t cell)
        {
            using G = mfem::Geometry;
            static const std::map<int, CellSpec> specs = []
            {
                std::map<int, CellSpec> s;
                const std::vector<HighOrderNode> tri6 = {{3, {0, 1}}, {4, {1, 2}}, {5, {2, 0}}};
                const std::vector<HighOrderNode> quad8 = {{4, {0, 1}}, {5, {1, 2}}, {6, {2, 3}}, {7, {3, 0}}};
                const std::vector<HighOrderNode> tet10 = {{4, {0, 1}}, {5, {1, 2}}, {6, {2, 0}},
                                                          {7, {0, 3}}, {8, {1, 3}}, {9, {2, 3}}};
                const std::vector<HighOrderNode> hex20 = {
                    {8, {0, 1}},  {9, {1, 2}},  {10, {2, 3}}, {11, {3, 0}}, {12, {4, 5}}, {13, {5, 6}},
                    {14, {6, 7}}, {15, {7, 4}}, {16, {0, 4}}, {17, {1, 5}}, {18, {2, 6}}, {19, {3, 7}}};
                const std::vector<HighOrderNode> wedge15 = {{6, {0, 1}},  {7, {1, 2}},  {8, {2, 0}},
                                                            {9, {3, 4}},  {10, {4, 5}}, {11, {5, 3}},
                                                            {12, {0, 3}}, {13, {1, 4}}, {14, {2, 5}}};
                const std::vector<HighOrderNode> pyramid13 = {{5, {0, 1}},  {6, {1, 2}},  {7, {2, 3}},  {8, {3, 0}},
                                                              {9, {0, 4}},  {10, {1, 4}}, {11, {2, 4}}, {12, {3, 4}}};
                auto join = [](std::vector<HighOrderNode> a, const std::vector<HighOrderNode> &b)
                {
                    a.insert(a.end(), b.begin(), b.end());
                    return a;
                };
                const std::vector<int> tri_flip = {0, 2, 1}, quad_flip = {0, 3, 2, 1}, tet_flip = {1, 0, 2, 3},
                                       hex_flip = {4, 5, 6, 7, 0, 1, 2, 3}, wedge_flip = {3, 4, 5, 0, 1, 2},
                                       pyramid_flip = {0, 3, 2, 1, 4};
                s[cell_type::TRIANGLE] = {cell_type::TRIANGLE, G::TRIANGLE, 2, 3, 3, 1, 0, {}, tri_flip};
                s[cell_type::QUADRATIC_TRIANGLE] = {cell_type::QUADRATIC_TRIANGLE, G::TRIANGLE, 2, 3, 6, 2, 0, tri6, tri_flip};
                s[cell_type::QUAD] = {cell_type::QUAD, G::SQUARE, 2, 4, 4, 1, 0, {}, quad_flip};
                s[cell_type::QUADRATIC_QUAD] = {cell_type::QUADRATIC_QUAD, G::SQUARE, 2, 4, 8, 2, 1, quad8, quad_flip};
                s[cell_type::BIQUADRATIC_QUAD] = {cell_type::BIQUADRATIC_QUAD, G::SQUARE, 2, 4, 9, 2, 2,
                                                  join(quad8, {{8, {0, 1, 2, 3}}}), quad_flip};
                s[cell_type::TETRA] = {cell_type::TETRA, G::TETRAHEDRON, 3, 4, 4, 1, 0, {}, tet_flip};
                s[cell_type::QUADRATIC_TETRA] = {cell_type::QUADRATIC_TETRA, G::TETRAHEDRON, 3, 4, 10, 2, 0, tet10, tet_flip};
                s[cell_type::HEXAHEDRON] = {cell_type::HEXAHEDRON, G::CUBE, 3, 8, 8, 1, 0, {}, hex_flip};
                s[cell_type::QUADRATIC_HEXAHEDRON] = {cell_type::QUADRATIC_HEXAHEDRON, G::CUBE, 3, 8, 20, 2, 1, hex20, hex_flip};
                // VTK_TRIQUADRATIC_HEXAHEDRON: face centres x-, x+, y-, y+, z-, z+, then the cell centre.
                s[cell_type::TRIQUADRATIC_HEXAHEDRON] = {
                    cell_type::TRIQUADRATIC_HEXAHEDRON, G::CUBE, 3, 8, 27, 2, 2,
                    join(hex20, {{20, {0, 3, 7, 4}}, {21, {1, 2, 6, 5}}, {22, {0, 1, 5, 4}}, {23, {3, 2, 6, 7}},
                                 {24, {0, 1, 2, 3}}, {25, {4, 5, 6, 7}}, {26, {0, 1, 2, 3, 4, 5, 6, 7}}}),
                    hex_flip};
                s[cell_type::WEDGE] = {cell_type::WEDGE, G::PRISM, 3, 6, 6, 1, 0, {}, wedge_flip};
                s[cell_type::QUADRATIC_WEDGE] = {cell_type::QUADRATIC_WEDGE, G::PRISM, 3, 6, 15, 2, 1, wedge15, wedge_flip};
                // VTK_BIQUADRATIC_QUADRATIC_WEDGE: centres of the quadrilateral faces 0-1-4-3, 1-2-5-4, 2-0-3-5.
                s[cell_type::BIQUADRATIC_QUADRATIC_WEDGE] = {
                    cell_type::BIQUADRATIC_QUADRATIC_WEDGE, G::PRISM, 3, 6, 18, 2, 2,
                    join(wedge15, {{15, {0, 1, 4, 3}}, {16, {1, 2, 5, 4}}, {17, {2, 0, 3, 5}}}), wedge_flip};
                s[cell_type::PYRAMID] = {cell_type::PYRAMID, G::PYRAMID, 3, 5, 5, 1, 0, {}, pyramid_flip};
                s[cell_type::QUADRATIC_PYRAMID] = {cell_type::QUADRATIC_PYRAMID, G::PYRAMID, 3, 5, 13, 2, 1, pyramid13, pyramid_flip};
                return s;
            }();
            const auto it = specs.find(vtk);
            if (it == specs.end())
                fail("cell " + std::to_string(cell) + " has the VTK cell type " + std::to_string(vtk) +
                     ", which the solver does not support (2D: triangles and quadrilaterals, 3D: tetrahedra, "
                     "hexahedra, wedges and pyramids, linear or quadratic).");
            return it->second;
        }

        using P3 = std::array<double, 3>;

        P3 sub(const P3 &a, const P3 &b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
        double triple(const P3 &a, const P3 &b, const P3 &c)
        {
            return a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
        }

        /** @brief libMesh's Elem::is_flipped(): the sign of a corner triple product (2D: of the z of a cross product). */
        bool is_flipped(mfem::Geometry::Type geometry, const std::vector<P3> &x)
        {
            using G = mfem::Geometry;
            const P3 ez = {0.0, 0.0, 1.0};
            switch (geometry)
            {
            case G::TRIANGLE: return triple(sub(x[1], x[0]), sub(x[2], x[0]), ez) < 0.0;
            case G::SQUARE: return triple(sub(x[1], x[0]), sub(x[3], x[0]), ez) < 0.0;
            case G::TETRAHEDRON: return triple(sub(x[1], x[0]), sub(x[2], x[0]), sub(x[3], x[0])) < 0.0;
            case G::CUBE: return triple(sub(x[1], x[0]), sub(x[3], x[0]), sub(x[4], x[0])) < 0.0;
            case G::PRISM: return triple(sub(x[1], x[0]), sub(x[2], x[0]), sub(x[3], x[0])) < 0.0;
            case G::PYRAMID: return triple(sub(x[1], x[0]), sub(x[3], x[0]), sub(x[4], x[0])) < 0.0;
            default: return false;
            }
        }

        /** @brief A set of points as a key: sorted ids. */
        std::vector<std::int64_t> key_of(std::vector<std::int64_t> ids)
        {
            std::sort(ids.begin(), ids.end());
            return ids;
        }

        struct KeyHash
        {
            std::size_t operator()(const std::vector<std::int64_t> &k) const
            {
                std::size_t h = 1469598103934665603ull;
                for (const auto v : k)
                    h = (h ^ static_cast<std::size_t>(v)) * 1099511628211ull;
                return h;
            }
        };
    } // namespace

    /** @brief The name of a VTK cell type, as the solver reports it. */
    const char *cell_type_name(int vtk)
    {
        switch (vtk)
        {
        case cell_type::TRIANGLE: return "TRI3";
        case cell_type::QUADRATIC_TRIANGLE: return "TRI6";
        case cell_type::QUAD: return "QUAD4";
        case cell_type::QUADRATIC_QUAD: return "QUAD8";
        case cell_type::BIQUADRATIC_QUAD: return "QUAD9";
        case cell_type::TETRA: return "TET4";
        case cell_type::QUADRATIC_TETRA: return "TET10";
        case cell_type::HEXAHEDRON: return "HEX8";
        case cell_type::QUADRATIC_HEXAHEDRON: return "HEX20";
        case cell_type::TRIQUADRATIC_HEXAHEDRON: return "HEX27";
        case cell_type::WEDGE: return "PRISM6";
        case cell_type::QUADRATIC_WEDGE: return "PRISM15";
        case cell_type::BIQUADRATIC_QUADRATIC_WEDGE: return "PRISM18";
        case cell_type::PYRAMID: return "PYRAMID5";
        case cell_type::QUADRATIC_PYRAMID: return "PYRAMID13";
        default: return "UNKNOWN";
        }
    }

    FEMesh::FEMesh(MPI_Comm comm, const MeshView &view, const IdArray &type_ids) : comm_(comm)
    {
        MPI_Comm_rank(comm, &rank_);
        MPI_Comm_size(comm, &n_ranks_);
        n_points_ = view.n_points();
        n_cells_ = view.n_cells();
        if (n_cells_ == 0)
            fail("the mesh has no cell.");

        // ---- Cells: types, orders, orientation --------------------------------------------------------------------
        std::map<int, std::size_t> type_counts;
        int order = 0;
        int quad_faces = 0;
        std::vector<std::int64_t> conn(view.connectivity.n_tuples());
        std::vector<char> is_corner(n_points_, 0);
        for (std::size_t i = 0; i < n_cells_; ++i)
        {
            const CellSpec &spec = cell_spec(view.cell_types(i), i);
            ++type_counts[spec.vtk];
            const std::int64_t begin = view.offsets(i);
            const std::int64_t n = view.offsets(i + 1) - begin;
            if (n != spec.n_nodes)
                fail("cell " + std::to_string(i) + " has " + std::to_string(n) + " points, its cell type " +
                     std::to_string(spec.n_nodes) + ".");
            if (dimension_ == 0)
                dimension_ = spec.dimension;
            else if (dimension_ != spec.dimension)
                fail("the cells have different dimensions (cell " + std::to_string(i) +
                     "): a mesh has 2D or 3D cells only.");
            if (order == 0)
                order = spec.order;
            else if (order != spec.order)
                fail("linear and quadratic cells in one mesh (cell " + std::to_string(i) + "): all the cells must have "
                     "the same geometric order.");
            if (spec.quad_faces != 0)
            {
                if (quad_faces != 0 && quad_faces != spec.quad_faces)
                    fail("quadratic cells with and without face-centre points in one mesh (cell " + std::to_string(i) +
                         "): e.g. HEX27 with HEX20, or PRISM18 with PYRAMID13.");
                quad_faces = spec.quad_faces;
            }

            std::vector<P3> x(static_cast<std::size_t>(spec.n_corners));
            for (int j = 0; j < spec.n_nodes; ++j)
            {
                const std::int64_t point = view.connectivity(static_cast<std::size_t>(begin + j));
                if (point < 0 || static_cast<std::size_t>(point) >= n_points_)
                    fail("cell " + std::to_string(i) + " refers to the point " + std::to_string(point) + ", out of range.");
                if (j < spec.n_corners)
                    x[static_cast<std::size_t>(j)] = {view.points(static_cast<std::size_t>(point), 0),
                                                      view.points(static_cast<std::size_t>(point), 1),
                                                      view.points(static_cast<std::size_t>(point), 2)};
            }
            const bool flipped = is_flipped(spec.geometry, x);
            flipped_elements_ += flipped ? 1 : 0;

            std::int64_t *out = conn.data() + begin;
            const std::int64_t *in = view.connectivity.data() + begin;
            if (!flipped)
                std::copy(in, in + n, out);
            else
            {
                // Corners permuted; each other point goes where the permuted corners of its edge, face or cell are.
                for (int j = 0; j < spec.n_corners; ++j)
                    out[j] = in[spec.flip[static_cast<std::size_t>(j)]];
                for (const HighOrderNode &target : spec.high)
                {
                    std::set<int> image;
                    for (const int c : target.corners)
                        image.insert(spec.flip[static_cast<std::size_t>(c)]);
                    for (const HighOrderNode &source : spec.high)
                        if (std::set<int>(source.corners.begin(), source.corners.end()) == image)
                            out[target.node] = in[source.node];
                }
            }
            for (int j = 0; j < spec.n_corners; ++j)
                is_corner[static_cast<std::size_t>(out[j])] = 1;
        }
        geometric_order_ = order;
        serendipity_ = quad_faces == 1;
        for (const auto &[type, count] : type_counts)
            cell_type_counts_.emplace_back(type, count);

        coordinates_.assign(view.points.data(), view.points.data() + 3 * n_points_);
        // The edges of every cell by their midpoints: the couplings of edge DOFs read them on any rank.
        for (std::size_t i = 0; i < n_cells_; ++i)
        {
            const CellSpec &spec = cell_spec(view.cell_types(i), i);
            const std::int64_t *pts = conn.data() + view.offsets(i);
            for (const HighOrderNode &node : spec.high)
                if (node.corners.size() == 2)
                    edge_corners_.emplace(pts[node.node], std::make_pair(pts[node.corners[0]], pts[node.corners[1]]));
        }

        if (dimension_ == 2)
            for (std::size_t p = 0; p < n_points_; ++p)
                if (view.points(p, 2) != 0.0)
                    fail("the 2D mesh has the point " + std::to_string(p) + " out of the xy plane (z = " +
                         std::to_string(view.points(p, 2)) + "): 2D models lie in the xy plane.");

        // ---- The serial mesh of the corners --------------------------------------------------------------------------
        std::vector<int> vertex_of_point(n_points_, -1);
        std::vector<std::int64_t> point_of_vertex;
        for (std::size_t p = 0; p < n_points_; ++p)
            if (is_corner[p])
            {
                vertex_of_point[p] = static_cast<int>(point_of_vertex.size());
                point_of_vertex.push_back(static_cast<std::int64_t>(p));
            }

        // Attributes: 1 + the rank of the TYPE id among those of the mesh (MFEM attributes are positive).
        std::set<int> types;
        for (std::size_t i = 0; i < n_cells_; ++i)
            types.insert(type_ids(i));
        subdomain_ids_ = types;
        std::map<int, int> attribute_of_type;
        for (const int t : types)
            attribute_of_type.emplace(t, static_cast<int>(attribute_of_type.size()) + 1);

        auto serial = std::make_unique<mfem::Mesh>(dimension_, static_cast<int>(point_of_vertex.size()),
                                                   static_cast<int>(n_cells_), 0, dimension_);
        for (const std::int64_t p : point_of_vertex)
        {
            const double xyz[3] = {view.points(static_cast<std::size_t>(p), 0), view.points(static_cast<std::size_t>(p), 1),
                                   view.points(static_cast<std::size_t>(p), 2)};
            serial->AddVertex(xyz);
        }
        for (std::size_t i = 0; i < n_cells_; ++i)
        {
            const CellSpec &spec = cell_spec(view.cell_types(i), i);
            int v[8];
            for (int j = 0; j < spec.n_corners; ++j)
                v[j] = vertex_of_point[static_cast<std::size_t>(conn[static_cast<std::size_t>(view.offsets(i)) + j])];
            const int attribute = attribute_of_type.at(type_ids(i));
            switch (spec.geometry)
            {
            case mfem::Geometry::TRIANGLE: serial->AddTri(v, attribute); break;
            case mfem::Geometry::SQUARE: serial->AddQuad(v, attribute); break;
            case mfem::Geometry::TETRAHEDRON: serial->AddTet(v, attribute); break;
            case mfem::Geometry::CUBE: serial->AddHex(v, attribute); break;
            case mfem::Geometry::PRISM: serial->AddWedge(v, attribute); break;
            case mfem::Geometry::PYRAMID: serial->AddPyramid(v, attribute); break;
            default: fail("unsupported geometry.");
            }
        }
        serial->FinalizeTopology();
        // No reordering of the element vertices (refine = false), and the orientation was fixed above.
        serial->Finalize(false, false);

        // ---- Partition (METIS, as libMesh) and the parallel mesh --------------------------------------------------
        std::vector<int> partitioning(n_cells_, 0);
        if (n_ranks_ > 1)
        {
            int *part = serial->GeneratePartitioning(n_ranks_, 1);
            std::copy(part, part + n_cells_, partitioning.begin());
            delete[] part;
        }
        mesh_ = std::make_unique<mfem::ParMesh>(comm_, *serial, partitioning.data());

        // The local elements are the serial ones of this rank, in serial order, with their vertices in the same order.
        std::vector<std::int64_t> local_vertex_point(static_cast<std::size_t>(mesh_->GetNV()), -1);
        for (std::size_t i = 0; i < n_cells_; ++i)
        {
            if (partitioning[i] != rank_)
                continue;
            const int e = static_cast<int>(cells_.size());
            if (e >= mesh_->GetNE())
                fail("the parallel mesh has fewer elements than its partition.");
            cells_.push_back(static_cast<std::int64_t>(i));
            subdomains_.push_back(type_ids(i));
            cell_types_.push_back(view.cell_types(i));
            point_offsets_.push_back(static_cast<std::int64_t>(point_ids_.size()));
            point_ids_.insert(point_ids_.end(), conn.begin() + view.offsets(i), conn.begin() + view.offsets(i + 1));
            const int *local = mesh_->GetElement(e)->GetVertices();
            const int *global = serial->GetElement(static_cast<int>(i))->GetVertices();
            for (int j = 0; j < mesh_->GetElement(e)->GetNVertices(); ++j)
                local_vertex_point[static_cast<std::size_t>(local[j])] = point_of_vertex[static_cast<std::size_t>(global[j])];
        }
        point_offsets_.push_back(static_cast<std::int64_t>(point_ids_.size()));
        if (static_cast<int>(cells_.size()) != mesh_->GetNE())
            fail("the parallel mesh and its partition have different elements.");
        serial.reset();
        vertex_points_ = std::move(local_vertex_point);
        for (const auto p : vertex_points_)
            if (p < 0)
                fail("a vertex of the parallel mesh has no point.");

        // ---- The points of the edges, faces and cells ----------------------------------------------------------
        std::unordered_map<std::vector<std::int64_t>, std::int64_t, KeyHash> high_points;
        element_points_center_.assign(cells_.size(), -1);
        for (int e = 0; e < n_local_elements(); ++e)
        {
            int n = 0;
            const std::int64_t *pts = element_points(e, n);
            const CellSpec &spec = cell_spec(cell_types_[static_cast<std::size_t>(e)], static_cast<std::size_t>(cells_[static_cast<std::size_t>(e)]));
            for (const HighOrderNode &node : spec.high)
            {
                if (static_cast<int>(node.corners.size()) == spec.n_corners && spec.dimension == 3)
                {
                    element_points_center_[static_cast<std::size_t>(e)] = pts[node.node];
                    continue;
                }
                if (static_cast<int>(node.corners.size()) == spec.n_corners && spec.dimension == 2 && node.corners.size() > 2)
                {
                    element_points_center_[static_cast<std::size_t>(e)] = pts[node.node];
                    continue;
                }
                std::vector<std::int64_t> corners;
                for (const int c : node.corners)
                    corners.push_back(pts[c]);
                high_points.emplace(key_of(corners), pts[node.node]);
            }
        }
        const auto lookup = [&](const mfem::Array<int> &vertices) -> std::int64_t
        {
            std::vector<std::int64_t> corners;
            for (int j = 0; j < vertices.Size(); ++j)
                corners.push_back(vertex_points_[static_cast<std::size_t>(vertices[j])]);
            const auto it = high_points.find(key_of(corners));
            return it == high_points.end() ? -1 : it->second;
        };
        mfem::Array<int> vertices;
        edge_points_.assign(static_cast<std::size_t>(mesh_->GetNEdges()), -1);
        for (int edge = 0; edge < mesh_->GetNEdges(); ++edge)
        {
            mesh_->GetEdgeVertices(edge, vertices);
            edge_points_[static_cast<std::size_t>(edge)] = lookup(vertices);
        }
        if (dimension_ == 3)
        {
            face_points_.assign(static_cast<std::size_t>(mesh_->GetNFaces()), -1);
            for (int face = 0; face < mesh_->GetNFaces(); ++face)
            {
                mesh_->GetFaceVertices(face, vertices);
                face_points_[static_cast<std::size_t>(face)] = lookup(vertices);
            }
        }

        // The points of the cell types of the local elements, built now: later lookups (on several threads) only read.
        for (int e = 0; e < n_local_elements(); ++e)
            point_rule(e);

        // ---- Curved geometry: the points of the quadratic cells as the mesh nodes -------------------------------
        if (geometric_order_ == 2)
        {
            auto *fec = new QuadraticH1Collection(dimension_, serendipity_);
            auto *fes = new mfem::ParFiniteElementSpace(mesh_.get(), fec, dimension_, mfem::Ordering::byVDIM);
            mesh_->SetNodalFESpace(fes);
            mfem::GridFunction &nodes = *mesh_->GetNodes();
            nodes.MakeOwner(fec); // The nodes own the space and the collection.
            const std::vector<std::int64_t> points = dof_points(*fes);
            for (int dof = 0; dof < fes->GetNDofs(); ++dof)
            {
                const std::int64_t p = points[static_cast<std::size_t>(dof)];
                if (p < 0)
                    fail("a node of the quadratic geometry has no point.");
                for (int c = 0; c < dimension_; ++c)
                    nodes(fes->DofToVDof(dof, c)) = view.points(static_cast<std::size_t>(p), static_cast<std::size_t>(c));
            }
        }
    }

    FEMesh::~FEMesh() = default;

    const std::int64_t *FEMesh::element_points(int e, int &n) const
    {
        const auto begin = point_offsets_[static_cast<std::size_t>(e)];
        n = static_cast<int>(point_offsets_[static_cast<std::size_t>(e) + 1] - begin);
        return point_ids_.data() + begin;
    }

    bool FEMesh::edge_corners(std::int64_t p, std::int64_t &a, std::int64_t &b) const
    {
        const auto it = edge_corners_.find(p);
        if (it == edge_corners_.end())
            return false;
        a = it->second.first;
        b = it->second.second;
        return true;
    }

    const std::vector<mfem::IntegrationPoint> &FEMesh::reference_points(int e) const
    {
        const int vtk = cell_types_[static_cast<std::size_t>(e)];
        auto it = reference_points_.find(vtk);
        if (it != reference_points_.end())
            return it->second;
        // Corners in the vertex order of the MFEM geometry (the VTK one); the other points are the centres of their
        // corners on the reference element (edge midpoints, face and cell centres).
        const CellSpec &spec = cell_spec(vtk, 0);
        const mfem::IntegrationRule &corners = *mfem::Geometries.GetVertices(spec.geometry);
        std::vector<mfem::IntegrationPoint> points(static_cast<std::size_t>(spec.n_nodes));
        for (int j = 0; j < spec.n_corners; ++j)
            points[static_cast<std::size_t>(j)] = corners.IntPoint(j);
        for (const HighOrderNode &node : spec.high)
        {
            double x[3] = {0.0, 0.0, 0.0};
            for (const int c : node.corners)
            {
                x[0] += corners.IntPoint(c).x / static_cast<double>(node.corners.size());
                x[1] += corners.IntPoint(c).y / static_cast<double>(node.corners.size());
                x[2] += corners.IntPoint(c).z / static_cast<double>(node.corners.size());
            }
            points[static_cast<std::size_t>(node.node)].Set3(x[0], x[1], x[2]);
        }
        return reference_points_.emplace(vtk, std::move(points)).first->second;
    }

    const mfem::IntegrationRule &FEMesh::point_rule(int e) const
    {
        const int vtk = cell_types_[static_cast<std::size_t>(e)];
        auto it = point_rules_.find(vtk);
        if (it != point_rules_.end())
            return *it->second;
        const std::vector<mfem::IntegrationPoint> &points = reference_points(e);
        auto rule = std::make_unique<mfem::IntegrationRule>(static_cast<int>(points.size()));
        for (std::size_t k = 0; k < points.size(); ++k)
            rule->IntPoint(static_cast<int>(k)) = points[k];
        return *point_rules_.emplace(vtk, std::move(rule)).first->second;
    }

    std::int64_t FEMesh::face_point(int face) const
    {
        if (dimension_ == 2)
            return edge_point(face);
        return face_points_[static_cast<std::size_t>(face)];
    }

    std::vector<std::int64_t> FEMesh::dof_points(const mfem::ParFiniteElementSpace &fes) const
    {
        std::vector<std::int64_t> points(static_cast<std::size_t>(fes.GetNDofs()), -1);
        mfem::Array<int> dofs;
        const auto assign = [&](const mfem::Array<int> &entity_dofs, std::int64_t point)
        {
            // One DOF per entity has the point of the entity; several (high order) have none.
            if (entity_dofs.Size() == 1)
                points[static_cast<std::size_t>(entity_dofs[0] >= 0 ? entity_dofs[0] : -1 - entity_dofs[0])] = point;
        };
        for (int v = 0; v < mesh_->GetNV(); ++v)
        {
            fes.GetVertexDofs(v, dofs);
            assign(dofs, vertex_point(v));
        }
        for (int edge = 0; edge < mesh_->GetNEdges(); ++edge)
        {
            fes.GetEdgeInteriorDofs(edge, dofs);
            assign(dofs, edge_point(edge));
        }
        if (dimension_ == 3)
            for (int face = 0; face < mesh_->GetNFaces(); ++face)
            {
                fes.GetFaceInteriorDofs(face, dofs);
                assign(dofs, face_point(face));
            }
        for (int e = 0; e < mesh_->GetNE(); ++e)
        {
            fes.GetElementInteriorDofs(e, dofs);
            assign(dofs, element_point(e));
        }
        return points;
    }
} // namespace pytec3d::fem
