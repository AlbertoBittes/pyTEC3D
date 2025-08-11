#include <libmesh/mesh_base.h>
#include <libmesh/elem.h>
#include <libmesh/node.h>
#include <libmesh/remote_elem.h>
#include <libmesh/libmesh_logging.h>
#include <libmesh/mesh_tools.h>
#include <libmesh/face_inf_quad4.h>
#include <libmesh/face_inf_quad6.h>
#include <libmesh/cell_inf_prism6.h>
#include <libmesh/cell_inf_prism12.h>
#include <libmesh/cell_inf_hex8.h>
#include <libmesh/cell_inf_hex16.h>
#include <libmesh/cell_inf_hex18.h>
#include <libmesh/boundary_info.h>

#include "inf_elem_builder_updated.h"

void InfElemBuilderUpdated::build_inf_elem(const Point &origin,
                                           const bool x_sym = false,
                                           const bool y_sym = false,
                                           const bool z_sym = false,
                                           const bool be_verbose = false,
                                           std::set<std::pair<dof_id_type, unsigned int>> *inner_faces = nullptr)
{
    //     if (be_verbose)
    //     {
    // #ifdef DEBUG
    //         libMesh::out << " Building Infinite Elements:" << std::endl;
    //         libMesh::out << "  updating element neighbor tables..." << std::endl;
    // #else
    //         libMesh::out << " Verbose mode disabled in non-debug mode." << std::endl;
    // #endif
    //     }

    // update element neighbors
    this->mesh_ptr->find_neighbors();

    LOG_SCOPE("build_inf_elem()", "InfElemBuilder");

    // A set for storing element number, side number pairs.
    // pair.first == element number, pair.second == side number
    std::set<std::pair<dof_id_type, unsigned int>> faces;
    std::set<std::pair<dof_id_type, unsigned int>> ofaces;

    // A set for storing node numbers on the outer faces.
    std::set<dof_id_type> onodes;

    // The distance to the farthest point in the mesh from the origin
    Real max_r = 0.;

    // The index of the farthest point in the mesh from the origin
    int max_r_node = -1;

    // #ifdef DEBUG
    //     if (be_verbose)
    //     {
    //         libMesh::out << "  collecting boundary sides";
    //         if (x_sym || y_sym || z_sym)
    //             libMesh::out << ", skipping sides in symmetry planes..." << std::endl;
    //         else
    //             libMesh::out << "..." << std::endl;
    //     }
    // #endif

    auto &boundary_info = mesh_ptr->get_boundary_info();
    std::vector<boundary_id_type> side_boundary_ids;
    // Iterate through all elements and sides, collect indices of all active
    // boundary sides in the faces set. Skip sides which lie in symmetry planes.
    // Later, sides of the inner boundary will be sorted out.
    for (const auto &elem : mesh_ptr->active_element_ptr_range())
        for (auto s : elem->side_index_range())
            if (elem->neighbor_ptr(s) == nullptr)
            {
                boundary_info.boundary_ids(elem, s, side_boundary_ids);
                if (std::find(side_boundary_ids.begin(), side_boundary_ids.end(), 110) != side_boundary_ids.end())
                {
                    // note that it is safe to use the Elem::side() method,
                    // which gives a non-full-ordered element
                    std::unique_ptr<Elem> side(elem->build_side_ptr(s));

                    // bool flags for symmetry detection
                    bool sym_side = false;
                    bool on_x_sym = true;
                    bool on_y_sym = true;
                    bool on_z_sym = true;

                    // Loop over the nodes to check whether they are on the symmetry planes,
                    // and therefore sufficient to use a non-full-ordered side element
                    for (const Node &node : side->node_ref_range())
                    {
                        const dof_id_type node_id = node.id();
                        const Point dist_from_origin =
                            this->mesh_ptr->point(node_id) - origin;

                        if (x_sym)
                            if (std::abs(dist_from_origin(0)) > 1.e-3)
                                on_x_sym = false;

                        if (y_sym)
                            if (std::abs(dist_from_origin(1)) > 1.e-3)
                                on_y_sym = false;

                        if (z_sym)
                            if (std::abs(dist_from_origin(2)) > 1.e-3)
                                on_z_sym = false;

                        //       if (x_sym)
                        // if (std::abs(dist_from_origin(0)) > 1.e-6)
                        //   on_x_sym=false;

                        //       if (y_sym)
                        // if (std::abs(dist_from_origin(1)) > 1.e-6)
                        //   on_y_sym=false;

                        //       if (z_sym)
                        // if (std::abs(dist_from_origin(2)) > 1.e-6)
                        //   on_z_sym=false;

                        // find the node most distant from origin

                        Real r = dist_from_origin.norm();
                        if (r > max_r)
                        {
                            max_r = r;
                            max_r_node = node_id;
                        }
                    }

                    sym_side = (x_sym && on_x_sym) || (y_sym && on_y_sym) || (z_sym && on_z_sym);

                    if (!sym_side)
                        faces.emplace(elem->id(), s);
                }
            } // neighbor(s) == nullptr

    //  If a boundary side has one node on the outer boundary,
    //  all points of this side are on the outer boundary.
    //  Start with the node most distant from origin, which has
    //  to be on the outer boundary, then recursively find all
    //  sides and nodes connected to it. Found sides are moved
    //  from faces to ofaces, nodes are collected in onodes.
    //  Here, the search is done iteratively, because, depending on
    //  the mesh, a very high level of recursion might be necessary.
    if (max_r_node >= 0)
        // include the possibility of the 1st element being most far away.
        // Only the case of no outer boundary is to be excluded.
        onodes.insert(max_r_node);

    {
        auto face_it = faces.begin();
        auto face_end = faces.end();
        unsigned int facesfound = 0;
        while (face_it != face_end)
        {
            std::pair<dof_id_type, unsigned int> p = *face_it;

            // This has to be a full-ordered side element,
            // since we need the correct n_nodes,
            std::unique_ptr<Elem> side(this->mesh_ptr->elem_ref(p.first).build_side_ptr(p.second));

            bool found = false;
            for (const Node &node : side->node_ref_range())
                if (onodes.count(node.id()))
                {
                    found = true;
                    break;
                }

            // If a new oface is found, include its nodes in onodes
            if (found)
            {
                for (const Node &node : side->node_ref_range())
                    onodes.insert(node.id());

                ofaces.insert(p);
                face_it = faces.erase(face_it); // increment is done here

                facesfound++;
            }

            else
                ++face_it; // increment is done here

            // If at least one new oface was found in this cycle,
            // do another search cycle.
            if (facesfound > 0 && face_it == faces.end())
            {
                facesfound = 0;
                face_it = faces.begin();
            }
        }
    }

#ifdef DEBUG
    if (be_verbose)
        libMesh::out << "  found "
                     << faces.size()
                     << " inner and "
                     << ofaces.size()
                     << " outer boundary faces"
                     << std::endl;
#endif

    // When the user provided a non-null pointer to
    // inner_faces, that implies he wants to have
    // this std::set.  For now, simply copy the data.
    if (inner_faces != nullptr)
        *inner_faces = faces;

    // free memory, clear our local variable, no need
    // for it any more.
    faces.clear();

    // outer_nodes maps onodes to their duplicates
    std::map<dof_id_type, Node *> outer_nodes;

    // We may need to pick our own object ids in parallel
    dof_id_type old_max_node_id = mesh_ptr->max_node_id();
    dof_id_type old_max_elem_id = mesh_ptr->max_elem_id();

    // Likewise with our unique_ids
#ifdef LIBMESH_ENABLE_UNIQUE_ID
    unique_id_type old_max_unique_id = mesh_ptr->parallel_max_unique_id();
#endif

    // for each boundary node, add an outer_node with
    // double distance from origin.
    for (const auto &dof : onodes)
    {
        Point p = (Point(this->mesh_ptr->point(dof)) * 2) - origin;
        if (mesh_ptr->is_serial())
        {
            // Add with a default id in serial
            outer_nodes[dof] = this->mesh_ptr->add_point(p);
        }
        else
        {
            // Pick a unique id in parallel
            Node &bnode = mesh_ptr->node_ref(dof);
            dof_id_type new_id = bnode.id() + old_max_node_id;
            std::unique_ptr<Node> new_node = Node::build(p, new_id);
            new_node->processor_id() = bnode.processor_id();
#ifdef LIBMESH_ENABLE_UNIQUE_ID
            new_node->set_unique_id(old_max_unique_id + bnode.id());
#endif

            outer_nodes[dof] =
                this->mesh_ptr->add_node(std::move(new_node));
        }
    }

#ifdef DEBUG
    // for verbose, remember n_elem
    dof_id_type n_conventional_elem = this->mesh_ptr->n_elem();
#endif

    // build Elems based on boundary side type
    for (auto &p : ofaces)
    {
        Elem &belem = this->mesh_ptr->elem_ref(p.first);

        // build a full-ordered side element to get the base nodes
        std::unique_ptr<Elem> side(belem.build_side_ptr(p.second));

        // create cell depending on side type, assign nodes,
        // use braces to force scope.
        bool is_higher_order_elem = false;

        std::unique_ptr<Elem> el;
        switch (side->type())
        {
            // 3D infinite elements
            // TRIs
        case TRI3:
            el = Elem::build(INFPRISM6);
            break;

        case TRI6:
            el = Elem::build(INFPRISM12);
            is_higher_order_elem = true;
            break;

            // QUADs
        case QUAD4:
            el = Elem::build(INFHEX8);
            break;

        case QUAD8:
            el = Elem::build(INFHEX16);
            is_higher_order_elem = true;
            break;

        case QUAD9:
            el = Elem::build(INFHEX18);

            // the method of assigning nodes (which follows below)
            // omits in the case of QUAD9 the bubble node; therefore
            // we assign these first by hand here.
            el->set_node(16, side->node_ptr(8));
            el->set_node(17, outer_nodes[side->node_id(8)]);
            is_higher_order_elem = true;
            break;

            // 2D infinite elements
        case EDGE2:
            el = Elem::build(INFQUAD4);
            break;

        case EDGE3:
            el = Elem::build(INFQUAD6);
            el->set_node(4, side->node_ptr(2));
            break;

            // 1D infinite elements not supported
        default:
            libMesh::out << "InfElemBuilder::build_inf_elem(Point, bool, bool, bool, bool): "
                         << "invalid face element "
                         << std::endl;
            continue;
        }

        const unsigned int n_base_vertices = side->n_vertices();

        // On a distributed mesh, manually assign unique ids to the new
        // element, and make sure any RemoteElem neighbor links are set.
        if (!mesh_ptr->is_serial())
        {
            el->processor_id() = belem.processor_id();

            // We'd better not have elements with more than 6 sides
            const unsigned int max_sides = 6;
            libmesh_assert_less_equal(el->n_sides(), max_sides);
            el->set_id(belem.id() * max_sides + p.second + old_max_elem_id);

#ifdef LIBMESH_ENABLE_UNIQUE_ID
            el->set_unique_id(old_max_unique_id + old_max_node_id +
                              belem.id() * max_sides + p.second);
#endif

            // If we have a remote neighbor on a boundary element side
            if (belem.dim() > 1)
                for (auto s : belem.side_index_range())
                    if (belem.neighbor_ptr(s) == remote_elem)
                    {
                        // Find any corresponding infinite element side
                        std::unique_ptr<const Elem> remote_side(belem.build_side_ptr(s));

                        for (auto inf_s : el->side_index_range())
                        {
                            // The base side 0 shares all vertices but isn't
                            // remote
                            if (!inf_s)
                                continue;

                            // But another side, one which shares enough
                            // vertices to show it's the same side, is.
                            unsigned int n_shared_vertices = 0;
                            for (unsigned int i = 0; i != n_base_vertices; ++i)
                                for (auto &node : remote_side->node_ref_range())
                                    if (side->node_ptr(i) == &node &&
                                        el->is_node_on_side(i, inf_s))
                                        ++n_shared_vertices;

                            if (n_shared_vertices + 1 >= belem.dim())
                            {
                                el->set_neighbor(inf_s, const_cast<RemoteElem *>(remote_elem));
                                break;
                            }
                        }
                    }
        }

        // assign vertices to the new infinite element
        for (unsigned int i = 0; i < n_base_vertices; i++)
        {
            el->set_node(i, side->node_ptr(i));
            el->set_node(i + n_base_vertices, outer_nodes[side->node_id(i)]);
        }

        // when this is a higher order element,
        // assign also the nodes in between
        if (is_higher_order_elem)
        {
            // n_safe_base_nodes is the number of nodes in \p side
            // that may be safely assigned using below for loop.
            // Actually, n_safe_base_nodes is _identical_ with el->n_vertices(),
            // since for QUAD9, the 9th node was already assigned above
            const unsigned int n_safe_base_nodes = el->n_vertices();

            for (unsigned int i = n_base_vertices; i < n_safe_base_nodes; i++)
            {
                el->set_node(i + n_base_vertices, side->node_ptr(i));
                el->set_node(i + n_safe_base_nodes,
                             outer_nodes[side->node_id(i)]);
            }
        }

        // add infinite element to mesh
        this->mesh_ptr->add_elem(std::move(el));
    } // for

#ifdef DEBUG
    mesh_ptr->libmesh_assert_valid_parallel_ids();

    if (be_verbose)
        libMesh::out << "  added "
                     << this->mesh_ptr->n_elem() - n_conventional_elem
                     << " infinite elements and "
                     << onodes.size()
                     << " nodes to the mesh"
                     << std::endl
                     << std::endl;
#endif
}