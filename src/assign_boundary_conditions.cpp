// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

// <h1>Introduction Example 1 - Creation of a Mesh Object</h1>
// \author Benjamin S. Kirk
// \date 2003
//
// This is the fourth example program.  It builds on
// the third example program by showing how to formulate
// the code in a dimension-independent way.  Very minor
// changes to the example will allow the problem to be
// solved in one, two or three dimensions.
//
// This example will also introduce the PerfLog class
// as a way to monitor your code's performance.  We will
// use it to instrument the matrix assembly code and look
// for bottlenecks where we should focus optimization efforts.
//
// This example also shows how to extend example 3 to run in
// parallel.  Notice how little has changed!

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/linear_implicit_system.h"
// Define the Finite Element object.
#include "libmesh/boundary_info.h"
// Define Gauss quadrature rules.

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"

// Define the PerfLog, a performance logging utility.
// It is useful for timing events in a code and giving
// you an idea where bottlenecks lie.

// The definition of a geometric element
#include "libmesh/elem.h"

// To impose Dirichlet boundary conditions
#include "libmesh/dirichlet_boundaries.h"

#include "libmesh/string_to_enum.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/zero_function.h"

using namespace libMesh;
// Define a new function that receives as input a pointer to mesh object and returns void
void assign_boundary_conditions(LinearImplicitSystem &system, unsigned int a_var)
{
    const MeshBase &mesh = system.get_mesh();

    const unsigned int dim = mesh.mesh_dimension();

    ZeroFunction<Number> zero;

    std::set<boundary_id_type> boundary_ids;
    if (dim == 2)
    {
        boundary_ids.insert(0);

        DirichletBoundary dirichlet_bc(boundary_ids, {a_var}, zero);
        system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);
    }
    // if (dim == 3)
    // {
    //     boundary_ids.insert(0);
    //     DirichletBoundary dirichlet_bc(boundary_ids(0), {a_var}, zero);

    //     boundary_ids.insert(1);
    //     boundary_ids.insert(2);

    //     system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);
    // }
}

void assign_boundary_ids(Mesh &mesh)
{
    const unsigned int dim = mesh.mesh_dimension();
    const boundary_id_type bound_id = 0;
    // Loop over all the nodes in the mesh
    for (const auto &node : mesh.node_ptr_range())
    {

        // If the node is on the boundary, set the boundary ID
        if (dim == 2)
        {
            // Put the variable pointed by node into a variable

            const dof_id_type node_id = node->id();
            Node nod = *node;
            if (std::abs(nod(0)) == 1 || std::abs(nod(1)) == 1)
            {
                mesh.get_boundary_info().add_node(node_id, bound_id);
            }
        }
        // if (dim == 3)
        // {
        //     // Put the variable pointed by node into a variable

        //     const dof_id_type node_id = node->id();
        //     Node nod = *node;
        //     if (std::abs(nod(0)) == 0 || std::abs(nod(0)) == 0.25)
        //     {
        //         const boundary_id_type bound_id = 0;
        //         mesh.get_boundary_info().add_node(node_id, bound_id);
        //     }
        //     else if (std::abs(nod(1)) == 0 || std::abs(nod(1)) == 0.25)
        //     {
        //         const boundary_id_type bound_id = 1;
        //         mesh.get_boundary_info().add_node(node_id, bound_id);
        //     }
        //     else if (std::abs(nod(2)) == -0.125 || std::abs(nod(1)) == 0.125)
        //     {
        //         const boundary_id_type bound_id = 2;
        //         mesh.get_boundary_info().add_node(node_id, bound_id);
        //     }
    }
}
