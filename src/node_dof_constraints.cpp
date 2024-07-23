#include "libmesh/equation_systems.h"
#include "libmesh/node.h"
#include "libmesh/dof_map.h"
#include "libmesh/mesh.h"
#include "libmesh/linear_implicit_system.h" // Include the header file for LinearImplicitSystem
#include "../include/node_dof_constraints.h"
#include "../include/global_resources.h"
#include "../include/loads_bc_object.h"

using namespace libMesh;

void node_dof_constraints(libMesh::EquationSystems &equation_systems, const std::string &system_name)
{
    // Example implementation
    // This function should implement the logic to apply constraints to the nodes of the system
    // identified by system_name within the provided equation_systems object.
    // The actual implementation would depend on the specific constraints to be applied.

    // Retrieve the system by name
    LinearImplicitSystem &system = equation_systems.get_system<LinearImplicitSystem>(system_name);
    LoadsBCObject loads_bc_object = GlobalResources::getLoadsBCObject();
    const NodeBCMap &node_bc_map = loads_bc_object.getNodeBCMap(1);
    const MeshBase &mesh = equation_systems.get_mesh();
    for (const auto &node : mesh.local_node_ptr_range())
    {
        // dof_id_type dof_loop = node->dof_number(const unsigned int s,
        //                                         const unsigned int var,
        //                                         const unsigned int comp) const

        auto row_bc = node_bc_map.find(node->id() + 1);
        if (row_bc != node_bc_map.end())
        { // Check if the element exists
            const NodeData &values = row_bc->second;
            Real ax_bc = values.AX;
            if (!std::isnan(ax_bc))
            {
                dof_id_type dof_id_ax = node->dof_number(0, 0, 0);
                DofConstraintRow constraint_row;
                if (ax_bc == 0.0)
                {
                    system.get_dof_map().add_constraint_row(dof_id_ax, constraint_row, false);
                }
                else
                {
                    system.get_dof_map().add_constraint_row(dof_id_ax, constraint_row, ax_bc, false);
                }
                // sffdh
            }
            Real ay_bc = values.AY;
            if (!std::isnan(ay_bc))
            {
                dof_id_type dof_id_ay = node->dof_number(0, 0, 1);
                DofConstraintRow constraint_row;
                if (ay_bc == 0.0)
                {
                    system.get_dof_map().add_constraint_row(dof_id_ay, constraint_row, false);
                }
                else
                {
                    system.get_dof_map().add_constraint_row(dof_id_ay, constraint_row, ay_bc, false);
                }
            }
            Real az_bc = values.AZ;
            if (!std::isnan(az_bc))
            {
                dof_id_type dof_id_az = node->dof_number(0, 0, 2);
                DofConstraintRow constraint_row;
                if (az_bc == 0.0)
                {
                    system.get_dof_map().add_constraint_row(dof_id_az, constraint_row, false);
                }
                else
                {
                    system.get_dof_map().add_constraint_row(dof_id_az, constraint_row, az_bc, false);
                }
            }
        }
    }
}