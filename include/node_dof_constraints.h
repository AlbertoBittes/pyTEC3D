#ifndef NODE_DOF_CONSTRAINTS_H
#define NODE_DOF_CONSTRAINTS_H

#include "vtkPointData.h"
#include <libmesh/equation_systems.h>
#include <string>
#include "../include/time_stepper.h"
#include "../include/loads_bc_io.h" // Add this line to include the definition of BoundaryConditionsMap

// Forward declaration from libMesh
using namespace libMesh;

class NodalBCConstraint : public System::Constraint
{
private:
public:
    bool update_boundary_conditions;
    // Modified constructor to take a vtkPointData pointer
    NodalBCConstraint(System &system)
        : system(system),
          update_boundary_conditions(true),
          old_bc_map(nullptr),
          active_bc_map(nullptr),
          bc_node(nullptr),
          bc_node_old(nullptr),
          time_step(nullptr)
    {
    }

    virtual ~NodalBCConstraint() {}

    void constrain();

    void check_boundary_conditions();

    void update_time_step(TimeStep &time_step)
    {
        this->time_step = &time_step;
    }

    // void set_bc_quantity(std::string load_name)
    // {
    //     // Get the load quantity from the cell data
    //     bc_node = &(*active_bc_map)[load_name];
    //     bc_node_old = &(*old_bc_map)[load_name];
    // }

    void update_bc_map(std::unique_ptr<BoundaryConditionsMap> new_bc_map)
    {
        // Swap the old load with the new load vector
        std::swap(old_bc_map, active_bc_map);
        // Move the new load to the new load vector
        active_bc_map = std::move(new_bc_map);
        // After this, oldLoad contains the previous contents of newLoadVector,
        // and newLoadVector contains the new load.
    }

    std::unique_ptr<BoundaryConditionsMap> old_bc_map;
    std::unique_ptr<BoundaryConditionsMap> active_bc_map;

    NodeBoundaryConditions *bc_node;
    NodeBoundaryConditions *bc_node_old;
    TimeStep *time_step;

    bool update_boundary_conditions;

    System &system;
};

#endif // NODE_DOF_CONSTRAINTS_H