#ifndef NODE_DOF_CONSTRAINTS_H
#define NODE_DOF_CONSTRAINTS_H

#include <libmesh/equation_systems.h>
#include <string>
#include "../include/time_stepper.h"
#include <fstream>
#include <memory>
#include "../include/assembler.h" // Include the header for AssemblerLinear
// Forward declaration from libMesh
using namespace libMesh;

typedef std::unordered_map<std::string, Real> NodeBoundaryConditions;
typedef std::unordered_map<dof_id_type, NodeBoundaryConditions> BoundaryConditionsMap;

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
          time_step(nullptr)
    {
    }

    virtual ~NodalBCConstraint() {}

    void constrain();

    void check_boundary_conditions();

    void reinit_constraints_map()
    {
        active_bc_map->clear();
    }

    void update_time_step(TimeStep &time_step)
    {
        this->time_step = &time_step;
    }

    void set_boundary_conditions();

    void stash_bc_map()
    {
        // Swap the old load with the new load vector
        std::swap(old_bc_map, active_bc_map);
    }

    void read_boundary_conditions(std::string &file_name);

    void read_update_boundary_conditions(std::string file_name, TimeStep &time_step);

    void set_update_boundary_conditions(bool update)
    {
        update_boundary_conditions = update;
    }

    bool get_update_boundary_conditions() const
    {
        return update_boundary_conditions;
    }

    std::unique_ptr<BoundaryConditionsMap> old_bc_map;
    std::unique_ptr<BoundaryConditionsMap> active_bc_map;

private:
    TimeStep *time_step;
    bool update_boundary_conditions;
    System &system;
};

#endif // NODE_DOF_CONSTRAINTS_H