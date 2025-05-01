#ifndef NODE_BC_OBJECT_H
#define NODE_BC_OBJECT_H

#include <libmesh/system.h>
#include <string>
#include "../include/load_step.h"
#include <fstream>
#include <memory>
// Forward declaration from libMesh
using namespace libMesh;

typedef std::unordered_map<std::string, std::pair<unsigned int, Real>> NodeBoundaryConditions;
typedef std::unordered_map<dof_id_type, NodeBoundaryConditions> BoundaryConditionsMap;

class NodeBCObject : public System::Constraint
{
public:
    // Default constructor
    NodeBCObject()
        : active_bc_map(nullptr), system(nullptr), parent_load_step(nullptr), update_boundary_conditions(true)
    {
    }

    // Constructor with System pointer
    NodeBCObject(System *sys)
        : active_bc_map(nullptr), system(sys), parent_load_step(nullptr), update_boundary_conditions(true)
    {
    }

    // Destructor
    virtual ~NodeBCObject() {}

    // Setters
    void set_system(System *sys) { system = sys; }
    void set_active_bc_map(std::shared_ptr<BoundaryConditionsMap> &bc_map) { active_bc_map = std::move(bc_map); }
    void set_update_boundary_conditions(bool update) { update_boundary_conditions = update; }
    void set_parent_load_step(std::shared_ptr<LoadStep> parent) { parent_load_step = parent; }

    // Getters
    System *get_system() const { return system; }
    std::shared_ptr<BoundaryConditionsMap> get_active_bc_map() const { return active_bc_map; }
    bool get_update_boundary_conditions() const { return update_boundary_conditions; }
    std::shared_ptr<LoadStep> get_parent_load_step() const { return parent_load_step; }

    // Methods
    void constrain();
    void check_boundary_conditions(std::shared_ptr<LoadStep> prev_load_step);
    void set_boundary_conditions(std::shared_ptr<LoadStep> prev_load_step);
    void read_boundary_conditions(std::istream &file);

private:
    bool update_boundary_conditions;
    System *system;
    std::shared_ptr<BoundaryConditionsMap> active_bc_map;
    std::shared_ptr<LoadStep> parent_load_step; // Pointer to the parent LoadStep
};

#endif // NODE_BC_OBJECT_H