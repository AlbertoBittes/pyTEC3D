#ifndef NODE_BC_OBJECT_H
#define NODE_BC_OBJECT_H

#include "fe_data_object.h" // Include the base class
#include "globals.h"        // Include the base class
#include <libmesh/system.h>
#include <string>
#include <fstream>
#include <memory>
#include <unordered_map>

using namespace libMesh;

// Type aliases for boundary conditions
typedef std::unordered_map<Globals::GlobalVariableName, Real> NodeBoundaryConditions;
typedef std::unordered_map<dof_id_type, NodeBoundaryConditions> BoundaryConditionsMap;

class NodeBCObject : public FEDataObjectMap<BoundaryConditionsMap>
{
public:
    // Default constructor
    using FEDataObjectMap<BoundaryConditionsMap>::FEDataObjectMap;

    // Destructor
    virtual ~NodeBCObject() {}

    // // Setters
    // void set_system(System *sys) { system = sys; }
    // void set_update_boundary_conditions(bool update) { update_boundary_conditions = update; }

    // // Getters
    // System *get_system() const { return system; }
    // bool get_update_boundary_conditions() const { return update_boundary_conditions; }

    void read_data(std::ifstream &file) override;

    // Methods
    // void constrain() override; // Override from System::Constraint
    // void check_boundary_conditions(std::shared_ptr<LoadStep> prev_load_step);
    // void set_boundary_conditions(std::shared_ptr<LoadStep> prev_load_step);

private:
    // bool update_boundary_conditions; // Flag to update boundary conditions
    // System *system;                  // Pointer to the libMesh System
};

#endif // NODE_BC_OBJECT_H