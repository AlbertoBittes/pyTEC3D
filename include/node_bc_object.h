#ifndef NODE_BC_OBJECT_H
#define NODE_BC_OBJECT_H

#include "fe_data_object.h" // Include the base class
#include "globals.h"        // Include the globals header for GlobalVariableName
// using namespace libMesh;

// Type aliases for boundary conditions
typedef std::unordered_map<Globals::GlobalVariableName, Real> NodeBoundaryConditions;
typedef std::unordered_map<dof_id_type, NodeBoundaryConditions> BoundaryConditionsMap;

// The NodeBCObject holds the nodal DOF imposed by the user (D in APDL), they will be used by the boundary
// conditions object to set the constraints in the system matrix
class NodeBCObject : public FEDataObjectMap<BoundaryConditionsMap>
{
public:
    // Default constructor
    NodeBCObject(std::ifstream &file)
        : FEDataObjectMap<BoundaryConditionsMap>()
    {
        this->read_data(file); // Call read_data to initialize the data map
    }

    // Destructor
    virtual ~NodeBCObject() {}

    void read_data(std::ifstream &file) override;
};

#endif // NODE_BC_OBJECT_H