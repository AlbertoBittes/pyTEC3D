#ifndef NODE_BC_OBJECT_H
#define NODE_BC_OBJECT_H

#include <libmesh/equation_systems.h>
#include "fe_data_object.h" // Include the base class
#include "globals.h"        // Include the globals header for GlobalVariableName
#include "apec_m_eq_system_object.h"
// using namespace libMesh;

// The NodeBCObject holds the nodal DOF imposed by the user (D in APDL), they will be used by the boundary
// conditions object to set the constraints in the system matrix
class NodeLoadsDOFConstraints : public FEDataObjectMap<ApecMEquationSystems>
{
public:
    // Default constructor
    NodeLoadsDOFConstraints(std::ifstream &file, std::shared_ptr<ApecMEquationSystems> external_map)
        : FEDataObjectMap<ApecMEquationSystems>(external_map)
    {
        this->read_data(file); // Call read_data to initialize the data map
    }

    // Destructor
    virtual ~NodeLoadsDOFConstraints() {}

    void read_data(std::ifstream &file) override;
};

#endif // NODE_BC_OBJECT_H