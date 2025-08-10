#ifndef LOADS_OBJECT_H
#define LOADS_OBJECT_H

#include <libmesh/elem.h>
#include <libmesh/system.h>
#include <libmesh/explicit_system.h>
#include <map>
#include <fstream> // Include this header for std::ifstream
#include "fe_data_object.h"
#include "apec_m_system_objects.h"

using namespace libMesh;

// Type aliases for loads
typedef std::unordered_map<dof_id_type, std::vector<Real>> ElementLoad;
typedef std::unordered_map<std::string, ElementLoad> LoadArraysMap;

// The LoadsObject holds the body loads (Js) in a useful structure, to be easily interfaced with the ElementBuilder
// (BFE loads in APDL)
class LoadsObject : public FEDataObjectMap<ApecMLoadsSystem>
{
public:
    // Constructor that calls the exact FEDataObjectMap constructor you specified
    LoadsObject(std::ifstream &file, std::shared_ptr<ApecMLoadsSystem> external_map)
        : FEDataObjectMap<ApecMLoadsSystem>(external_map)
    {
        this->read_data(file); // Call read_data to initialize the data map
        // No additional initialization
    }

    // Read loads array from a binary file stream
    void read_data(std::ifstream &file) override;
};

#endif // LOADS_OBJECT_H