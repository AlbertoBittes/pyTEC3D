#ifndef LOADS_OBJECT_H
#define LOADS_OBJECT_H

#include <libmesh/elem.h>
#include <libmesh/system.h>
#include <map>
#include <fstream>         // Include this header for std::ifstream
#include "element_types.h" // Include the header file for ElementTypesMap
#include "fe_data_object.h"
#include "globals.h"

class LoadStep;

// Type aliases for loads
typedef std::unordered_map<dof_id_type, std::vector<Real>> ElementLoad;
typedef std::unordered_map<std::string, ElementLoad> LoadArraysMap;

class LoadsObject : public FEDataObjectMap<ExplicitSystem>
{
public:
    // Constructor that calls the exact FEDataObjectMap constructor you specified
    LoadsObject(std::ifstream &file, std::shared_ptr<ExplicitSystem> external_map)
        : FEDataObjectMap<ExplicitSystem>(file, external_map)
    {
        // No additional initialization
    }

    // // Build load data vectors
    // void build_load_data_vectors(const MeshBase &mesh);

    // // Zero out the loads array
    // void zero_loads_array();

    // Read loads array from a binary file stream
    void read_data(std::ifstream &file) override;

    // Swap the active and old loads vectors
    // void swap_loads_vector(std::unique_ptr<LoadArraysMap> loads_vector)
    // {
    //     *oldLoadsVector = *data_map;
    //     data_map = std::move(loads_vector);
    // }

    // // Getters for old loads vector
    // LoadArraysMap *get_old_loads_vector() const { return oldLoadsVector.get(); }

    // Setter for old loads vector

    // private:
    //     std::unique_ptr<LoadArraysMap> oldLoadsVector; // Old loads vector
};

#endif // LOADS_OBJECT_H