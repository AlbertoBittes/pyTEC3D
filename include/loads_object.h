#ifndef LOADS_OBJECT_H
#define LOADS_OBJECT_H

#include <libmesh/elem.h>
#include <libmesh/system.h>
#include <map>
#include <fstream>                    // Include this header for std::ifstream
#include "../include/element_types.h" // Include the header file for ElementTypesMap
#include "../include/fe_data_object.h"
#include "../include/globals.h"

class LoadStep;

// Type aliases for loads
typedef std::unordered_map<dof_id_type, std::vector<Real>> ElementLoad;
typedef std::unordered_map<std::string, ElementLoad> LoadArraysMap;

class LoadsObject : public FEDataObject<LoadArraysMap, LoadStep>
{
public:
    // Default constructor
    LoadsObject(LoadStep *parent_object = nullptr,
                std::shared_ptr<LoadArraysMap> map = nullptr)
        : FEDataObject(parent_object, map),
          oldLoadsVector(std::make_shared<LoadArraysMap>()) {}

    // Constructor with MeshBase and load labels
    LoadsObject(const MeshBase &mesh,
                LoadStep *parent_object = nullptr,
                std::shared_ptr<LoadArraysMap> map = nullptr)
        : FEDataObject(parent_object, map),
          oldLoadsVector(std::make_shared<LoadArraysMap>())
    {
        build_load_data_vectors(mesh);
    }

    // Build load data vectors
    void build_load_data_vectors(const MeshBase &mesh);

    // Zero out the loads array
    void zero_loads_array();

    // Set the load quantity
    void set_load_quantity(std::string &load_name);

    // Read loads array from a binary file stream
    void read_data(std::ifstream &file) override;

    // Swap the active and old loads vectors
    void swap_loads_vector(std::shared_ptr<LoadArraysMap> loads_vector)
    {
        *oldLoadsVector = *data_map;
        data_map = std::move(loads_vector);
    }

    // Getters for old loads vector
    std::shared_ptr<LoadArraysMap> get_old_loads_vector() const { return oldLoadsVector; }

    // Setter for old loads vector
    void set_old_loads_vector(std::shared_ptr<LoadArraysMap> loads_vector) { oldLoadsVector = std::move(loads_vector); }

private:
    std::shared_ptr<LoadArraysMap> oldLoadsVector; // Old loads vector
};

#endif // LOADS_OBJECT_H