#ifndef MATERIALS_PROPERTIES_H
#define MATERIALS_PROPERTIES_H

#include <unordered_map>
#include <string>
#include <libmesh/libmesh.h>

using namespace libMesh;

// Structure to hold material properties
struct MaterialsProperties
{
    Real murx, mury, murz;
    Real rsvx, rsvy, rsvz;
    Real perx, pery, perz;
};

typedef std::unordered_map<int, MaterialsProperties> MaterialsMap;

// Function to read material properties from a CSV file
MaterialsMap readMaterialsProperties(const std::string &filename);

#endif // MATERIAL_PROPERTIES_H