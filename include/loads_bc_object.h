// MatrixTypes.h
#ifndef LOADS_BC_OBJECT_H
#define LOADS_BC_OBJECT_H

#include <vector>
#include <string>
#include <cmath>
#include <map>
#include "libmesh/libmesh.h"

using namespace libMesh;
// Define CurrentLoadsMatrix type
typedef struct NodeData
{
    double AX = std::nan("");
    double AY = std::nan("");
    double AZ = std::nan("");
    double VOLT = std::nan("");
} NodeData;

typedef std::map<dof_id_type, std::vector<double>> CurrentLoadsMap;
typedef std::vector<std::tuple<dof_id_type, std::string, double>> NodeDataMatrix;
typedef std::map<dof_id_type, NodeData> NodeBCMap;

#endif // LOADS_BC_OBJECT_H