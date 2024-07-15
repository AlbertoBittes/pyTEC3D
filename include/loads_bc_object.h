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
    Real AX = std::nan("");
    Real AY = std::nan("");
    Real AZ = std::nan("");
    Real VOLT = std::nan("");
} NodeData;

typedef std::map<dof_id_type, std::vector<Real>> CurrentLoadsMap;
typedef std::vector<std::tuple<dof_id_type, std::string, Real>> NodeDataMatrix;
typedef std::map<dof_id_type, NodeData> NodeBCMap;

class LoadsBCObject
{
private:
    std::map<int, CurrentLoadsMap> timeStepCurrentLoadsMap;
    std::map<int, NodeBCMap> timeStepNodeBCMap;
    std::string loads_folderpath;
    std::string bc_folderpath;
    std::string inp_name;

public:
    LoadsBCObject();
    // Constructor
    LoadsBCObject(std::string loads_folderpath, std::string bc_folderpath, std::string inp_name)
        : loads_folderpath(loads_folderpath),
          bc_folderpath(bc_folderpath),
          inp_name(inp_name)
    {
        setCurrentLoadsMap(1, buildCurrentLoadsMapBFElist(loads_folderpath + "/" + inp_name + "_loads_" + std::to_string(1) + ".txt"));
        setNodeBCMap(1, buildNodeBCMapDlist(bc_folderpath + "/" + inp_name + "_bcs_" + std::to_string(1) + ".txt"));
    }

    // Method to set CurrentLoadsMap for a specific time step
    void setCurrentLoadsMap(int timeStep, const CurrentLoadsMap &currentLoadsMap)
    {
        timeStepCurrentLoadsMap[timeStep] = currentLoadsMap;
    }

    // Method to get CurrentLoadsMap for a specific time step
    // Method to get CurrentLoadsMap for a specific time step
    const CurrentLoadsMap &getCurrentLoadsMap(int timeStep) const
    {
        static const CurrentLoadsMap emptyMap; // Static empty map to return by reference if not found
        auto it = timeStepCurrentLoadsMap.find(timeStep);
        if (it != timeStepCurrentLoadsMap.end())
        {
            return it->second;
        }
        return emptyMap; // Return reference to the static empty map if not found
    }

    // Method to set NodeBCMap for a specific time step
    void setNodeBCMap(int timeStep, const NodeBCMap &nodeBCMap)
    {
        timeStepNodeBCMap[timeStep] = nodeBCMap;
    }

    const NodeBCMap &getNodeBCMap(int timeStep) const
    {
        static const NodeBCMap emptyMap; // Static empty map to return by reference if not found
        auto it = timeStepNodeBCMap.find(timeStep);
        if (it != timeStepNodeBCMap.end())
        {
            return it->second;
        }
        return emptyMap; // Return reference to the static empty map if not found
    }

    NodeDataMatrix readInputDlist(const std::string &filename);

    CurrentLoadsMap buildCurrentLoadsMapBFElist(const std::string &filename);

    NodeBCMap buildNodeBCMapDlist(const std::string &dlist_filepath);
};

#endif // LOADS_BC_OBJECT_H