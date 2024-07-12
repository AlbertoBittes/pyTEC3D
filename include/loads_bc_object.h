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

template <typename Sys>
class LoadsBCObject
{
private:
    Sys system;
    std::map<int, CurrentLoadsMap> timeStepCurrentLoadsMap;
    std::map<int, NodeBCMap> timeStepNodeBCMap;

public:
    // Constructor
    LoadsBCObject(const Sys &system) : system(system) {}

    // Method to set CurrentLoadsMap for a specific time step
    void setCurrentLoadsMap(int timeStep, const CurrentLoadsMap &currentLoadsMap)
    {
        timeStepCurrentLoadsMap[timeStep] = currentLoadsMap;
    }

    // Method to get CurrentLoadsMap for a specific time step
    CurrentLoadsMap getCurrentLoadsMap(int timeStep) const
    {
        auto it = timeStepCurrentLoadsMap.find(timeStep);
        if (it != timeStepCurrentLoadsMap.end())
        {
            return it->second;
        }
        return CurrentLoadsMap(); // Return empty map if not found
    }

    // Method to set NodeBCMap for a specific time step
    void setNodeBCMap(int timeStep, const NodeBCMap &nodeBCMap)
    {
        timeStepNodeBCMap[timeStep] = nodeBCMap;
    }

    // Method to get NodeBCMap for a specific time step
    NodeBCMap getNodeBCMap(int timeStep) const
    {
        auto it = timeStepNodeBCMap.find(timeStep);
        if (it != timeStepNodeBCMap.end())
        {
            return it->second;
        }
        return NodeBCMap(); // Return empty map if not found
    }

    NodeDataMatrix readInputDlist(const std::string &filename);

    CurrentLoadsMap buildCurrentLoadsMapBFElist(const std::string &filename);

    NodeBCMap buildNodeBCMapDlist(const std::string &dlist_filepath);
};

#endif // LOADS_BC_OBJECT_H