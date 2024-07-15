#include "../include/loads_bc_object.h"
#include "../include/global_resources.h" // Include the full definition of LoadsBCObject

// Definition of static member
LoadsBCObject GlobalResources::loadsBCObject;

void GlobalResources::setLoadsBCObject(const LoadsBCObject &obj)
{
    loadsBCObject = obj;
}

LoadsBCObject &GlobalResources::getLoadsBCObject()
{
    return loadsBCObject;
}