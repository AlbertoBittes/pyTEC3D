#include "../include/loads_bc_mat_object.h"
#include "../include/global_resources.h" // Include the full definition of LoadsBCObject

// Definition of static member
LoadsBCMatObject GlobalResources::loadsBCObject;

void GlobalResources::setLoadsBCObject(const LoadsBCMatObject &obj)
{
    loadsBCObject = obj;
}

LoadsBCMatObject &GlobalResources::getLoadsBCObject()
{
    return loadsBCObject;
}