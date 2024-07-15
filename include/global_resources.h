#ifndef GLOBAL_RESOURCES_H
#define GLOBAL_RESOURCES_H

#include <string>
// Forward declaration of LoadsBCObject if its full definition is not needed here
class LoadsBCObject;

class GlobalResources
{
public:
    static LoadsBCObject loadsBCObject;

    static void setLoadsBCObject(const LoadsBCObject &obj);
    static LoadsBCObject &getLoadsBCObject();
};

#endif // GLOBAL_RESOURCES_H