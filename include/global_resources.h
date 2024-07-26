#ifndef GLOBAL_RESOURCES_H
#define GLOBAL_RESOURCES_H

#include <string>
// Forward declaration of LoadsBCObject if its full definition is not needed here
class LoadsBCMatObject;

class GlobalResources
{
public:
    static LoadsBCMatObject loadsBCObject;

    static void setLoadsBCObject(const LoadsBCMatObject &obj);
    static LoadsBCMatObject &getLoadsBCObject();
};

#endif // GLOBAL_RESOURCES_H