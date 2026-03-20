#ifndef VTK_DATA_MAP_OBJECT_H
#define VTK_DATA_MAP_OBJECT_H

#include <libmesh/elem.h>
#include <libmesh/system.h>
#include <libmesh/explicit_system.h>
#include <map>
#include <fstream> // Include this header for std::ifstream
#include "fe_data_object.h"
#include "apec_m_system.h"

using namespace libMesh;

class VTKDataMapFEDataObject : public FEDataObjectMap<EquationSystems>
{
public:
    VTKDataMapFEDataObject(std::shared_ptr<EquationSystems> external_map,
                           bool is_element_data_flag = false)
        : FEDataObjectMap<EquationSystems>(external_map), is_element_data(is_element_data_flag)
    {
    }

    bool is_element_data = false;

    // User sets this before reading
    void set_is_element_data(bool flag) { is_element_data = flag; }

    // The read_data method uses the current configuration
    void read_data(std::ifstream &file);
};

#endif // VTK_DATA_MAP_OBJECT_H
