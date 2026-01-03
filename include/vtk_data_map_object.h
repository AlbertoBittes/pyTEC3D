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
    struct SystemDataRepo
    {
        unsigned int sys_number;
        NumericVector<Number> *vec;
        NumericVector<Number> *old_vec;
        ApecMSystemBase::SystemType *system_type;
        Real initial_val;
    };

    VTKDataMapFEDataObject(std::shared_ptr<EquationSystems> external_map,
                           bool is_element_data_flag = false)
        : FEDataObjectMap<EquationSystems>(external_map), is_element_data(is_element_data_flag)
    {
    }

    std::vector<SystemDataRepo> systems;
    bool is_element_data = false;

    // User can call these to set up the systems before reading
    void add_dof_system()
    {
        ApecMLinearSystem *system = &data_map->get_system<ApecMLinearSystem>("ApecMLinearSystem");
        systems.push_back({system->number(),
                           system->request_vector("node_dof_constraints"),
                           system->request_vector("old_node_dof_constraints"),
                           &system->system_type,
                           std::numeric_limits<Real>::quiet_NaN()});
    }

    void add_nodal_loads_system()
    {
        ApecMExplicitSystem *nodal_loads_system = &data_map->get_system<ApecMExplicitSystem>("NodalLoadsSystem");
        systems.push_back({nodal_loads_system->number(),
                           nodal_loads_system->request_vector("new_loads"),
                           nodal_loads_system->request_vector("old_loads"),
                           &nodal_loads_system->system_type,
                           0});
    }

    void add_element_loads_system()
    {
        ApecMExplicitSystem *element_system = &data_map->get_system<ApecMExplicitSystem>("ElementLoadsSystem");
        systems.push_back({element_system->number(),
                           element_system->request_vector("new_loads"),
                           element_system->request_vector("old_loads"),
                           &element_system->system_type,
                           0});
    }

    // User sets this before reading
    void set_is_element_data(bool flag) { is_element_data = flag; }

    // The read_data method uses the current configuration
    void read_data(std::ifstream &file);
};

#endif // VTK_DATA_MAP_OBJECT_H
