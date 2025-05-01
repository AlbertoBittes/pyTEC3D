#ifndef LOADS_OBJECT_H
#define LOADS_OBJECT_H

#include <libmesh/elem.h>
#include <libmesh/system.h>
#include <map>
#include <fstream>                    // Include this header for std::ifstream
#include "../include/element_types.h" // Include the header file for ElementTypesMap
#include "../include/load_step.h"     // Add this line to include the definition of MaterialsProperties

class LoadStep;

typedef std::unordered_map<dof_id_type, std::vector<Real>> ElementLoad;
typedef std::unordered_map<std::string, ElementLoad> LoadArraysMap;

class LoadsObject
{
public:
    // Constructor
    LoadsObject(std::shared_ptr<LoadStep> parent_load_step = nullptr)
        : activeLoadsVector(std::make_shared<LoadArraysMap>()),
          oldLoadsVector(std::make_shared<LoadArraysMap>()),
          parent_load_step(parent_load_step)
    {
    }

    LoadsObject(const MeshBase &mesh, std::set<std::string> &loads_labels, std::shared_ptr<LoadStep> parent_load_step = nullptr)
        : activeLoadsVector(std::make_shared<LoadArraysMap>()),
          oldLoadsVector(std::make_shared<LoadArraysMap>()),
          parent_load_step(parent_load_step)
    {
        build_load_data_vectors(mesh, loads_labels);
    }

    void build_load_data_vectors(const MeshBase &mesh, std::set<std::string> &loads_labels);

    void zero_loads_array();

    // inline Real compute_load(dof_id_type elem_id, std::shared_ptr<LoadArraysMap> &loads_vector, std::string &load_label, unsigned int comp, std::shared_ptr<LoadStep> &prev_load_step, Substep *time_substep = nullptr)
    // {
    //     // Get the value from the "Js" array using the element ID
    //     if (parent_load_step->get_time_step()->get_ramped_load())
    //         return this->ramp_load(elem_id, loads_vector, load_label, comp, prev_load_step, time_substep);
    //     else
    //         return (*loads_vector)[load_label][elem_id][comp];
    // }

    // inline Real ramp_load(dof_id_type elem_id, std::shared_ptr<LoadArraysMap> &loads_vector, std::string &load_label, unsigned int comp, std::shared_ptr<LoadStep> &prev_load_step, Substep *time_substep)
    // {
    //     Real y1;
    //     Real y2;

    //     if (prev_load_step != nullptr)
    //         y1 = (*prev_load_step->get_loads()->get_active_loads_vector())[load_label][elem_id][comp];
    //     else
    //         y1 = 0;

    //     y2 = (*loads_vector)[load_label][elem_id][comp];
    //     return y1 + ((y2 - y1) / (time_substep->get_parent().get_time() - time_substep->get_parent().get_time_beginning())) * (time_substep->get_time() - time_substep->get_parent().get_time_beginning());
    // }

    void set_load_quantity(std::string &load_name);

    void read_loads_array(std::istream &file_stream);

    // Getters
    std::shared_ptr<LoadStep> get_parent_load_step() const { return parent_load_step; }
    std::shared_ptr<LoadArraysMap> &get_active_loads_vector() { return activeLoadsVector; }
    std::shared_ptr<LoadArraysMap> &get_old_loads_vector() { return oldLoadsVector; }

    // Setters
    void set_parent_load_step(std::shared_ptr<LoadStep> parent) { parent_load_step = parent; }
    void set_active_loads_vector(std::shared_ptr<LoadArraysMap> loads_vector) { activeLoadsVector = std::move(loads_vector); }
    void set_old_loads_vector(std::shared_ptr<LoadArraysMap> loads_vector) { oldLoadsVector = std::move(loads_vector); }

    void swap_loads_vector(std::shared_ptr<LoadArraysMap> loads_vector)
    {
        *oldLoadsVector = *activeLoadsVector;
        activeLoadsVector = std::move(loads_vector);
    }

private:
    // Reference to the parent LoadStep
    std::shared_ptr<LoadStep> parent_load_step;
    std::shared_ptr<LoadArraysMap> activeLoadsVector;
    std::shared_ptr<LoadArraysMap> oldLoadsVector;
};

#endif // LOADS_OBJECT_H