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
    LoadsObject()
        : activeLoadsVector(nullptr),
          parent_load_step(nullptr)
    {
    }

    LoadsObject(const MeshBase &mesh, std::set<std::string> loads_labels)
        : activeLoadsVector(nullptr),
          parent_load_step(nullptr)
    {
        build_load_data_vectors(mesh, loads_labels);
    }

    void build_load_data_vectors(const MeshBase &mesh, std::set<std::string> loads_labels);

    void zero_loads_array();

    inline Real compute_load(dof_id_type elem_id, std::string load_label, unsigned int comp, LoadStep *prev_load_step = nullptr, Substep *time_substep = nullptr)
    {
        // Get the value from the "Js" array using the element ID
        if (parent_load_step->get_time_step()->get_ramped_load())
            return this->ramp_load(elem_id, load_label, comp, prev_load_step, time_substep);
        else
            return (*activeLoadsVector)[load_label][elem_id][comp];
    }

    inline Real ramp_load(dof_id_type elem_id, std::string load_label, unsigned int comp, LoadStep *prev_load_step, Substep *time_substep)
    {
        Real y1;
        Real y2;

        if (prev_load_step != nullptr)
            y1 = (*prev_load_step->get_loads()->get_active_loads_vector())[load_label][elem_id][comp];
        else
            y1 = 0;

        y2 = (*activeLoadsVector)[load_label][elem_id][comp];
        return y1 + ((y2 - y1) / (time_substep->get_parent().get_time() - time_substep->get_parent().get_time_beginning())) * (time_substep->get_time() - time_substep->get_parent().get_time_beginning());
    }

    void set_load_quantity(std::string load_name);

    void read_loads_array(std::istream &file_stream);

    // Getters
    LoadStep *get_parent_load_step() const { return parent_load_step; }
    const std::shared_ptr<LoadArraysMap> &get_active_loads_vector() const { return activeLoadsVector; }

    // Setters
    void set_parent_load_step(LoadStep *parent) { parent_load_step = parent; }
    void set_active_loads_vector(std::shared_ptr<LoadArraysMap> loads_vector) { activeLoadsVector = std::move(loads_vector); }

private:
    // Reference to the parent LoadStep
    LoadStep *parent_load_step;
    std::shared_ptr<LoadArraysMap> activeLoadsVector;
};

#endif // LOADS_OBJECT_H