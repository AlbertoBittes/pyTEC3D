#ifndef ELEMENT_DATA_MANAGER_H
#define ELEMENT_DATA_MANAGER_H

#include <libmesh/elem.h>
#include <libmesh/system.h>
#include <libmesh/dense_matrix.h>
#include <libmesh/dense_vector.h>
#include <libmesh/dense_submatrix.h>
#include <libmesh/dense_subvector.h>
#include <libmesh/dof_map.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/perf_log.h>
#include <libmesh/fe.h>
#include <vector>
#include <map>
#include <fstream> // Include this header for std::ifstream
#include "../include/time_stepper.h"
#include "../include/element_types.h"        // Include the header file for ElementTypesMap
#include "../include/materials_properties.h" // Add this line to include the definition of MaterialsProperties

typedef std::unordered_map<dof_id_type, std::vector<Real>> ElementLoad;
typedef std::unordered_map<std::string, ElementLoad> LoadArraysMap;

class ElementDataManager
{
public:
    // Constructor
    ElementDataManager()
        : oldLoadsVector(nullptr),
          activeLoadsVector(nullptr),
          materials(nullptr),
          load_array(nullptr),
          old_load_array(nullptr),
          time_substep(nullptr)
    {
    }

    ElementDataManager(const MeshBase &mesh, std::set<std::string> loads_labels)
        : oldLoadsVector(nullptr),
          activeLoadsVector(nullptr),
          materials(nullptr),
          load_array(nullptr),
          old_load_array(nullptr),
          time_substep(nullptr)
    {
        build_load_data_vectors(mesh, loads_labels);
    }

    void build_load_data_vectors(const MeshBase &mesh, std::set<std::string> loads_labels);

    void zero_load_arrays();

    inline void ElementDataManager::fill_mu_matrix(dof_id_type elem_id)
    {
        // get ansys_material_type from the materials map
        const MaterialsProperties &mat_properties = materials->at(elem_id);

        mu_matrix[0] = 4 * libMesh::pi * 1e-7 * mat_properties.murx;
        mu_matrix[1] = 4 * libMesh::pi * 1e-7 * mat_properties.mury;
        mu_matrix[2] = 4 * libMesh::pi * 1e-7 * mat_properties.murz;
    }

    inline void ElementDataManager::fill_rs_matrix(dof_id_type elem_id)
    {
        const MaterialsProperties &mat_properties = materials->at(elem_id);

        rs_matrix[0] = mat_properties.rsvx;
        rs_matrix[1] = mat_properties.rsvy;
        rs_matrix[2] = mat_properties.rsvz;
    }

    inline void ElementDataManager::compute_material_properties(dof_id_type elem_id, bool conductive)
    {
        // Call fill_mu_matrix and fill_rs_matrix to compute the material properties
        this->fill_mu_matrix(elem_id);
        if (conductive)
            this->fill_rs_matrix(elem_id);
    }

    inline Real ElementDataManager::compute_js(dof_id_type elem_id, unsigned int comp)
    {
        // Get the value from the "Js" array using the element ID
        if (time_substep->get_parent().get_ramped_load())
            return this->ramp_current(elem_id, comp);
        else
            return (*load_array)[elem_id][comp];
    }

    inline Real ElementDataManager::ramp_current(dof_id_type elem_id, unsigned int comp)
    {
        Real y1;
        Real y2;
        if (oldLoadsVector != nullptr)
        {
            y1 = (*old_load_array)[elem_id][comp];
            y2 = (*load_array)[elem_id][comp];
        }
        else
        {
            y1 = 0;
            y2 = (*load_array)[elem_id][comp];
        }
        return y1 + ((y2 - y1) / (time_substep->get_parent().get_time() - time_substep->get_parent().get_time_beginning())) * (time_substep->get_time() - time_substep->get_parent().get_time_beginning());
    }

    void update_time_step(const Substep &time_substep)
    {
        this->time_substep = &time_substep;
    };

    void set_load_quantity(std::string load_name);

    void stashLoadsVector();

    void read_loads_arrays(std::string &file_name);

    void read_update_loads_vectors(std::string file_name);

    std::unique_ptr<LoadArraysMap> oldLoadsVector;
    std::unique_ptr<LoadArraysMap> activeLoadsVector;
    std::shared_ptr<MaterialsMap> materials;

    Real mu_matrix[3];
    Real rs_matrix[3];

    ElementLoad *load_array;
    ElementLoad *old_load_array;
    const Substep *time_substep;
    // std::shared_ptr<TimeStepper> timeStepper;
    // Real time_step_starting_time;
    // Real time_step_ending_time;
    // Real current_time;
    // bool ramped;
    // std::map<dof_id_type, DenseMatrix<Number>> ke_matrix_map;
    // std::map<dof_id_type, DenseMatrix<Number>> ce_matrix_map;
    // std::map<dof_id_type, DenseVector<Number>> rhs_map;

    // std::map<dof_id_type, DenseMatrix<Number>> aae_matrix_derived_map;
    // std::map<dof_id_type, DenseMatrix<Number>> me_matrix_derived_map;
    // std::map<dof_id_type, DenseMatrix<Number>> jae_matrix_derived_map;
    // std::map<dof_id_type, DenseMatrix<Number>> jve_matrix_derived_map;
    // Here I assume static dofs
    // std::map<dof_id_type, std::vector<dof_id_type>> dof_indices;
};

#endif // ELEMENT_DATA_MANAGER_H