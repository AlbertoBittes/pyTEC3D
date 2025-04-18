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

    inline void fill_mu_matrix(dof_id_type elem_id);

    inline void fill_rs_matrix(dof_id_type elem_id);

    inline void compute_material_properties(dof_id_type elem_id, bool conductive);

    inline Real compute_js(dof_id_type elem_id, unsigned int comp);

    inline Real ramp_current(dof_id_type elem_id, unsigned int comp);

    inline void update_time_step(const Substep &time_substep);

    inline void set_load_quantity(std::string load_name);

    inline void stashLoadsVector();

    void read_loads_arrays(std::string &file_name);

    void read_update_loads_vectors(std::string file_name);

    std::unique_ptr<LoadArraysMap> oldLoadsVector;
    std::unique_ptr<LoadArraysMap> activeLoadsVector;
    std::unique_ptr<MaterialsMap> materials;

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