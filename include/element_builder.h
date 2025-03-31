#ifndef ELEMENT_BUILDER_H
#define ELEMENT_BUILDER_H

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
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vector>
#include <map>
#include "../include/time_stepper.h"
#include "../include/loads_bc_io.h"
#include "../include/materials_properties.h" // Add this line to include the definition of MaterialsProperties

using namespace libMesh;

class ElementDataManager;

class ElementBuilder
{
public:
    // Constructor
    ElementBuilder(std::shared_ptr<ElementDataManager> element_data_manager);

    void build_element(const ElementType &element_type, bool build_matrix);

    void build_solid97(const ElementType &element_type, bool build_matrix);
    void build_solid236(const ElementType &element_type, bool build_matrix);
    void build_plane13(const ElementType &element_type, bool build_matrix);

private:
    DenseMatrix<Number> Ke;
    DenseMatrix<Number> Ce;
    DenseVector<Number> Fe;

    DenseMatrix<Number> AAe;
    DenseMatrix<Number> Me;
    DenseMatrix<Number> JAe;
    DenseMatrix<Number> JVe;

    std::shared_ptr<ElementDataManager> element_data_manager;
};

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
          time_step(nullptr),
          time_substep(nullptr)
    {
    }

    inline void fill_mu_matrix(dof_id_type elem_id)
    {
        // get ansys_material_type from the materials map
        const MaterialsProperties &mat_properties = materials->at(elem_id);

        mu_matrix[0] = 4 * libMesh::pi * 1e-7 * mat_properties.murx;
        mu_matrix[1] = 4 * libMesh::pi * 1e-7 * mat_properties.mury;
        mu_matrix[2] = 4 * libMesh::pi * 1e-7 * mat_properties.murz;
    }

    inline void fill_rs_matrix(dof_id_type elem_id)
    {
        const MaterialsProperties &mat_properties = materials->at(elem_id);

        rs_matrix[0] = mat_properties.rsvx;
        rs_matrix[1] = mat_properties.rsvy;
        rs_matrix[2] = mat_properties.rsvz;
    }

    inline void compute_material_properties(dof_id_type elem_id, bool conductive)
    {
        // Call fill_mu_matrix and fill_rs_matrix to compute the material properties
        this->fill_mu_matrix(elem_id);
        if (conductive)
            this->fill_rs_matrix(elem_id);
    }

    inline Real compute_js(dof_id_type elem_id, unsigned int comp)
    {
        // Get the value from the "Js" array using the element ID
        if (time_step->get_ramped_load())
            return this->ramp_current(elem_id, comp);
        else
            return (*load_array)[elem_id][comp];
    }

    inline Real ramp_current(dof_id_type elem_id, unsigned int comp)
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
        return y1 + ((y2 - y1) / (time_step->get_time() - time_step->get_time_beginning())) * (time_substep->get_time() - time_step->get_time_beginning());
    }

    void ElementDataManager::update_time_step(Substep &time_substep)
    {
        this->time_substep = &time_substep;
        time_step = time_substep.get_parent();
    }

    void set_materials_map(std::unique_ptr<MaterialsMap> materials)
    {
        this->materials = std::move(materials);
    }

    void set_load_quantity(std::string load_name)
    {
        // Get the load quantity from the cell data
        load_array = &(*activeLoadsVector)[load_name];
        load_array = &(*oldLoadsVector)[load_name];
    }

    void updateLoad(std::unique_ptr<LoadArraysMap> newLoadsVector)
    {
        // Swap the old load with the new load vector
        std::swap(oldLoadsVector, activeLoadsVector);
        // Move the new load to the new load vector
        activeLoadsVector = std::move(newLoadsVector);
        // After this, oldLoad contains the previous contents of newLoadVector,
        // and newLoadVector contains the new load.
    }

    std::unique_ptr<LoadArraysMap> oldLoadsVector;
    std::unique_ptr<LoadArraysMap> activeLoadsVector;

    std::unique_ptr<MaterialsMap> materials;

    Real mu_matrix[3];
    Real rs_matrix[3];
    Real js[3];

    ElementLoad *load_array;
    ElementLoad *old_load_array;
    TimeStep *time_step;
    Substep *time_substep;
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
#endif // ELEMENT_BUILDER_H