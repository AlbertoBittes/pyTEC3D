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
#include "../include/materials_properties.h" // Add this line to include the definition of MaterialsProperties

using namespace libMesh;

class ElementDataManager;

class ElementBuilder
{
public:
    // Constructor
    ElementBuilder(DenseMatrix<Number> &Ke, DenseMatrix<Number> &Ce, DenseVector<Number> &Fe, std::shared_ptr<ElementDataManager> element_data_manager, const std::unique_ptr<FEVectorBase> &fe_vec, const std::unique_ptr<FEBase> &fe_scalar);
    void build_element_matrix(unsigned int n_dofs_var, bool conductive, const std::string &formulation, const unsigned int dim, const std::string &fe_family);
    void build_element_rhs(unsigned int n_dofs_var, bool conductive, const std::string &formulation, const unsigned int dim, const std::string &fe_family);
    void build_solid97_cond();
    void build_solid97();
    void build_solid97_rhs();
    void build_solid236();
    void build_solid236_cond();
    void build_plane13();
    void build_plane13_rhs();
    void build_plane13_cond();

    void build_element_matrix_derived(DenseMatrix<Number> &AAe, DenseMatrix<Number> &Me, DenseMatrix<Number> &JAe, DenseMatrix<Number> &JVe, bool conductive, const std::string &formulation, const unsigned int dim, const std::string &fe_family, unsigned int n_A_dofs);
    // std::vector<std::vector<dof_id_type>> *get_dof_indices_var() { return &dof_indices_var; }
    // unsigned int get_n_var_dofs() const { return n_var_dofs; }
    // unsigned int get_n_var_dofs_old() const { return n_var_dofs_old; }
    // void set_n_var_dofs_old(unsigned int n) { n_var_dofs_old = n; }
    // void set_n_var_dofs(unsigned int n) { n_var_dofs = n; }

private:
    DenseMatrix<Number> &Ke;
    DenseMatrix<Number> &Ce;
    DenseVector<Number> &Fe;
    const std::vector<Real> &JxW_A;
    const std::vector<Real> &JxW_V;
    const std::vector<std::vector<RealGradient>> &phi_A;
    const std::vector<std::vector<RealTensor>> &dphi_A;
    const std::vector<std::vector<RealGradient>> &curl_phi_A;
    // const std::vector<std::vector<Real>> &div_phi_A;
    const std::vector<std::vector<Real>> &phi_V;
    const std::vector<std::vector<RealGradient>> &dphi_V;
    std::shared_ptr<ElementDataManager> element_data_manager;
    unsigned int n_var_dofs;
    // std::vector<std::vector<dof_id_type>> dof_indices_var;
};

class ElementDataManager
{
public:
    // Constructor
    ElementDataManager(MaterialsMap &materials, vtkCellData *cell_data, std::shared_ptr<TimeStepper> timeStepper);

    inline void fill_mu_matrix()
    {
        // get ansys_material_type from the materials map
        const MaterialsProperties &mat_properties = materials.at(mat_array->GetTuple1(local_elem_id));

        mu_matrix[0] = 4 * libMesh::pi * 1e-7 * mat_properties.murx;
        mu_matrix[1] = 4 * libMesh::pi * 1e-7 * mat_properties.mury;
        mu_matrix[2] = 4 * libMesh::pi * 1e-7 * mat_properties.murz;
    }

    inline void fill_rs_matrix()
    {
        const MaterialsProperties &mat_properties = materials.at(mat_array->GetTuple1(local_elem_id));

        rs_matrix[0] = mat_properties.rsvx;
        rs_matrix[1] = mat_properties.rsvy;
        rs_matrix[2] = mat_properties.rsvz;
    }

    inline void compute_material_properties(bool conductive)
    {
        // Call fill_mu_matrix and fill_rs_matrix to compute the material properties
        this->fill_mu_matrix();
        if (conductive)
            this->fill_rs_matrix();
    }

    inline Real compute_js(unsigned int comp)
    {
        // Get the value from the "Js" array using the element ID
        if (ramped)
            return this->ramp_current(comp);
        else
            return js_array->GetComponent(local_elem_id, comp);
    }

    inline Real ramp_current(unsigned int comp)
    {
        Real y1;
        Real y2;
        if (js_array_old)
        {
            y1 = js_array_old->GetComponent(local_elem_id, comp);
            y2 = js_array->GetComponent(local_elem_id, comp);
        }
        else
        {
            y1 = 0;
            y2 = js_array->GetComponent(local_elem_id, comp);
        }
        return y1 + ((y2 - y1) / (time_step_ending_time - time_step_starting_time)) * (current_time - time_step_starting_time);
    }

    void update_time_step();

    dof_id_type local_elem_id;
    MaterialsMap &materials;
    Real mu_matrix[3];
    Real rs_matrix[3];
    Real js[3];
    vtkCellData *cell_data;
    vtkDataArray *mat_array;
    vtkDataArray *js_array;
    vtkDataArray *js_array_old;
    std::shared_ptr<TimeStepper> timeStepper;
    Real time_step_starting_time;
    Real time_step_ending_time;
    Real current_time;
    bool ramped;
    std::map<dof_id_type, DenseMatrix<Number>> ke_matrix_map;
    std::map<dof_id_type, DenseMatrix<Number>> ce_matrix_map;
    std::map<dof_id_type, DenseVector<Number>> rhs_map;

    std::map<dof_id_type, DenseMatrix<Number>> aae_matrix_derived_map;
    std::map<dof_id_type, DenseMatrix<Number>> me_matrix_derived_map;
    std::map<dof_id_type, DenseMatrix<Number>> jae_matrix_derived_map;
    std::map<dof_id_type, DenseMatrix<Number>> jve_matrix_derived_map;
    // Here I assume static dofs
    std::map<dof_id_type, std::vector<dof_id_type>> dof_indices;
};
#endif // ELEMENT_BUILDER_H