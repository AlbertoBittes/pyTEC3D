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

class ElementBuilder
{
public:
    // Constructor
    ElementBuilder(DenseMatrix<Number> &Ke, DenseMatrix<Number> &Ce, DenseVector<Number> &Fe, MaterialsMap &materials, const std::unique_ptr<FEVectorBase> &fe_vec, const std::unique_ptr<FEBase> &fe_scalar, vtkCellData *cell_data, const TimeStep &ts);
    void build_element_matrix(const Elem *elem, const DofMap &dof_map, subdomain_id_type subdomain_id, const std::string &formulation, const unsigned int dim, const std::string &fe_family);
    void build_element_rhs(const Elem *elem, const DofMap &dof_map, subdomain_id_type subdomain_id, const std::string &formulation, const unsigned int dim, const std::string &fe_family);
    void build_solid97_cond(const Elem *elem, const DofMap &dof_map);
    void build_solid97(const Elem *elem, const DofMap &dof_map);
    void build_solid97_rhs(const Elem *elem, const DofMap &dof_map);
    void build_solid236(const Elem *elem, const DofMap &dof_map);
    void build_solid236_cond(const Elem *elem, const DofMap &dof_map);
    void build_plane13(const Elem *elem, const DofMap &dof_map);
    void build_plane13_rhs(const Elem *elem, const DofMap &dof_map);
    void build_plane13_cond(const Elem *elem, const DofMap &dof_map);
    void fill_mu_matrix(dof_id_type elem_id);
    void fill_rs_matrix(dof_id_type elem_id);
    void compute_material_properties(dof_id_type elem_id, subdomain_id_type subdomain_id);

    void build_element_matrix_derived(DenseMatrix<Number> &AAe, DenseMatrix<Number> &Me, const Elem *elem, const DofMap &dof_map, subdomain_id_type subdomain_id, const std::string &formulation, const unsigned int dim, const std::string &fe_family, unsigned int n_A_dofs);
    // std::vector<std::vector<dof_id_type>> *get_dof_indices_var() { return &dof_indices_var; }
    // unsigned int get_n_var_dofs() const { return n_var_dofs; }
    // unsigned int get_n_var_dofs_old() const { return n_var_dofs_old; }
    // void set_n_var_dofs_old(unsigned int n) { n_var_dofs_old = n; }
    // void set_n_var_dofs(unsigned int n) { n_var_dofs = n; }

private:
    DenseMatrix<Number> &Ke;
    DenseMatrix<Number> &Ce;
    DenseVector<Number> &Fe;
    MaterialsMap &materials;
    std::vector<std::vector<double>> mu_matrix;
    std::vector<std::vector<double>> rs_matrix;
    const std::vector<Real> &JxW_A;
    const std::vector<Real> &JxW_V;
    const std::vector<std::vector<RealGradient>> &phi_A;
    const std::vector<std::vector<RealTensor>> &dphi_A;
    const std::vector<std::vector<RealGradient>> &curl_phi_A;
    // const std::vector<std::vector<Real>> &div_phi_A;
    const std::vector<std::vector<Real>> &phi_V;
    const std::vector<std::vector<RealGradient>> &dphi_V;
    const TimeStep &ts;
    vtkDataArray *mat_array;
    unsigned int n_var_dofs;
    // std::vector<std::vector<dof_id_type>> dof_indices_var;
};

#endif // ELEMENT_BUILDER_H