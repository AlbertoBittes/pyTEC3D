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
    ElementBuilder(DenseMatrix<Number> &Ke, DenseMatrix<Number> &Ce, DenseVector<Number> &Fe, MaterialsMap &materials, const std::unique_ptr<FEBase> &fe, vtkCellData *cell_data, const TimeStep &ts, QGauss &qrule);

    void build_element_matrix(const Elem *elem, const DofMap &dof_map, const unsigned int subdomain_id, const std::string &formulation, const unsigned int dim, const std::string &fe_family);
    void build_element_rhs(const Elem *elem, const DofMap &dof_map, const unsigned int subdomain_id, const std::string &formulation, const unsigned int dim, const std::string &fe_family);
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
    void compute_material_properties(dof_id_type elem_id, const unsigned int subdomain_id);
    void initialize(const Elem *elem, const DofMap &dof_map, unsigned int num_vars, bool volt);
    void initialize_rhs(const Elem *elem, const DofMap &dof_map, unsigned int num_vars);

    const std::vector<Real> &get_JxW() const { return JxW; }
    const std::vector<std::vector<Real>> &get_phi() const { return phi; }
    const std::vector<std::vector<RealGradient>> &get_dphi() const { return dphi; }

    const std::vector<std::vector<double>> &get_mu_matrix() const { return mu_matrix; }
    const std::vector<std::vector<double>> &get_rs_matrix() const { return rs_matrix; }

    vtkDataArray *get_js_array() const { return js_array; }
    vtkDataArray *get_mat_array() const { return mat_array; }

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
    DenseSubVector<Number> F[4];
    DenseSubMatrix<Number> KAA[3][3];
    DenseSubMatrix<Number> KAV[3];
    DenseSubMatrix<Number> KVA[3];
    DenseSubMatrix<Number> CAA[3][3];
    DenseSubMatrix<Number> CAV[3];
    DenseSubMatrix<Number> CVA[3];
    DenseSubMatrix<Number> CVV;
    std::vector<std::vector<double>> mu_matrix;
    std::vector<std::vector<double>> rs_matrix;
    const std::vector<Real> &JxW;
    const std::vector<std::vector<Real>> &phi;
    const std::vector<std::vector<RealGradient>> &dphi;
    const TimeStep &ts;
    QGauss &qrule;
    vtkSmartPointer<vtkDataArray> js_array;
    vtkSmartPointer<vtkDataArray> mat_array;
    unsigned int n_var_dofs;
    unsigned int n_var_dofs_old;
    // std::vector<std::vector<dof_id_type>> dof_indices_var;
};

#endif // ELEMENT_BUILDER_H