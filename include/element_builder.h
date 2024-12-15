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
    ElementBuilder(DenseMatrix<Number> &Ke, DenseMatrix<Number> &Ce, DenseVector<Number> &Fe, MaterialsMap &materials, const std::unique_ptr<FEBase> &fe, vtkCellData *cell_data, TimeStep &ts, const DofMap &dof_map, QGauss &qrule);

    // Method
    void build_solid97_cond(const Elem *elem, DofMap &dof_map);
    void build_solid97(const Elem *elem, DofMap &dof_map);
    // Methods to fill matrices
    void fill_mu_matrix(dof_id_type elem_id);
    void fill_rs_matrix(dof_id_type elem_id);
    void compute_material_properties(dof_id_type elem_id, const unsigned int subdomain_id);
    void initialize(const Elem *elem, DofMap &dof_map, unsigned int num_vars, bool volt);

    const std::vector<Real> &get_JxW() const { return JxW; }
    const std::vector<std::vector<Real>> &get_phi() const { return phi; }
    const std::vector<std::vector<RealGradient>> &get_dphi() const { return dphi; }

    const std::vector<std::vector<double>> &get_mu_matrix() const { return mu_matrix; }
    const std::vector<std::vector<double>> &get_rs_matrix() const { return rs_matrix; }

    vtkDataArray *get_js_array() const { return js_array; }
    vtkDataArray *get_mat_array() const { return mat_array; }

    DenseMatrix<Number> &Ke;
    DenseMatrix<Number> &Ce;
    DenseVector<Number> &Fe;

    // Allocate memory for 3D nodal elements matrices with transient effects
    // I can have up to 4 variables per node/edge
    // References to momentum equation right hand sides
    DenseSubVector<Number> F[4];

    // References to velocity-velocity coupling blocks
    DenseSubMatrix<Number> KAA[3][3];

    // References to velocity-pressure coupling blocks
    DenseSubMatrix<Number> KAV[3];
    DenseSubMatrix<Number> KVA[3];
    DenseSubMatrix<Number> KVV;

    // References to velocity-velocity coupling blocks
    DenseSubMatrix<Number> CAA[3][3];

    // References to velocity-pressure coupling blocks
    DenseSubMatrix<Number> CAV[3];
    DenseSubMatrix<Number> CVA[3];
    DenseSubMatrix<Number> CVV;

    MaterialsMap &materials;
    std::vector<std::vector<double>> mu_matrix;
    std::vector<std::vector<double>> rs_matrix;

    // Attributes for JxW, phi, and dphi
    const std::vector<Real> &JxW;
    const std::vector<std::vector<Real>> &phi;
    const std::vector<std::vector<RealGradient>> &dphi;

    // Attributes for time step, time substep, and cell data
    TimeStep &ts;
    vtkDataArray *js_array;
    vtkDataArray *mat_array;
    std::vector<Real> old_solution_values;

    // Attributes for dof map and quadrature rule
    QGauss &qrule;
    unsigned int num_vars;
    unsigned int n_var_dofs;
};

#endif // ELEMENT_BUILDER_H