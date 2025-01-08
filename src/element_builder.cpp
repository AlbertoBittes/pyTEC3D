#include "../include/element_builder.h"

// Constructor
ElementBuilder::ElementBuilder(DenseMatrix<Number> &Ke, DenseMatrix<Number> &Ce, DenseVector<Number> &Fe, MaterialsMap &materials, const std::unique_ptr<FEVectorBase> &fe_A, const std::unique_ptr<FEBase> &fe_V, vtkCellData *cell_data, const TimeStep &ts, QGauss &qrule)
    : Ke(Ke), Ce(Ce), Fe(Fe), materials(materials),
      F(DenseSubVector<Number>(Fe)),
      KAA(DenseSubMatrix<Number>(Ke)),
      CAA(DenseSubMatrix<Number>(Ce)),
      CAV(DenseSubMatrix<Number>(Ce)),
      CVA(DenseSubMatrix<Number>(Ce)),
      CVV(DenseSubMatrix<Number>(Ce)),
      mu_matrix(3, std::vector<double>(3, 0.0)),
      rs_matrix(3, std::vector<double>(3, 0.0)),
      JxW(fe_A->get_JxW()),
      phi_A(fe_A->get_phi()),
      dphi_A(fe_A->get_dphi()),
      curl_phi_A(fe_A->get_curl_phi()),
      div_phi_A(fe_A->get_div_phi()),
      phi_V(fe_V->get_phi()),
      dphi_V(fe_V->get_dphi()),
      ts(ts),
      qrule(qrule),
      mat_array(nullptr),
      n_var_dofs(0),
      n_var_dofs_old(0)
{
    std::string array_name = "ansys_material_type";
    mat_array = cell_data->GetArray(array_name.c_str());
    if (mat_array == nullptr)
    {
        libmesh_error_msg("Error: 'ansys_material_type' array not found in vtkCellData.");
    }
}

void ElementBuilder::initialize(const Elem *elem, const DofMap &dof_map, unsigned int num_vars, bool volt)
{
    // this->get_dof_indices_var()->resize(num_vars);

    // for (unsigned int var = 0; var < num_vars; var++)
    //     dof_map.dof_indices(elem, this->get_dof_indices_var()[var], var);
    // This will have to be fixed once you introduce non-nodal elements!

    // if (this->n_var_dofs == this->n_var_dofs_old)
    // {
    //     return;
    // }
    // else
    // {
    // Here you can get num_vars from the size of the variable (3 for A)
    // in the future, make it more general with whatever number of variables and their size
    // KAA.reposition(0, 0, n_var_dofs * num_vars, n_var_dofs * num_vars);

    // if (volt)
    // {
    //     CAA.reposition(0, 0, n_var_dofs * num_vars, n_var_dofs * num_vars);

    //     // KAV[var_i].reposition(var_i * n_var_dofs, num_vars * n_var_dofs, n_var_dofs, n_var_dofs);
    //     // KVA[var_i].reposition(num_vars * n_var_dofs, var_i * n_var_dofs, n_var_dofs, n_var_dofs);

    //     CAV.reposition(0, num_vars * n_var_dofs, n_var_dofs, n_var_dofs);
    //     CVA.reposition(num_vars * n_var_dofs, 0, n_var_dofs, n_var_dofs);

    //     CVV.reposition(num_vars * n_var_dofs, num_vars * n_var_dofs, n_var_dofs, n_var_dofs);
    // }
    // // KVV.reposition(num_vars * n_var_dofs, num_vars * n_var_dofs, n_var_dofs, n_var_dofs);
    // if (ts.get_transient_effect())
    // this->n_var_dofs_old = this->n_var_dofs;
    // }
}

void ElementBuilder::initialize_rhs(const Elem *elem, const DofMap &dof_map, unsigned int num_vars)
{
    // this->get_dof_indices_var()->resize(num_vars);

    // for (unsigned int var = 0; var < num_vars; var++)
    //     dof_map.dof_indices(elem, this->get_dof_indices_var()[var], var);
    // This will have to be fixed once you introduce non-nodal elements!

    this->n_var_dofs = elem->n_nodes();
    // for (unsigned int var = 0; var < num_vars; var++)
    //     F[var].reposition(var * n_var_dofs, n_var_dofs);
    F.reposition(0, n_var_dofs * num_vars);
}
// Method
void ElementBuilder::build_solid97_cond(const Elem *elem, const DofMap &dof_map)
{
    // Numeric ids corresponding to each variable in the system
    unsigned int num_vars = 3;
    this->initialize(elem, dof_map, num_vars, true);
    // Now we will build the element matrix.  This involves
    // a double loop to integrate the test functions (i) against
    // the trial functions (j).
    //
    // We have split the numeric integration into two loops
    // so that we can log the matrix and right-hand-side
    // computation separately.
    //
    // Now start logging the element matrix computation
    Real CAV_value;
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
        for (unsigned int dof_i = 0; dof_i < n_var_dofs * num_vars; dof_i++)
        {
            for (unsigned int dof_j = 0; dof_j < n_var_dofs * num_vars; dof_j++)
            {
                for (unsigned int k = 0; k < 3; k++)
                {

                    Ke(dof_i, dof_j) += JxW[qp] / mu_matrix[k][k] * (dphi_A[dof_j][qp].row(k) * dphi_A[dof_i][qp].row(k));
                    Ce(dof_i, dof_j) += JxW[qp] / rs_matrix[k][k] * (phi_A[dof_i][qp](k) * phi_A[dof_j][qp](k));
                    if (dof_j < n_var_dofs)
                    {
                        CAV_value = JxW[qp] / rs_matrix[k][k] * phi_A[dof_i][qp](k) * dphi_V[dof_j][qp](k);
                        Ce(dof_i, num_vars * n_var_dofs + dof_j) += CAV_value;
                        Ce(num_vars * n_var_dofs + dof_j, dof_i) += CAV_value;
                        if (dof_i < n_var_dofs)
                        {
                            Ce(num_vars * n_var_dofs + dof_i, num_vars * n_var_dofs + dof_j) += JxW[qp] * dphi_V[dof_i][qp](k) / rs_matrix[k][k] * dphi_V[dof_j][qp](k);
                        }
                    }
                }
            }
        }
    }
}
void ElementBuilder::build_solid97(const Elem *elem, const DofMap &dof_map)
{
    const unsigned int num_vars = 3;
    this->initialize(elem, dof_map, num_vars, false);

    // Now we will build the element matrix.  This involves
    // a double loop to integrate the test functions (i) against
    // the trial functions (j).
    //
    // We have split the numeric integration into two loops
    // so that we can log the matrix and right-hand-side
    // computation separately.
    //
    // Now start logging the element matrix computation
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {

        for (unsigned int dof_i = 0; dof_i < n_var_dofs * num_vars; dof_i++)
        {
            for (unsigned int dof_j = 0; dof_j < n_var_dofs * num_vars; dof_j++)
            {
                for (unsigned int k = 0; k < 3; k++)
                {
                    Ke(dof_i, dof_j) += JxW[qp] / mu_matrix[k][k] * (dphi_A[dof_j][qp].row(k) * dphi_A[dof_i][qp].row(k));
                }
            }
        }
    }
    // libmesh_error_msg("build_solid97 is not implemented yet.");
}
void ElementBuilder::build_solid97_rhs(const Elem *elem, const DofMap &dof_map)
{
    const unsigned int num_vars = 3;
    this->initialize_rhs(elem, dof_map, num_vars);
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
        for (unsigned int dof_i = 0; dof_i < n_var_dofs * num_vars; dof_i++)
        {
            for (unsigned int k = 0; k < num_vars; k++)
            {
                Fe(dof_i) += JxW[qp] * (phi_A[dof_i][qp](k));
            }
        }
    }
    // libmesh_error_msg("build_solid97 is not implemented yet.");
}
// Temporary placeholder for build_solid236
void ElementBuilder::build_solid236(const Elem *elem, const DofMap &dof_map)
{
    libmesh_error_msg("build_solid236 is not implemented yet.");
}

// Temporary placeholder for build_solid236_cond
void ElementBuilder::build_solid236_cond(const Elem *elem, const DofMap &dof_map)
{
    libmesh_error_msg("build_solid236_cond is not implemented yet.");
}

// Temporary placeholder for build_plane13
void ElementBuilder::build_plane13(const Elem *elem, const DofMap &dof_map)
{
    // Stop logging the shape function initialization.
    // If you forget to stop logging an event the PerfLog
    // object will probably catch the error and abort.
    const unsigned int num_vars = 1;
    this->initialize(elem, dof_map, num_vars, false);
    // Now we will build the element matrix.  This involves
    // a double loop to integrate the test functions (i) against
    // the trial functions (j).
    //
    // We have split the numeric integration into two loops
    // so that we can log the matrix and right-hand-side
    // computation separately.
    //
    // Now start logging the element matrix computation
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
        for (unsigned int i = 0; i != n_var_dofs; i++)
        {
            for (unsigned int j = 0; j != n_var_dofs; j++)
                Ke(i, j) += JxW[qp] * (dphi_V[i][qp] * dphi_V[j][qp]) / mu_matrix[0][0];
        }
    }
}
// Temporary placeholder for build_plane13
void ElementBuilder::build_plane13_rhs(const Elem *elem, const DofMap &dof_map)
{
    // Stop logging the shape function initialization.
    // If you forget to stop logging an event the PerfLog
    // object will probably catch the error and abort.
    const unsigned int num_vars = 1;
    this->initialize_rhs(elem, dof_map, num_vars);
    // Now we will build the element matrix.  This involves
    // a double loop to integrate the test functions (i) against
    // the trial functions (j).
    //
    // We have split the numeric integration into two loops
    // so that we can log the matrix and right-hand-side
    // computation separately.
    //
    // Now start logging the element matrix computation

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
        for (unsigned int i = 0; i != n_var_dofs; i++)
        {
            Fe(i) += JxW[qp] * phi_V[i][qp];
        }
    }
}
// Temporary placeholder for build_plane13_cond
void ElementBuilder::build_plane13_cond(const Elem *elem, const DofMap &dof_map)
{
    // Numeric ids corresponding to each variable in the system
    const unsigned int num_vars = 1;
    this->initialize(elem, dof_map, num_vars, false);

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
        // Real js[3];
        // jsArray->GetTuple(el_id, js);
        // // Here check that js is 0

        for (unsigned int dof_i = 0; dof_i < n_var_dofs; dof_i++)
            for (unsigned int dof_j = 0; dof_j < n_var_dofs; dof_j++)
            {
                // Here you will put the time factor as well
                Ke(dof_i, dof_j) += JxW[qp] * (dphi_V[dof_i][qp] * dphi_V[dof_j][qp]) / mu_matrix[0][0];

                Ce(dof_i, dof_j) += JxW[qp] / rs_matrix[0][0] * (phi_V[dof_i][qp] * phi_V[dof_j][qp]);
            }
    }
    // libmesh_error_msg("build_solid97 is not implemented yet.");
}

void ElementBuilder::build_element_matrix(const Elem *elem, const DofMap &dof_map, subdomain_id_type subdomain_id, const std::string &formulation, const unsigned int dim, const std::string &fe_family)
{
    this->n_var_dofs = elem->n_nodes();
    if (subdomain_id == 1 && formulation == "AV" && dim == 3 && fe_family == "LAGRANGE")
    {
        this->build_solid97_cond(elem, dof_map);
    }
    if (subdomain_id == 0 && formulation == "AV" && dim == 3 && fe_family == "LAGRANGE")
    {
        this->build_solid97(elem, dof_map);
    }
    if (subdomain_id == 1 && formulation == "AV" && fe_family == "NEDELEC_ONE")
    {
        this->build_solid236_cond(elem, dof_map);
    }
    if (subdomain_id == 0 && formulation == "AV" && fe_family == "NEDELEC_ONE")
    {
        this->build_solid236(elem, dof_map);
    }
    if (subdomain_id == 1 && formulation == "AV" && dim == 2)
    {
        this->build_plane13_cond(elem, dof_map);
    }
    if (subdomain_id == 0 && formulation == "AV" && dim == 2)
    {
        this->build_plane13(elem, dof_map);
    }
}
void ElementBuilder::build_element_rhs(const Elem *elem, const DofMap &dof_map, subdomain_id_type subdomain_id, const std::string &formulation, const unsigned int dim, const std::string &fe_family)
{
    if (subdomain_id == 0 && formulation == "AV" && dim == 3 && fe_family == "LAGRANGE")
    {
        this->build_solid97_rhs(elem, dof_map);
    }
    if (subdomain_id == 0 && formulation == "AV" && dim == 2)
    {
        this->build_plane13_rhs(elem, dof_map);
    }
}
// Method to fill mu_matrix
void ElementBuilder::fill_mu_matrix(dof_id_type elem_id)
{
    // get ansys_material_type from the materials map

    const unsigned int material_number = mat_array->GetTuple1(elem_id); // Retrieve the material number using GetTuple1
    const MaterialsProperties &mat_properties = materials.at(material_number);

    mu_matrix[0][0] = 4 * libMesh::pi * 1e-7 * mat_properties.murx;
    mu_matrix[1][1] = 4 * libMesh::pi * 1e-7 * mat_properties.mury;
    mu_matrix[2][2] = 4 * libMesh::pi * 1e-7 * mat_properties.murz;
}

// Method to fill rs_matrix
void ElementBuilder::fill_rs_matrix(dof_id_type elem_id)
{
    const unsigned int material_number = mat_array->GetTuple1(elem_id); // Retrieve the material number using GetTuple1
    const MaterialsProperties &mat_properties = materials.at(material_number);

    rs_matrix[0][0] = mat_properties.rsvx;
    rs_matrix[1][1] = mat_properties.rsvy;
    rs_matrix[2][2] = mat_properties.rsvz;
}

void ElementBuilder::compute_material_properties(dof_id_type elem_id, subdomain_id_type subdomain_id)
{
    // Call fill_mu_matrix and fill_rs_matrix to compute the material properties
    this->fill_mu_matrix(elem_id);
    if (subdomain_id == 1)
        this->fill_rs_matrix(elem_id);
}