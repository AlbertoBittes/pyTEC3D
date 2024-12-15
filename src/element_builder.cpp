#include "../include/element_builder.h"

// Constructor
ElementBuilder::ElementBuilder(DenseMatrix<Number> &Ke, DenseMatrix<Number> &Ce, DenseVector<Number> &Fe, MaterialsMap &materials, const std::unique_ptr<FEBase> &fe, vtkCellData *cell_data, TimeStep &ts, const DofMap &dof_map, QGauss &qrule)
    : Ke(Ke), Ce(Ce), Fe(Fe), materials(materials),
      F{DenseSubVector<Number>(Fe), DenseSubVector<Number>(Fe), DenseSubVector<Number>(Fe), DenseSubVector<Number>(Fe)},
      KAA{{DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
          {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
          {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)}},
      KAV{DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
      KVA{DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
      KVV(DenseSubMatrix<Number>(Ke)),
      CAA{{DenseSubMatrix<Number>(Ce), DenseSubMatrix<Number>(Ce), DenseSubMatrix<Number>(Ce)},
          {DenseSubMatrix<Number>(Ce), DenseSubMatrix<Number>(Ce), DenseSubMatrix<Number>(Ce)},
          {DenseSubMatrix<Number>(Ce), DenseSubMatrix<Number>(Ce), DenseSubMatrix<Number>(Ce)}},
      CAV{DenseSubMatrix<Number>(Ce), DenseSubMatrix<Number>(Ce), DenseSubMatrix<Number>(Ce)},
      CVA{DenseSubMatrix<Number>(Ce), DenseSubMatrix<Number>(Ce), DenseSubMatrix<Number>(Ce)},
      CVV(DenseSubMatrix<Number>(Ce)),
      mu_matrix(3, std::vector<double>(3, 0.0)),
      rs_matrix(3, std::vector<double>(3, 0.0)),
      JxW(fe->get_JxW()),
      phi(fe->get_phi()),
      dphi(fe->get_dphi()),
      ts(ts),
      dof_map(dof_map),
      qrule(qrule)
{
    std::string array_name = "Js_" + std::to_string(ts.get_step());
    js_array = cell_data->GetArray(array_name.c_str());
    array_name = "ansys_material_type";
    mat_array = cell_data->GetArray(array_name.c_str());
    if (!js_array)
    {
        // Get the value from the "Js" array using the element ID
        libmesh_error_msg("Error: Js_" + std::to_string(ts.get_step()) + " array not found in vtkCellData.");
    }
}

void ElementBuilder::initialize(const Elem *elem, DofMap &dof_map, unsigned int num_vars, bool volt)
{
    std::vector<std::vector<dof_id_type>> dof_indices_var(num_vars);

    for (unsigned int var = 0; var < num_vars; var++)
        dof_map.dof_indices(elem, dof_indices_var[var], var);

    n_var_dofs = dof_indices_var[0].size();

    for (unsigned int var_i = 0; var_i < num_vars - 1; var_i++)
        for (unsigned int var_j = 0; var_j < num_vars - 1; var_j++)
        {
            KAA[var_i][var_j].reposition(var_i * n_var_dofs, var_j * n_var_dofs, n_var_dofs, n_var_dofs);
            CAA[var_i][var_j].reposition(var_i * n_var_dofs, var_j * n_var_dofs, n_var_dofs, n_var_dofs);
        }

    if (volt)
    {
        for (unsigned int var_i = 0; var_i < num_vars - 1; var_i++)
        {
            KAV[var_i].reposition(var_i * n_var_dofs, (num_vars - 1) * n_var_dofs, n_var_dofs, n_var_dofs);
            KVA[var_i].reposition((num_vars - 1) * n_var_dofs, var_i * n_var_dofs, n_var_dofs, n_var_dofs);
            CAV[var_i].reposition(var_i * n_var_dofs, (num_vars - 1) * n_var_dofs, n_var_dofs, n_var_dofs);
            CVA[var_i].reposition((num_vars - 1) * n_var_dofs, var_i * n_var_dofs, n_var_dofs, n_var_dofs);
        }
        for (unsigned int var = 0; var < num_vars; var++)
            F[var].reposition(var * n_var_dofs, n_var_dofs);

        KVV.reposition((num_vars - 1) * n_var_dofs, (num_vars - 1) * n_var_dofs, n_var_dofs, n_var_dofs);
        CVV.reposition((num_vars - 1) * n_var_dofs, (num_vars - 1) * n_var_dofs, n_var_dofs, n_var_dofs);
    }
}
// Method
void ElementBuilder::build_solid97_cond(const Elem *elem, DofMap &dof_map)
{
    // Numeric ids corresponding to each variable in the system
    unsigned int num_vars = 4;
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
        // Real js[3];
        // jsArray->GetTuple(el_id, js);
        // // Here check that js is 0

        for (unsigned int dof_i = 0; dof_i < n_var_dofs; dof_i++)
            for (unsigned int dof_j = 0; dof_j < n_var_dofs; dof_j++)
                for (unsigned int k = 0; k < 3; k++)
                {
                    // Here you will put the time factor as well
                    CAV_value = JxW[qp] / rs_matrix[k][k] * phi[dof_i][qp] * dphi[dof_j][qp](k); // Not sure
                    CAV[k](dof_i, dof_j) += CAV_value;
                    CVA[k](dof_j, dof_i) += CAV_value;
                    CVV(dof_i, dof_j) += JxW[qp] * dphi[dof_i][qp](k) / rs_matrix[k][k] * dphi[dof_j][qp](k);
                    KAA[k][k](dof_i, dof_j) += JxW[qp] / (4 * libMesh::pi * 1e-7 * mu_matrix[k][k]) *
                                               (dphi[dof_i][qp] * dphi[dof_j][qp]);
                    CAA[k][k](dof_i, dof_j) += JxW[qp] / rs_matrix[k][k] *
                                               (phi[dof_i][qp] * phi[dof_j][qp]);
                }

        Ce.scale(1 / ts.get_substeps()[0].get_dt());
        Ke += Ce;
        Ce.vector_mult(Fe, old_solution_values);
    }
}

void ElementBuilder::build_solid97(const Elem *elem, DofMap &dof_map)
{
    const unsigned int num_vars = 3;
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
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
        Real js[3];
        js_array->GetTuple(elem->id(), js);
        // // Here check that js is 0

        for (unsigned int dof_i = 0; dof_i < n_var_dofs; dof_i++)
            for (unsigned int dof_j = 0; dof_j < n_var_dofs; dof_j++)
                for (unsigned int k = 0; k < 3; k++)
                {
                    KAA[k][k](dof_i, dof_j) += JxW[qp] / (4 * libMesh::pi * 1e-7 * mu_matrix[k][k]) *
                                               (dphi[dof_i][qp] * dphi[dof_j][qp]);

                    F[k](dof_i) += JxW[qp] * js[k] * phi[dof_i][qp];
                }
    }
    // libmesh_error_msg("build_solid97 is not implemented yet.");
}
// Method to fill mu_matrix
void ElementBuilder::fill_mu_matrix(dof_id_type elem_id)
{
    const unsigned int material_number = static_cast<unsigned int>(materials.at(elem_id).murx); // Retrieve the material number using GetTuple1
    const MaterialsProperties &mat_properties = materials.at(material_number);

    mu_matrix[0][0] = mat_properties.murx;
    mu_matrix[1][1] = mat_properties.mury;
    mu_matrix[2][2] = mat_properties.murz;
}

// Method to fill rs_matrix
void ElementBuilder::fill_rs_matrix(dof_id_type elem_id)
{
    const unsigned int material_number = static_cast<unsigned int>(materials.at(elem_id).rsvx); // Retrieve the material number using GetTuple1
    const MaterialsProperties &mat_properties = materials.at(material_number);

    rs_matrix[0][0] = mat_properties.rsvx;
    rs_matrix[1][1] = mat_properties.rsvy;
    rs_matrix[2][2] = mat_properties.rsvz;
}

void ElementBuilder::compute_material_properties(dof_id_type elem_id, const unsigned int subdomain_id)
{
    // Call fill_mu_matrix and fill_rs_matrix to compute the material properties
    this->fill_mu_matrix(elem_id);
    if (subdomain_id == 1)
        this->fill_rs_matrix(elem_id);
}