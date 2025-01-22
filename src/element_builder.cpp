#include "../include/element_builder.h"
#include <cmath>
// Constructor
ElementBuilder::ElementBuilder(DenseMatrix<Number> &Ke, DenseMatrix<Number> &Ce, DenseVector<Number> &Fe, MaterialsMap &materials, const std::unique_ptr<FEVectorBase> &fe_vec, const std::unique_ptr<FEBase> &fe_scalar, vtkCellData *cell_data, const TimeStep &ts)
    : Ke(Ke), Ce(Ce), Fe(Fe), materials(materials),
      mu_matrix(3, std::vector<double>(3, 0.0)),
      rs_matrix(3, std::vector<double>(3, 0.0)),
      JxW_A(fe_vec->get_JxW()),
      JxW_V(fe_scalar->get_JxW()),
      phi_A(fe_vec->get_phi()),
      dphi_A(fe_vec->get_dphi()),
      curl_phi_A(fe_vec->get_curl_phi()),
      //   div_phi_A(fe_vec->get_div_phi()),
      phi_V(fe_scalar->get_phi()),
      dphi_V(fe_scalar->get_dphi()),
      ts(ts),
      mat_array(nullptr),
      n_var_dofs(0)
{
    std::string array_name = "ansys_material_type";
    mat_array = cell_data->GetArray(array_name.c_str());
    if (mat_array == nullptr)
    {
        libmesh_error_msg("Error: 'ansys_material_type' array not found in vtkCellData.");
    }
}

// Method
void ElementBuilder::build_solid97_cond(const Elem *elem, const DofMap &dof_map)
{
    // Numeric ids corresponding to each variable in the system
    unsigned int num_vars = 3;
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
    for (unsigned int qp = 0; qp < phi_A[0].size(); qp++)
    {
        for (unsigned int dof_i = 0; dof_i < n_var_dofs * num_vars; dof_i++)
        {
            for (unsigned int dof_j = 0; dof_j < n_var_dofs * num_vars; dof_j++)
            {
                for (unsigned int k = 0; k < 3; k++)
                {

                    Ke(dof_i, dof_j) += JxW_A[qp] / mu_matrix[k][k] * (dphi_A[dof_j][qp].row(k) * dphi_A[dof_i][qp].row(k));
                    Ce(dof_i, dof_j) += JxW_A[qp] / rs_matrix[k][k] * (phi_A[dof_i][qp](k) * phi_A[dof_j][qp](k));
                    if (dof_j < n_var_dofs)
                    {
                        CAV_value = JxW_A[qp] / rs_matrix[k][k] * phi_A[dof_i][qp](k) * dphi_V[dof_j][qp](k);
                        Ce(dof_i, num_vars * n_var_dofs + dof_j) += CAV_value;
                        Ce(num_vars * n_var_dofs + dof_j, dof_i) += CAV_value;
                        if (dof_i < n_var_dofs)
                        {
                            Ce(num_vars * n_var_dofs + dof_i, num_vars * n_var_dofs + dof_j) += JxW_A[qp] * dphi_V[dof_i][qp](k) / rs_matrix[k][k] * dphi_V[dof_j][qp](k);
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
    // Now we will build the element matrix.  This involves
    // a double loop to integrate the test functions (i) against
    // the trial functions (j).
    //
    // We have split the numeric integration into two loops
    // so that we can log the matrix and right-hand-side
    // computation separately.
    //
    // Now start logging the element matrix computation
    for (unsigned int qp = 0; qp < phi_A[0].size(); qp++)
    {

        for (unsigned int dof_i = 0; dof_i < n_var_dofs * num_vars; dof_i++)
        {
            for (unsigned int dof_j = 0; dof_j < n_var_dofs * num_vars; dof_j++)
            {
                for (unsigned int k = 0; k < 3; k++)
                {
                    Ke(dof_i, dof_j) += JxW_A[qp] / mu_matrix[k][k] * (dphi_A[dof_j][qp].row(k) * dphi_A[dof_i][qp].row(k));
                }
            }
        }
    }
    // libmesh_error_msg("build_solid97 is not implemented yet.");
}
void ElementBuilder::build_solid97_rhs(const Elem *elem, const DofMap &dof_map)
{
    const unsigned int num_vars = 3;
    for (unsigned int qp = 0; qp < phi_A[0].size(); qp++)
    {
        for (unsigned int dof_i = 0; dof_i < n_var_dofs * num_vars; dof_i++)
        {
            for (unsigned int k = 0; k < num_vars; k++)
            {
                Fe(dof_i) += JxW_A[qp] * (phi_A[dof_i][qp](k));
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
    // Now we will build the element matrix.  This involves
    // a double loop to integrate the test functions (i) against
    // the trial functions (j).
    //
    // We have split the numeric integration into two loops
    // so that we can log the matrix and right-hand-side
    // computation separately.
    //
    // Now start logging the element matrix computation
    for (unsigned int qp = 0; qp < phi_V[0].size(); qp++)
    {
        for (unsigned int i = 0; i != n_var_dofs; i++)
        {
            for (unsigned int j = 0; j != n_var_dofs; j++)
                Ke(i, j) += JxW_V[qp] * (dphi_V[i][qp] * dphi_V[j][qp]) / mu_matrix[0][0];
        }
    }
}
// Temporary placeholder for build_plane13
void ElementBuilder::build_plane13_rhs(const Elem *elem, const DofMap &dof_map)
{
    // Stop logging the shape function initialization.
    // If you forget to stop logging an event the PerfLog
    // object will probably catch the error and abort.
    // Now we will build the element matrix.  This involves
    // a double loop to integrate the test functions (i) against
    // the trial functions (j).
    //
    // We have split the numeric integration into two loops
    // so that we can log the matrix and right-hand-side
    // computation separately.
    //
    // Now start logging the element matrix computation

    for (unsigned int qp = 0; qp < phi_V[0].size(); qp++)
    {
        for (unsigned int i = 0; i != n_var_dofs; i++)
        {
            Fe(i) += JxW_V[qp] * phi_V[i][qp];
        }
    }
}
// Temporary placeholder for build_plane13_cond
void ElementBuilder::build_plane13_cond(const Elem *elem, const DofMap &dof_map)
{

    for (unsigned int qp = 0; qp < phi_V[0].size(); qp++)
    {
        // Real js[3];
        // jsArray->GetTuple(el_id, js);
        // // Here check that js is 0

        for (unsigned int dof_i = 0; dof_i < n_var_dofs; dof_i++)
            for (unsigned int dof_j = 0; dof_j < n_var_dofs; dof_j++)
            {
                // Here you will put the time factor as well
                Ke(dof_i, dof_j) += JxW_V[qp] * (dphi_V[dof_i][qp] * dphi_V[dof_j][qp]) / mu_matrix[0][0];

                Ce(dof_i, dof_j) += JxW_V[qp] / rs_matrix[0][0] * (phi_V[dof_i][qp] * phi_V[dof_j][qp]);
            }
    }
    // libmesh_error_msg("build_solid97 is not implemented yet.");
}

void ElementBuilder::build_element_matrix(const Elem *elem, const DofMap &dof_map, bool conductive, const std::string &formulation, const unsigned int dim, const std::string &fe_family)
{
    this->n_var_dofs = elem->n_nodes();
    if (conductive && formulation == "AV" && dim == 3 && fe_family == "LAGRANGE")
    {
        this->build_solid97_cond(elem, dof_map);
    }
    if (!conductive && formulation == "AV" && dim == 3 && fe_family == "LAGRANGE")
    {
        this->build_solid97(elem, dof_map);
    }
    if (conductive && formulation == "AV" && fe_family == "NEDELEC_ONE")
    {
        this->build_solid236_cond(elem, dof_map);
    }
    if (!conductive && formulation == "AV" && fe_family == "NEDELEC_ONE")
    {
        this->build_solid236(elem, dof_map);
    }
    if (conductive && formulation == "AV" && dim == 2)
    {
        this->build_plane13_cond(elem, dof_map);
    }
    if (!conductive && formulation == "AV" && dim == 2)
    {
        this->build_plane13(elem, dof_map);
    }
}
void ElementBuilder::build_element_rhs(const Elem *elem, const DofMap &dof_map, bool conductive, const std::string &formulation, const unsigned int dim, const std::string &fe_family)
{
    this->n_var_dofs = elem->n_nodes();

    if (!conductive && formulation == "AV" && dim == 3 && fe_family == "LAGRANGE")
    {
        this->build_solid97_rhs(elem, dof_map);
    }
    if (!conductive && formulation == "AV" && dim == 2)
    {
        this->build_plane13_rhs(elem, dof_map);
    }
}
// Method to fill mu_matrix
void ElementBuilder::fill_mu_matrix(dof_id_type local_elem_id)
{
    // get ansys_material_type from the materials map

    const unsigned int material_number = mat_array->GetTuple1(local_elem_id); // Retrieve the material number using GetTuple1
    const MaterialsProperties &mat_properties = materials.at(material_number);

    mu_matrix[0][0] = 4 * libMesh::pi * 1e-7 * mat_properties.murx;
    mu_matrix[1][1] = 4 * libMesh::pi * 1e-7 * mat_properties.mury;
    mu_matrix[2][2] = 4 * libMesh::pi * 1e-7 * mat_properties.murz;
}

// Method to fill rs_matrix
void ElementBuilder::fill_rs_matrix(dof_id_type local_elem_id)
{
    const unsigned int material_number = mat_array->GetTuple1(local_elem_id); // Retrieve the material number using GetTuple1
    const MaterialsProperties &mat_properties = materials.at(material_number);

    rs_matrix[0][0] = mat_properties.rsvx;
    rs_matrix[1][1] = mat_properties.rsvy;
    rs_matrix[2][2] = mat_properties.rsvz;
}

void ElementBuilder::compute_material_properties(dof_id_type local_elem_id, bool conductive)
{
    // Call fill_mu_matrix and fill_rs_matrix to compute the material properties
    this->fill_mu_matrix(local_elem_id);
    if (conductive)
        this->fill_rs_matrix(local_elem_id);
}

void ElementBuilder::build_element_matrix_derived(DenseMatrix<Number> &AAe, DenseMatrix<Number> &Me, DenseMatrix<Number> &JAe, DenseMatrix<Number> &JVe, const Elem *elem, const DofMap &dof_map, bool conductive, const std::string &formulation, const unsigned int dim, const std::string &fe_family, unsigned int n_A_dofs)
{

    if (dim == 2)
    {
        for (unsigned int qp = 0; qp < phi_V[0].size(); qp++)
        {
            for (unsigned int i = 0; i < 2 * n_A_dofs; i++)
            {
                if (conductive)
                    if (i < n_var_dofs)
                    {
                        // JAe(0, i) -= JxW_V[qp] * phi_V[i][qp] / phi_V[0].size();
                        JAe(0, i) -= phi_V[i][qp] / phi_V[0].size();
                    }
                for (unsigned int j = 0; j < 2 * n_A_dofs; j++)
                {
                    if (j < n_A_dofs)
                    {
                        AAe(i, j) += JxW_V[qp] * std::pow(-1, i) * (phi_V[i / 2][qp] * dphi_V[j][qp](1 - (i % 2)));
                    }
                    if ((i + j) % 2 == 0)
                        Me(i, j) += JxW_V[qp] * (phi_V[i / 2][qp] * phi_V[j / 2][qp]);
                }
            }
        }
    }
    // if (dim == 2)
    // {
    //     for (unsigned int qp = 0; qp < phi_V[0].size(); qp++)
    //     {
    //         for (unsigned int i = 0; i < n_A_dofs; i++)
    //         {
    //             for (unsigned int j = 0; j < n_A_dofs; j++)
    //             {
    //                 AAe(i, j) += JxW_V[qp] * (phi_V[i][qp] * dphi_V[j][qp](1));
    //                 Me(i, j) += JxW_V[qp] * (phi_V[i][qp] * phi_V[j][qp]);
    //             }
    //         }
    //     }
    // }
    else if (dim == 3)
    {
        for (unsigned int qp = 0; qp < phi_A[0].size(); qp++)
        {
            for (unsigned int i = 0; i < n_A_dofs; i++)
            {
                if (conductive)
                    for (unsigned int k = 0; k < 3; k++)
                    {
                        // JAe(k, i) += JxW_A[qp] * phi_A[i][qp](k);
                        JAe(k, i) -= phi_A[i][qp](k) / phi_A[0].size();
                        // JAe(k, i) -= phi_A[i][qp](k);
                        if (i < n_var_dofs)
                            // JVe(k, i) += JxW_A[qp] * dphi_V[i][qp](k);
                            JVe(k, i) += dphi_V[i][qp](k) / dphi_V[0].size();
                        // JVe(k, i) += JxW_A[qp] * dphi_V[i][qp](k) / dphi_V[0].size();
                    }

                for (unsigned int j = 0; j < n_A_dofs; j++)
                {
                    AAe(i, j) += JxW_A[qp] * (phi_A[i][qp] * curl_phi_A[j][qp]);
                    Me(i, j) += JxW_A[qp] * (phi_A[i][qp] * phi_A[j][qp]);
                }
            }
        }
    }
}

MaterialBuilder::MaterialBuilder(MaterialsMap &materials, vtkCellData *cell_data)
    : materials(materials)
{
    std::string array_name = "ansys_material_type";
    mat_array = cell_data->GetArray(array_name.c_str());
    if (mat_array == nullptr)
    {
        libmesh_error_msg("Error: 'ansys_material_type' array not found in vtkCellData.");
    }
    mu_matrix.resize(3);
    rs_matrix.resize(3);
}
// Method to fill mu_matrix
void MaterialBuilder::fill_mu_matrix(dof_id_type local_elem_id)
{
    // get ansys_material_type from the materials map

    const unsigned int material_number = mat_array->GetTuple1(local_elem_id); // Retrieve the material number using GetTuple1
    const MaterialsProperties &mat_properties = materials.at(material_number);

    mu_matrix(0) = 4 * libMesh::pi * 1e-7 * mat_properties.murx;
    mu_matrix(1) = 4 * libMesh::pi * 1e-7 * mat_properties.mury;
    mu_matrix(2) = 4 * libMesh::pi * 1e-7 * mat_properties.murz;
}

// Method to fill rs_matrix
void MaterialBuilder::fill_rs_matrix(dof_id_type local_elem_id)
{
    const unsigned int material_number = mat_array->GetTuple1(local_elem_id); // Retrieve the material number using GetTuple1
    const MaterialsProperties &mat_properties = materials.at(material_number);

    rs_matrix(0) = mat_properties.rsvx;
    rs_matrix(1) = mat_properties.rsvy;
    rs_matrix(2) = mat_properties.rsvz;
}

void MaterialBuilder::compute_material_properties(dof_id_type local_elem_id, bool conductive)
{
    // Call fill_mu_matrix and fill_rs_matrix to compute the material properties
    this->fill_mu_matrix(local_elem_id);
    if (conductive)
        this->fill_rs_matrix(local_elem_id);
}