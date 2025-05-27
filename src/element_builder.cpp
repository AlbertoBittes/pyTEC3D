#include "../include/element_builder.h"
#include "libmesh/equation_systems.h"
#include "libmesh/string_to_enum.h"
// Define the Finite Element object.
#include "libmesh/fe.h"

// Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"
#include "../include/session.h" // Ensure the Session class is included
#include <cmath>
#include "../include/element_types.h" // Include the header defining ElementType
#include "../include/loads_object.h"
#include "../include/assembler.h"

void ElementBuilder::buildVectors(const std::map<std::string, std::tuple<unsigned int, unsigned int>> &sizes_map,
                                  std::vector<DenseVector<Number>> &elem_vectors_vector)
{
    elem_vectors_vector.resize(sizes_map.size());
    for (const auto &[key, value] : sizes_map)
    {
        unsigned int index = std::distance(sizes_map.begin(), sizes_map.find(key));
        elem_vectors_vector[index] = DenseVector<Number>(std::get<1>(value));
    }
}

void ElementBuilder::buildMatrices(const std::map<std::string, std::tuple<unsigned int, unsigned int, unsigned int>> &sizes_map,
                                   std::vector<DenseMatrix<Number>> &elem_matrices_vector)
{
    elem_matrices_vector.resize(sizes_map.size());
    for (const auto &[key, value] : sizes_map)
    {
        unsigned int index = std::distance(sizes_map.begin(), sizes_map.find(key));
        elem_matrices_vector[index] = DenseMatrix<Number>(std::get<1>(value), std::get<2>(value));
    }
}

void ElementBuilder::writeElementMatrices(std::ofstream &outFile,
                                          std::vector<DenseMatrix<Number>> &elem_matrices_vector)
{
    // Write the number of containers as an unsigned int
    unsigned int numContainers = static_cast<unsigned int>(elem_matrices_vector.size());
    outFile.write(reinterpret_cast<const char *>(&numContainers), sizeof(numContainers));

    for (const auto &elem_matrix : elem_matrices_vector)
    {
        unsigned int elem_id = this->elem->id();
        // Write the element ID as an unsigned int
        outFile.write(reinterpret_cast<const char *>(&elem_id), sizeof(elem_id));

        const auto &values = elem_matrix.get_values();

        // Write the size of the values vector as an unsigned int
        // unsigned int valuesSize = static_cast<unsigned int>(values.size());
        unsigned int rows = elem_matrix.m();
        unsigned int cols = elem_matrix.n();

        // Write the number of rows and columns as unsigned ints
        outFile.write(reinterpret_cast<const char *>(&rows), sizeof(rows));
        outFile.write(reinterpret_cast<const char *>(&cols), sizeof(cols));

        // Write the actual values
        outFile.write(reinterpret_cast<const char *>(values.data()), values.size() * sizeof(Number));
    }
}

void ElementBuilder::readElementMatrices(std::ifstream &inFile, std::vector<DenseMatrix<Number>> &elem_matrices_vector)
{
    // Read the number of containers as an unsigned int
    unsigned int numContainers;
    inFile.read(reinterpret_cast<char *>(&numContainers), sizeof(numContainers));

    // Resize the containers vector to match the number of containers
    elem_matrices_vector.resize(numContainers);

    for (auto &elem_matrix : elem_matrices_vector)
    {
        // Read the element ID as an unsigned int
        unsigned int elem_id;
        inFile.read(reinterpret_cast<char *>(&elem_id), sizeof(elem_id));
        // Read the size of the values vector as an unsigned int
        unsigned int rows, cols;
        inFile.read(reinterpret_cast<char *>(&rows), sizeof(rows));
        inFile.read(reinterpret_cast<char *>(&cols), sizeof(cols));

        // Resize the container's values vector to match the size
        elem_matrix.resize(rows, cols);
        auto &values = elem_matrix.get_values();

        // Read the actual values into the container's values vector
        inFile.read(reinterpret_cast<char *>(values.data()), values.size() * sizeof(Number));
    }
}

// void ElementBuilder::readElementData(std::ifstream &inFile, std::vector<DenseMatrix<Number>> &elem_matrices_vector, std::vector<DenseVector<Number>> &elem_vectors_vector)
// {
//     // Read the matrices and vectors
//     readElementMatrices(inFile, elem_matrices_vector);
//     readElementVectors(inFile, elem_vectors_vector);
// }

// void ElementBuilder::writeElementData(std::ofstream &outFile, const std::vector<DenseMatrix<Number>> &elem_matrices_vector, const std::vector<DenseVector<Number>> &elem_vectors_vector)
// {
//     // Write the matrices and vectors
//     writeElementMatrices(outFile, elem_matrices_vector);
//     writeElementVectors(outFile, elem_vectors_vector);
// }

void ElementBuilder::writeElementVectors(std::ofstream &outFile, std::vector<DenseVector<Number>> &elem_vectors_vector)
{
    // Write the number of containers as an unsigned int
    unsigned int numContainers = static_cast<unsigned int>(elem_vectors_vector.size());
    outFile.write(reinterpret_cast<const char *>(&numContainers), sizeof(numContainers));

    for (const auto &elem_vector : elem_vectors_vector)
    {
        const auto &values = elem_vector.get_values();

        unsigned int elem_id = this->elem->id();
        // Write the element ID as an unsigned int
        outFile.write(reinterpret_cast<const char *>(&elem_id), sizeof(elem_id));
        // Write the size of the values vector as an unsigned int
        unsigned int valuesSize = elem_vector.size();
        outFile.write(reinterpret_cast<const char *>(&valuesSize), sizeof(valuesSize));

        // Write the actual values
        outFile.write(reinterpret_cast<const char *>(values.data()), values.size() * sizeof(Number));
    }
}

void ElementBuilder::readElementVectors(std::ifstream &inFile, std::vector<DenseVector<Number>> &elem_vectors_vector)
{
    // Read the number of containers as an unsigned int
    unsigned int numContainers;
    inFile.read(reinterpret_cast<char *>(&numContainers), sizeof(numContainers));

    // Resize the containers vector to match the number of containers
    elem_vectors_vector.resize(numContainers);

    for (auto &elem_vector : elem_vectors_vector)
    {
        // Read the element ID as an unsigned int
        unsigned int elem_id;
        inFile.read(reinterpret_cast<char *>(&elem_id), sizeof(elem_id));
        // Read the size of the values vector as an unsigned int
        unsigned int valuesSize;
        inFile.read(reinterpret_cast<char *>(&valuesSize), sizeof(valuesSize));

        // Resize the container's values vector to match the size
        elem_vector.resize(valuesSize);
        auto &values = elem_vector.get_values();

        // Read the actual values into the container's values vector
        inFile.read(reinterpret_cast<char *>(values.data()), values.size() * sizeof(Number));
    }
}

// {0: "Ke", 1: "Ce", 2: "Me", 3: "AAe", 4: "JAe", 5: "JVe"}
void ElementBuilder::build_solid97_matrix()
{
    // Now we will build the element matrix.  This involves
    // a double loop to integrate the test functions (i) against
    // the trial functions (j).
    //
    // We have split the numeric integration into two loops
    // so that we can log the matrix and right-hand-side
    // computation separately.
    //
    // Now start logging the element matrix computation
    // Numeric ids corresponding to each variable in the system
    // Now we will build the element matrix.  This involves
    // a double loop to integrate the test functions (i) against
    // the trial functions (j).
    //
    // We have split the numeric integration into two loops
    // so that we can log the matrix and right-hand-side
    // computation separately.
    //
    // Now start logging the element matrix computation
    unsigned int n_dofs_var = element_type->get_fe()->n_shape_functions();
    unsigned int n_dofs_A = n_dofs_var * 3;
    unsigned int nq = element_type->get_fe()->n_quadrature_points();

    auto eq_sys = this->get_parent_object()->get_parent_object();
    auto load_step = eq_sys->parameters.get<std::shared_ptr<LoadStepper>>("load_stepper")->get_active_load_step();

    auto materials = load_step->get_materials();

    materials->compute_material_properties(this->material_num, element_type->is_conductive());
    auto mu_matrix = materials->get_mu_matrix();
    auto rs_matrix = materials->get_rs_matrix();

    if (element_type->is_conductive())
    {
        // Now start logging the element matrix computation
        Real CAV_value;

        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            {
                // matrices_derived[2](qp, dof_i) += scalar_shape_functions->phi[dof_i][qp];

                // for (unsigned int k = 0; k < 3; k++)
                // matrices_derived[3](qp + k * n_dofs_var, dof_i) += scalar_shape_functions->dphi[dof_i][qp](k);

                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    for (unsigned int l = 0; l < 3; l++)
                    {
                        matrices[1](dof_i + l * n_dofs_var, dof_j + l * n_dofs_var) += scalar_shape_functions->JxW[qp] / rs_matrix[l] * (scalar_shape_functions->phi[dof_i][qp] * scalar_shape_functions->phi[dof_j][qp]);
                        CAV_value = scalar_shape_functions->JxW[qp] / rs_matrix[l] * scalar_shape_functions->phi[dof_i][qp] * scalar_shape_functions->dphi[dof_j][qp](l);
                        matrices[1](dof_i + l * n_dofs_var, n_dofs_A + dof_j) += CAV_value;
                        matrices[1](n_dofs_A + dof_j, dof_i + l * n_dofs_var) += CAV_value;
                        matrices[1](n_dofs_A + dof_i, n_dofs_A + dof_j) += scalar_shape_functions->JxW[qp] * scalar_shape_functions->dphi[dof_i][qp](l) / rs_matrix[l] * scalar_shape_functions->dphi[dof_j][qp](l);
                        for (unsigned int k = 0; k < 3; k++)
                        {
                            matrices[0](dof_i + l * n_dofs_var, dof_j + l * n_dofs_var) += scalar_shape_functions->JxW[qp] / mu_matrix[k] * (scalar_shape_functions->dphi[dof_i][qp](k) * scalar_shape_functions->dphi[dof_j][qp](k));
                            // matrices_derived[1](dof_i + k * n_dofs_var, dof_j + l * n_dofs_var) += scalar_shape_functions->JxW[qp] * (vector_shape_functions->phi_vec[dof_i * 3 + k][qp] * vector_shape_functions->curl_phi_vec[dof_j * 3 + l][qp]);
                            // matrices_derived[0](dof_i + k * n_dofs_var, dof_j + l * n_dofs_var) += scalar_shape_functions->JxW[qp] * (vector_shape_functions->phi_vec[dof_i * 3 + k][qp] * vector_shape_functions->phi_vec[dof_j * 3 + l][qp]);
                        }
                    }
            }
    }
    else
        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    for (unsigned int l = 0; l < 3; l++)
                        for (unsigned int k = 0; k < 3; k++)
                        {
                            matrices[0](dof_i + l * n_dofs_var, dof_j + l * n_dofs_var) += scalar_shape_functions->JxW[qp] / mu_matrix[k] * (scalar_shape_functions->dphi[dof_i][qp](k) * scalar_shape_functions->dphi[dof_j][qp](k));
                            // matrices_derived[1](dof_i + k * n_dofs_var, dof_j + l * n_dofs_var) += scalar_shape_functions->JxW[qp] * (vector_shape_functions->phi_vec[dof_i * 3 + k][qp] * vector_shape_functions->curl_phi_vec[dof_j * 3 + l][qp]);
                            // matrices_derived[0](dof_i + k * n_dofs_var, dof_j + l * n_dofs_var) += scalar_shape_functions->JxW[qp] * (vector_shape_functions->phi_vec[dof_i * 3 + k][qp] * vector_shape_functions->phi_vec[dof_j * 3 + l][qp]);
                        }
}

// Temporary placeholder for build_solid236
void ElementBuilder::build_solid236_matrix()
{
    libmesh_error_msg("build_solid236 is not implemented yet.");
}

// Temporary placeholder for build_plane13
void ElementBuilder::build_plane13_matrix()
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
    // Now we will build the element matrix.  This involves
    // a double loop to integrate the test functions (i) against
    // the trial functions (j).
    //
    // We have split the numeric integration into two loops
    // so that we can log the matrix and right-hand-side
    // computation separately.
    //
    // Now start logging the element matrix computation
    // Numeric ids corresponding to each variable in the system
    // Now we will build the element matrix.  This involves
    // a double loop to integrate the test functions (i) against
    // the trial functions (j).
    //
    // We have split the numeric integration into two loops
    // so that we can log the matrix and right-hand-side
    // computation separately.
    //
    unsigned int n_dofs_var = element_type->get_fe()->n_shape_functions();
    unsigned int n_dofs_A = n_dofs_var * 3;
    unsigned int nq = element_type->get_fe()->n_quadrature_points();

    auto eq_sys = this->get_parent_object()->get_parent_object();
    auto load_step = eq_sys->parameters.get<std::shared_ptr<LoadStepper>>("load_stepper")->get_active_load_step();
    auto materials = load_step->get_materials();

    materials->compute_material_properties(this->material_num, element_type->is_conductive());
    auto mu_matrix = materials->get_mu_matrix();
    auto rs_matrix = materials->get_rs_matrix();

    if (element_type->is_conductive())
    {
        // Now start logging the element matrix computation
        for (unsigned int qp = 0; qp < nq; qp++)
        {
            for (unsigned int dof_i = 0; dof_i < 2 * n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < 2 * n_dofs_var; dof_j++)
                {
                    if (dof_j < n_dofs_var)
                        matrices_solution[1](dof_i, dof_j) += scalar_shape_functions->JxW[qp] * std::pow(-1, dof_i) * (scalar_shape_functions->phi[dof_i / 2][qp] * scalar_shape_functions->dphi[dof_j][qp](1 - (dof_i % 2)));

                    if ((dof_i + dof_j) % 2 == 0)
                        matrices_solution[0](dof_i, dof_j) += scalar_shape_functions->JxW[qp] * (scalar_shape_functions->phi[dof_i / 2][qp] * scalar_shape_functions->phi[dof_j / 2][qp]);
                }
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            {
                // matrices[4](0, i) -= scalar_shape_functions->JxW[qp] * scalar_shape_functions->phi[i][qp] / scalar_shape_functions->phi[0].size();
                matrices_solution[2](qp, dof_i) -= scalar_shape_functions->phi[dof_i][qp];

                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                {
                    matrices[0](dof_i, dof_j) += scalar_shape_functions->JxW[qp] * (scalar_shape_functions->dphi[dof_i][qp] * scalar_shape_functions->dphi[dof_j][qp]) / mu_matrix[0];
                    matrices[1](dof_i, dof_j) += scalar_shape_functions->JxW[qp] / rs_matrix[0] * (scalar_shape_functions->phi[dof_i][qp] * scalar_shape_functions->phi[dof_j][qp]);
                }
            }
        }
    }
    else
        for (unsigned int qp = 0; qp < nq; qp++)
        {
            for (unsigned int dof_i = 0; dof_i < 2 * n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < 2 * n_dofs_var; dof_j++)
                {
                    if (dof_j < n_dofs_var)
                        matrices[3](dof_i, dof_j) += scalar_shape_functions->JxW[qp] * std::pow(-1, dof_i) * (scalar_shape_functions->phi[dof_i / 2][qp] * scalar_shape_functions->dphi[dof_j][qp](1 - (dof_i % 2)));

                    if ((dof_i + dof_j) % 2 == 0)
                        matrices[2](dof_i, dof_j) += scalar_shape_functions->JxW[qp] * (scalar_shape_functions->phi[dof_i / 2][qp] * scalar_shape_functions->phi[dof_j / 2][qp]);
                }

            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    for (unsigned int l = 0; l < 3; l++)
                        for (unsigned int k = 0; k < 3; k++)
                            matrices[0](dof_i, dof_j) += scalar_shape_functions->JxW[qp] * (scalar_shape_functions->dphi[dof_i][qp] * scalar_shape_functions->dphi[dof_j][qp]) / mu_matrix[0];
        }

    for (unsigned int qp = 0; qp < nq; qp++)
        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            vectors[0](dof_i) += scalar_shape_functions->JxW[qp] * scalar_shape_functions->phi[dof_i][qp]; // * js_loads[elem->id()][2];
}

void ElementBuilder::build_element_matrix()
{
    auto etype = element_type->get_properties()[0];

    if (etype == 97)
    {
        this->build_solid97_matrix();
    }
    else if (etype == 236)
    {
        this->build_solid236_matrix();
    }
    else if (etype == 13)
    {
        this->build_plane13_matrix();
    }
    else
    {
        libmesh_error_msg("Element type not supported: " + etype);
    }
}

void ElementBuilder::build_solid97_rhs()
{
    // Now we will build the element matrix.  This involves
    // a double loop to integrate the test functions (i) against
    // the trial functions (j).
    //
    // We have split the numeric integration into two loops
    // so that we can log the matrix and right-hand-side
    // computation separately.
    //
    // Now start logging the element matrix computation
    // Numeric ids corresponding to each variable in the system
    // Now we will build the element matrix.  This involves
    // a double loop to integrate the test functions (i) against
    // the trial functions (j).
    //
    // We have split the numeric integration into two loops
    // so that we can log the matrix and right-hand-side
    // computation separately.
    //
    // Now start logging the element matrix computation

    unsigned int n_dofs_var = element_type->get_fe()->n_shape_functions();
    unsigned int n_dofs_A = n_dofs_var * 3;
    unsigned int nq = element_type->get_fe()->n_quadrature_points();

    std::string js_label = "Js";
    auto eq_sys = this->get_parent_object()->get_parent_object();
    auto load_step = eq_sys->parameters.get<std::shared_ptr<LoadStepper>>("load_stepper")->get_active_load_step();
    auto js_loads = (*load_step->get_loads()->get_data_map())[js_label];

    for (unsigned int qp = 0; qp < nq; qp++)
        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            for (unsigned int k = 0; k < 3; k++)
                vectors[0](dof_i + k * n_dofs_var) += scalar_shape_functions->JxW[qp] * scalar_shape_functions->phi[dof_i][qp] * js_loads[elem->id()][k];
}

// Temporary placeholder for build_solid236
void ElementBuilder::build_solid236_rhs()
{
    libmesh_error_msg("build_solid236 is not implemented yet.");
}

// Temporary placeholder for build_plane13
void ElementBuilder::build_plane13_rhs()
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
    // Now we will build the element matrix.  This involves
    // a double loop to integrate the test functions (i) against
    // the trial functions (j).
    //
    // We have split the numeric integration into two loops
    // so that we can log the matrix and right-hand-side
    // computation separately.
    //
    // Now start logging the element matrix computation
    // Numeric ids corresponding to each variable in the system
    // Now we will build the element matrix.  This involves
    // a double loop to integrate the test functions (i) against
    // the trial functions (j).
    //
    // We have split the numeric integration into two loops
    // so that we can log the matrix and right-hand-side
    // computation separately.
    //
    // Now start logging the element matrix computation

    unsigned int n_dofs_var = element_type->get_fe()->n_shape_functions();
    unsigned int n_dofs_A = n_dofs_var * 3;
    unsigned int nq = element_type->get_fe()->n_quadrature_points();

    std::string js_label = "Js";
    auto eq_sys = this->get_parent_object()->get_parent_object();
    auto load_step = eq_sys->parameters.get<std::shared_ptr<LoadStepper>>("load_stepper")->get_active_load_step();

    auto js_loads = (*load_step->get_loads()->get_data_map())[js_label];

    for (unsigned int qp = 0; qp < nq; qp++)
        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            vectors[0](dof_i) += scalar_shape_functions->JxW[qp] * scalar_shape_functions->phi[dof_i][qp] * js_loads[elem->id()][2];
}

void ElementBuilder::build_element_rhs()
{
    auto etype = element_type->get_properties()[0];

    if (etype == 97)
    {
        this->build_solid97_rhs();
    }
    else if (etype == 236)
    {
        this->build_solid236_rhs();
    }
    else if (etype == 13)
    {
        this->build_plane13_rhs();
    }
    else
    {
        libmesh_error_msg("Element type not supported: " + etype);
    }
}
// void ElementBuilder::apply_material_properties()
// {
//     if (element_type->properties.element_type == "SOLID97")
//     {
//         this->apply_material_properties_solid97();
//     }
//     else if (element_type->properties.element_type == "SOLID236")
//     {
//         this->apply_material_properties_solid236();
//     }
//     else if (element_type->properties.element_type == "PLANE13")
//     {
//         this->apply_material_properties_plane13();
//     }
//     else
//     {
//         libmesh_error_msg("Element type not supported: " + element_type->properties.element_type);
//     }
// }
