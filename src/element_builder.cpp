#include "element_builder.h"
#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/string_to_enum.h"
// Define the Finite Element object.
#include "libmesh/fe.h"

// Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"
#include "session.h" // Ensure the Session class is included
#include <cmath>
#include "element_types.h" // Include the header defining ElementType
#include "loads_object.h"
#include "assembler.h"

void ElementBuilder::build_vectors(std::vector<std::vector<unsigned int>> *sizes_map, unsigned int vectors_set)
{
    this->vectors[vectors_set].resize(sizes_map->size());
    for (unsigned int i = 0; i < sizes_map->size(); ++i)
    {
        const auto &size_pair = (*sizes_map)[i];
        unsigned int rows = size_pair[0];
        this->vectors[vectors_set][i].resize(rows);
    }
}

void ElementBuilder::build_matrices(std::vector<std::vector<unsigned int>> *sizes_map, unsigned int matrices_set)
{
    this->matrices[matrices_set].resize(sizes_map->size());
    for (unsigned int i = 0; i < sizes_map->size(); ++i)
    {
        const auto &size_pair = (*sizes_map)[i];
        unsigned int rows = size_pair[0];
        unsigned int cols = size_pair[1];
        this->matrices[matrices_set][i].resize(rows, cols);
    }
}

void ElementBuilder::write_element_matrices(std::ofstream &outFile, unsigned int matrices_set)
{
    std::vector<DenseMatrix<Number>> &elem_matrices_vector = this->matrices[matrices_set];

    unsigned int elem_id = this->elem_data.elem->id();
    // Write the element ID as an unsigned int
    outFile.write(reinterpret_cast<const char *>(&elem_id), sizeof(elem_id));
    // Write the number of containers as an unsigned int
    unsigned int numContainers = static_cast<unsigned int>(elem_matrices_vector.size());
    outFile.write(reinterpret_cast<const char *>(&numContainers), sizeof(numContainers));

    for (const auto &elem_matrix : elem_matrices_vector)
    {
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

void ElementBuilder::read_element_matrices(std::ifstream &inFile, unsigned int matrices_set)
{
    std::vector<DenseMatrix<Number>> &elem_matrices_vector = this->matrices[matrices_set];

    // Read the element ID as an unsigned int
    unsigned int elem_id;
    inFile.read(reinterpret_cast<char *>(&elem_id), sizeof(elem_id));
    // Read the number of containers as an unsigned int
    unsigned int numContainers;
    inFile.read(reinterpret_cast<char *>(&numContainers), sizeof(numContainers));

    // Resize the containers vector to match the number of containers
    elem_matrices_vector.resize(numContainers);

    for (auto &elem_matrix : elem_matrices_vector)
    {
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

void ElementBuilder::write_element_vectors(std::ofstream &outFile, unsigned int vectors_set)
{
    // Write the number of containers as an unsigned int
    std::vector<DenseVector<Number>> &elem_vectors_vector = this->vectors[vectors_set];
    unsigned int elem_id = this->elem_data.elem->id();
    // Write the element ID as an unsigned int
    outFile.write(reinterpret_cast<const char *>(&elem_id), sizeof(elem_id));
    unsigned int numContainers = static_cast<unsigned int>(elem_vectors_vector.size());
    outFile.write(reinterpret_cast<const char *>(&numContainers), sizeof(numContainers));

    for (const auto &elem_vector : elem_vectors_vector)
    {
        const auto &values = elem_vector.get_values();
        // Write the size of the values vector as an unsigned int
        unsigned int valuesSize = elem_vector.size();
        outFile.write(reinterpret_cast<const char *>(&valuesSize), sizeof(valuesSize));
        // Write the actual values
        outFile.write(reinterpret_cast<const char *>(values.data()), values.size() * sizeof(Number));
    }
}

void ElementBuilder::read_element_vectors(std::ifstream &inFile, unsigned int vectors_set)
{
    std::vector<DenseVector<Number>> &elem_vectors_vector = this->vectors[vectors_set];
    // Read the number of containers as an unsigned int
    // Read the element ID as an unsigned int
    unsigned int elem_id;
    inFile.read(reinterpret_cast<char *>(&elem_id), sizeof(elem_id));
    unsigned int numContainers;
    inFile.read(reinterpret_cast<char *>(&numContainers), sizeof(numContainers));

    // Resize the containers vector to match the number of containers
    elem_vectors_vector.resize(numContainers);

    for (auto &elem_vector : elem_vectors_vector)
    {
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

void ElementBuilder::build_write_element_matrices_rhs(bool recompute_matrix = true)
{
    // Declare a performance log.  Give it a descriptive
    // string to identify what part of the code we are
    // logging, since there may be many PerfLogs in an
    // application.
    // PerfLog perf_log("Matrix Assembly");

    auto equation_systems = this->get_parent_object();

    System &system = equation_systems->get_system<System>(0);
    // Get a constant reference to the mesh object.
    const MeshBase &mesh = system.get_mesh();

    const DofMap &dof_map = system.get_dof_map();
    // Get a reference to the LinearImplicitSystem we are solving
    const auto &comm = system.comm();
    // A reference to the DofMap object for this system.  The DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the DofMap
    // in future examples.

    // Allocate memory for 3D nodal elements matrixes with transient effects
    // I can have up to 4 variables per node/edge
    // References to momentum equation right hand sides

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.

    unsigned int etype_num_index = mesh.get_elem_integer_index("ansys_etype");
    unsigned int mat_num_index = mesh.get_elem_integer_index("ansys_material_type");
    unsigned int real_constants_num_index = mesh.get_elem_integer_index("ansys_real_constant");

    dof_id_type etype_num;
    dof_id_type material_num;
    dof_id_type real_constant_num;

    auto context = std::dynamic_pointer_cast<ElementBuilderContext>(this->get_context());

    auto session = context->get_session();

    // if (old_time_substep == nullptr) // Here add the check for nonlinear elements
    // {
    // This function must be called prior to a check for the need of recomputing the matrix:
    // e.g. new materials, new real constants, nonlinear elements.

    const auto rank_str = std::to_string(comm.rank());
    std::vector<std::string> filenames = {
        "mtrx_" + rank_str + ".temp",
        "der_mtrx_" + rank_str + ".temp",
        "rhs_" + rank_str + ".temp"};

    auto element_types = context->get_element_types();
    std::filesystem::path emat_folderpath = session->get_emat_folderpath();
    // Handle renaming of rhs file if it exists
    std::filesystem::path rhs_path = emat_folderpath / ("rhs_" + rank_str + ".temp");
    if (std::filesystem::exists(rhs_path))
    {
        std::filesystem::path new_rhs_path = emat_folderpath / ("rhs_" + rank_str + "_prev.temp");
        std::filesystem::rename(rhs_path, new_rhs_path);
    }

    // Open all files using the session method
    auto streams = session->open_file_streams(emat_folderpath, filenames, FileMode::Write);

    // Downcast to use as ofstream
    std::ofstream *mtrx_output_file = dynamic_cast<std::ofstream *>(streams[0].get());
    std::ofstream *der_mtrx_output_file = dynamic_cast<std::ofstream *>(streams[1].get());
    std::ofstream *rhs_output_file = dynamic_cast<std::ofstream *>(streams[2].get());

    std::vector<dof_id_type> dof_indices;

    NumericVector<Number> &base_rhs = system.get_vector("base_rhs");

    if (recompute_matrix)
        for (const auto &elem : mesh.active_local_element_ptr_range())
        {
            dof_map.dof_indices(elem, dof_indices);
            // Start logging the shape function initialization.
            // This is done through a simple function call with
            // the name of the event to log.

            etype_num = elem->get_extra_integer(etype_num_index);
            material_num = elem->get_extra_integer(mat_num_index);
            real_constant_num = elem->get_extra_integer(real_constants_num_index);

            ElementType *elem_type = element_types->get_element_type(etype_num);

            // Do this elsewhere, logically must no stay here
            // load_step->get_materials()->compute_material_properties(elem->id(), conductive);
            this->set_element_data(elem, elem_type, material_num, real_constant_num);

            this->build_element();

            base_rhs.add_vector(vectors[0][0], dof_indices);

            // this->update_rhs();
            // Here the function to write the matrices and vectors down to file
            this->write_element_matrices(*mtrx_output_file, 0);

            this->write_element_matrices(*der_mtrx_output_file, 1);

            // this->write_element_vectors(*rhs_output_file, 0);
        }

    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        // Start logging the shape function initialization.
        // This is done through a simple function call with
        // the name of the event to log.
        if (!recompute_matrix)
        {
            dof_map.dof_indices(elem, dof_indices);

            etype_num = elem->get_extra_integer(etype_num_index);
            material_num = elem->get_extra_integer(mat_num_index);
            real_constant_num = elem->get_extra_integer(real_constants_num_index);

            ElementType *elem_type = element_types->get_element_type(etype_num);
            this->set_element_data(elem, elem_type, material_num, real_constant_num);
        }

        base_rhs.get(dof_indices, vectors[0][0].get_values());

        this->build_element_rhs();

        this->write_element_vectors(*rhs_output_file, 0);
    }
    session->close_file_streams(streams);
}

bool ElementBuilder::check_compute_matrix(std::shared_ptr<FEMDataContext> old_context = nullptr)
{
    // bool compute_matrix = false;
    auto eq_sys = this->get_parent_object();

    bool update_matrix;

    TransientLinearImplicitSystem &system = eq_sys->get_system<TransientLinearImplicitSystem>("ApecMLinearSystem");

    PetscLinearSolver<Number> *solver = dynamic_cast<PetscLinearSolver<Number> *>(system.get_linear_solver());

    std::shared_ptr<ElementBuilderContext> current_ctx = std::dynamic_pointer_cast<ElementBuilderContext>(this->get_context());
    std::shared_ptr<ElementBuilderContext> old_ctx = std::dynamic_pointer_cast<ElementBuilderContext>(old_context);

    if (old_context == nullptr)
    {
        // Here add the check for nonlinear elements
        // This function must be called prior to a check for the need of recomputing the matrix:
        // e.g. new materials, new real constants, nonlinear elements.
        update_matrix = true;
        solver->reuse_preconditioner(false);
    }
    // else if (active_time_substep == nullptr && prev_load_step != nullptr)
    // {
    //     // Here add the check for nonlinear elements
    //     // This function must be called prior to a check for the need of recomputing the matrix:
    //     // e.g. new materials, new real constants, nonlinear elements.
    //     update_matrix = true;
    //     solver->reuse_preconditioner(false);
    // }
    // in a new time step we could have a different dt, different bcs, different time integration effects
    // Here add the flag for ferromagnetic
    else
    {
        auto current_material_object = current_ctx->get_materials();
        auto old_material_object = old_ctx->get_materials();
        auto current_real_constant_object = current_ctx->get_real_constants();
        auto old_real_constant_object = old_ctx->get_real_constants();
        if (
            current_material_object->get_data_map() != old_material_object->get_data_map() ||
            current_real_constant_object->get_data_map() != old_real_constant_object->get_data_map())
        {
            update_matrix = true;
            solver->reuse_preconditioner(false);
        }
        else
        {
            update_matrix = false;
            solver->reuse_preconditioner(true);
        }
    }
    return update_matrix;
}
// {0: "Ke", 1: "Ce", 2: "Me", 3: "AAe", 4: "JAe", 5: "JVe"}
void ElementBuilder::build_solid97()
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
    FEBase unsigned int n_dofs_var = element_type->get_fe()->n_shape_functions();
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
    for (unsigned int qp = 0; qp < nq; qp++)
        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            for (unsigned int k = 0; k < 3; k++)
                vectors[0](dof_i + k * n_dofs_var) += scalar_shape_functions->JxW[qp] * scalar_shape_functions->phi[dof_i][qp]; //* js_loads[elem->id()][k];
}

// Temporary placeholder for build_solid236
void ElementBuilder::build_solid236()
{
    libmesh_error_msg("build_solid236 is not implemented yet.");
}

// Temporary placeholder for build_plane13
void ElementBuilder::build_plane13()
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

void ElementBuilder::build_element()
{
    auto etype = static_cast<ElementTypeEnum>(elem_data.element_type->get_keyopts()[0]);

    switch (etype)
    {
    case ElementTypeEnum::SOLID97:
        this->build_solid97();
        break;
    case ElementTypeEnum::SOLID236:
        this->build_solid236();
        break;
    case ElementTypeEnum::PLANE13:
        this->build_plane13();
        break;
    default:
        libmesh_error_msg("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
    }
}

void ElementBuilder::build_element_rhs()
{
    auto etype = static_cast<ElementTypeEnum>(elem_data.element_type->get_keyopts()[0]);

    switch (etype)
    {
    case ElementTypeEnum::SOLID97:
        this->build_solid97_rhs();
        break;
    case ElementTypeEnum::SOLID236:
        this->build_solid236_rhs();
        break;
    case ElementTypeEnum::PLANE13:
        this->build_plane13_rhs();
        break;
    default:
        libmesh_error_msg("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
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

    unsigned int n_dofs_var = elem_data.element_type->fe_array->n_shape_functions();
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
