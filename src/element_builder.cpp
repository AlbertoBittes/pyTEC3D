#include "element_builder.h"
// #include <libmesh/equation_systems.h>
#include <libmesh/dof_map.h>
#include <libmesh/numeric_vector.h>
// #include <libmesh/string_to_enum.h>
// Define the Finite Element object.
#include <libmesh/fe.h>
#include <libmesh/petsc_linear_solver.h>

// Define Gauss quadrature rules.
#include "session.h"       // Ensure the Session class is included
#include "element_types.h" // Include the header defining ElementType
#include "loads_object.h"
#include "apec_m_system_objects.h"

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

    unsigned int num_matrices = static_cast<unsigned int>(sizes_map->size());
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

void ElementBuilder::build_write_element_matrices_rhs(bool recompute_matrix)
{
    // Declare a performance log.  Give it a descriptive
    // string to identify what part of the code we are
    // logging, since there may be many PerfLogs in an
    // application.
    // PerfLog perf_log("Matrix Assembly");

    auto equation_systems = this->get_parent_object();
    // Get a constant reference to the mesh object.
    const MeshBase &mesh = equation_systems->get_mesh();
    // Get a reference to the LinearImplicitSystem we are solving
    const auto &comm = equation_systems->comm();
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

    auto session = this->context<Session>();

    // if (old_time_substep == nullptr) // Here add the check for nonlinear elements
    // {
    // This function must be called prior to a check for the need of recomputing the matrix:
    // e.g. new materials, new real constants, nonlinear elements.

    const auto rank_str = std::to_string(comm.rank());

    std::filesystem::path emat_folderpath = session->get_path("emat_folderpath");
    // Handle renaming of rhs file if it exists
    std::filesystem::path rhs_path = emat_folderpath / ("rhs_" + rank_str + ".temp");
    if (std::filesystem::exists(rhs_path))
    {
        std::filesystem::path new_rhs_path = emat_folderpath / ("rhs_" + rank_str + "_prev.temp");
        std::filesystem::rename(rhs_path, new_rhs_path);
    }

    std::vector<dof_id_type> dof_indices;

    if (recompute_matrix)
    {
        std::vector<std::string> filenames = {
            "mtrx_" + rank_str + ".temp",
            "nodal_sol_mtrx_" + rank_str + ".temp",
            "elem_sol_mtrx_" + rank_str + ".temp",
            "base_rhs_" + rank_str + ".temp"};

        // Open all files using the session method
        auto streams = session->open_file_streams(emat_folderpath, filenames, FileMode::Write);

        // Downcast to use as ofstream
        std::ofstream *mtrx_output_file = dynamic_cast<std::ofstream *>(streams[0].get());
        std::ofstream *nodal_sol_mtrx_output_file = dynamic_cast<std::ofstream *>(streams[1].get());
        std::ofstream *elem_sol_mtrx_output_file = dynamic_cast<std::ofstream *>(streams[2].get());
        std::ofstream *base_rhs_output_file = dynamic_cast<std::ofstream *>(streams[3].get());

        for (const auto &elem : mesh.active_local_element_ptr_range())
        {
            main_dof_map->dof_indices(elem, dof_indices);
            // Start logging the shape function initialization.
            // This is done through a simple function call with
            // the name of the event to log.
            // Do this elsewhere, logically must no stay here
            // load_step->get_materials()->compute_material_properties(elem->id(), conductive);
            this->set_element_data(elem);

            this->build_element();

            // this->update_rhs();
            // Here the function to write the matrices and vectors down to file
            this->write_element_matrices(*mtrx_output_file, 0);

            this->write_element_matrices(*nodal_sol_mtrx_output_file, 1);

            this->write_element_matrices(*elem_sol_mtrx_output_file, 2);

            this->write_element_vectors(*base_rhs_output_file, 0);
        }

        session->close_file_streams(streams);
    }
    std::vector<std::string> inp_filenames = {
        "base_rhs_" + rank_str + ".temp"};
    std::vector<std::string> out_filenames = {
        "rhs_" + rank_str + ".temp",
    };

    auto streams_in = session->open_file_streams(emat_folderpath, inp_filenames, FileMode::Read);
    auto streams_out = session->open_file_streams(emat_folderpath, out_filenames, FileMode::Write);
    // Downcast to use as ifstream
    std::ifstream *base_rhs_input_file = dynamic_cast<std::ifstream *>(streams_in[0].get());
    std::ofstream *rhs_output_file = dynamic_cast<std::ofstream *>(streams_out[0].get());

    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        // Start logging the shape function initialization.
        // This is done through a simple function call with
        // the name of the event to log.
        main_dof_map->dof_indices(elem, dof_indices);

        this->set_element_data(elem);

        this->read_element_vectors(*base_rhs_input_file, 0);

        this->build_element_rhs();

        this->write_element_vectors(*rhs_output_file, 0);
    }
    session->close_file_streams(streams_in);
    session->close_file_streams(streams_out);
}

bool ElementBuilder::check_compute_matrix(std::shared_ptr<TypedContextView<Session, ElementTypeFactory, MaterialsObject, RealConstantsObject, LoadsObject>> old_view)
{
    auto eq_sys = this->get_parent_object();
    bool update_matrix;

    ApecMLinearSystem &system = eq_sys->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    PetscLinearSolver<Number> *solver = dynamic_cast<PetscLinearSolver<Number> *>(system.get_linear_solver());

    // Get current and old materials/real constants from the views
    auto current_material_object = this->context<MaterialsObject>();
    auto old_material_object = old_view ? old_view->template get<MaterialsObject>() : nullptr;
    auto current_real_constant_object = this->context<RealConstantsObject>();
    auto old_real_constant_object = old_view ? old_view->template get<RealConstantsObject>() : nullptr;

    if (old_view == nullptr)
    {
        update_matrix = true;
        solver->reuse_preconditioner(false);
    }
    else if (!old_material_object || !old_real_constant_object)
    {
        update_matrix = true;
        solver->reuse_preconditioner(false);
    }
    else if (
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
    FEBase *fe = dynamic_cast<FEBase *>(elem_data.element_type->fe_array[0].get());
    FEVectorBase *fe_vec = dynamic_cast<FEVectorBase *>(elem_data.element_type->fe_array[1].get());
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real> &JxW = fe->get_JxW();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real>> &phi = fe->get_phi();

    // The element shape function gradients for the velocity
    // variables evaluated at the quadrature points.
    const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

    // The velocity shape function gradients at interior
    // quadrature points.
    const std::vector<std::vector<RealVectorValue>> &phi_vec = fe_vec->get_phi();

    const std::vector<std::vector<RealGradient>> &curl_phi_vec = fe_vec->get_curl_phi();

    fe->reinit(elem_data.elem);
    fe_vec->reinit(elem_data.elem);

    this->prepare_element_matrix_building();

    // Get the scalar shape functions and vector shape functions
    unsigned int n_dofs_var = fe->n_shape_functions();
    unsigned int n_dofs_A = n_dofs_var * 3;
    unsigned int nq = fe->n_quadrature_points();

    auto materials = this->context<MaterialsObject>();

    materials->compute_material_properties(this->elem_data.material_num, this->elem_data.element_type->is_conductive());
    auto &mu_matrix = materials->get_mu_matrix();
    auto &rs_matrix = materials->get_rs_matrix();

    auto &Ke = this->matrices[0][0];
    auto &Ce = this->matrices[0][1];
    auto &Me = this->matrices[1][0];
    auto &AAe = this->matrices[1][1];
    auto &JAe = this->matrices[2][0];
    auto &JVe = this->matrices[2][1];
    auto &Fe = this->vectors[0][0];

    if (this->elem_data.element_type->is_conductive())
    {
        Real CAV_value;

        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            {
                JAe(qp, dof_i) += phi[dof_i][qp];

                for (unsigned int k = 0; k < 3; k++)
                    JVe(qp + k * n_dofs_var, dof_i) += dphi[dof_i][qp](k);

                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    for (unsigned int l = 0; l < 3; l++)
                    {
                        Ce(dof_i + l * n_dofs_var, dof_j + l * n_dofs_var) += JxW[qp] / rs_matrix[l] * (phi[dof_i][qp] * phi[dof_j][qp]);
                        CAV_value = JxW[qp] / rs_matrix[l] * phi[dof_i][qp] * dphi[dof_j][qp](l);
                        Ce(dof_i + l * n_dofs_var, n_dofs_A + dof_j) += CAV_value;
                        Ce(n_dofs_A + dof_j, dof_i + l * n_dofs_var) += CAV_value;
                        Ce(n_dofs_A + dof_i, n_dofs_A + dof_j) += JxW[qp] * dphi[dof_i][qp](l) / rs_matrix[l] * dphi[dof_j][qp](l);
                        for (unsigned int k = 0; k < 3; k++)
                        {
                            Ke(dof_i + l * n_dofs_var, dof_j + l * n_dofs_var) += JxW[qp] / mu_matrix[k] * (dphi[dof_i][qp](k) * dphi[dof_j][qp](k));
                            AAe(dof_i + k * n_dofs_var, dof_j + l * n_dofs_var) += JxW[qp] * (phi_vec[dof_i * 3 + k][qp] * curl_phi_vec[dof_j * 3 + l][qp]);
                            Me(dof_i + k * n_dofs_var, dof_j + l * n_dofs_var) += JxW[qp] * (phi_vec[dof_i * 3 + k][qp] * phi_vec[dof_j * 3 + l][qp]);
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
                            Ke(dof_i + l * n_dofs_var, dof_j + l * n_dofs_var) += JxW[qp] / mu_matrix[k] * (dphi[dof_i][qp](k) * dphi[dof_j][qp](k));
                            AAe(dof_i + k * n_dofs_var, dof_j + l * n_dofs_var) += JxW[qp] * (phi_vec[dof_i * 3 + k][qp] * curl_phi_vec[dof_j * 3 + l][qp]);
                            Me(dof_i + k * n_dofs_var, dof_j + l * n_dofs_var) += JxW[qp] * (phi_vec[dof_i * 3 + k][qp] * phi_vec[dof_j * 3 + l][qp]);
                        }
    for (unsigned int qp = 0; qp < nq; qp++)
        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            for (unsigned int k = 0; k < 3; k++)
                Fe(dof_i + k * n_dofs_var) += JxW[qp] * phi[dof_i][qp]; //* js_loads[elem->id()][k];
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
    FEBase *fe = dynamic_cast<FEBase *>(elem_data.element_type->fe_array[0].get());
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real> &JxW = fe->get_JxW();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real>> &phi = fe->get_phi();

    // The element shape function gradients for the velocity
    // variables evaluated at the quadrature points.
    const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

    fe->reinit(elem_data.elem);

    this->prepare_element_matrix_building();

    // Get the scalar shape functions and vector shape functions
    unsigned int n_dofs_var = fe->n_shape_functions();
    unsigned int nq = fe->n_quadrature_points();

    auto materials = this->context<MaterialsObject>();

    materials->compute_material_properties(this->elem_data.material_num, this->elem_data.element_type->is_conductive());
    auto &mu_matrix = materials->get_mu_matrix();
    auto &rs_matrix = materials->get_rs_matrix();

    auto &Ke = this->matrices[0][0];
    auto &Ce = this->matrices[0][1];
    auto &Me = this->matrices[1][0];
    auto &AAe = this->matrices[1][1];
    auto &JAe = this->matrices[2][0];
    auto &Fe = this->vectors[0][0];

    if (elem_data.element_type->is_conductive())
    {
        // Now start logging the element matrix computation
        for (unsigned int qp = 0; qp < nq; qp++)
        {
            for (unsigned int dof_i = 0; dof_i < 2 * n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < 2 * n_dofs_var; dof_j++)
                {
                    if (dof_j < n_dofs_var)
                        AAe(dof_i, dof_j) += JxW[qp] * std::pow(-1, dof_i) * (phi[dof_i / 2][qp] * dphi[dof_j][qp](1 - (dof_i % 2)));

                    if ((dof_i + dof_j) % 2 == 0)
                        Me(dof_i, dof_j) += JxW[qp] * (phi[dof_i / 2][qp] * phi[dof_j / 2][qp]);
                }
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            {
                // matrices[4](0, i) -= scalar_shape_functions->JxW[qp] * scalar_shape_functions->phi[i][qp] / scalar_shape_functions->phi[0].size();
                JAe(qp, dof_i) -= phi[dof_i][qp];

                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                {
                    Ke(dof_i, dof_j) += JxW[qp] * (dphi[dof_i][qp] * dphi[dof_j][qp]) / mu_matrix[0];
                    Ce(dof_i, dof_j) += JxW[qp] / rs_matrix[0] * (phi[dof_i][qp] * phi[dof_j][qp]);
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
                        AAe(dof_i, dof_j) += JxW[qp] * std::pow(-1, dof_i) * (phi[dof_i / 2][qp] * dphi[dof_j][qp](1 - (dof_i % 2)));

                    if ((dof_i + dof_j) % 2 == 0)
                        Me(dof_i, dof_j) += JxW[qp] * (phi[dof_i / 2][qp] * phi[dof_j / 2][qp]);
                }

            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    for (unsigned int l = 0; l < 3; l++)
                        for (unsigned int k = 0; k < 3; k++)
                            Ke(dof_i, dof_j) += JxW[qp] * (dphi[dof_i][qp] * dphi[dof_j][qp]) / mu_matrix[0];
        }

    for (unsigned int qp = 0; qp < nq; qp++)
        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            Fe(dof_i) += JxW[qp] * phi[dof_i][qp]; // * js_loads[elem->id()][2];
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

    unsigned int n_dofs_var = elem_data.elem->n_nodes();

    auto loads_object = this->context<LoadsObject>();

    auto loads_map = loads_object->get_data_map();

    unsigned int system_num = loads_map->number();

    Real js[3];

    unsigned int jsx_var_num = loads_map->variable_number("Js_x");
    unsigned int jsy_var_num = loads_map->variable_number("Js_y");
    unsigned int jsz_var_num = loads_map->variable_number("Js_z");

    dof_id_type dof_jx = elem_data.elem->dof_number(system_num, jsx_var_num, 0);
    dof_id_type dof_jy = elem_data.elem->dof_number(system_num, jsy_var_num, 0);
    dof_id_type dof_jz = elem_data.elem->dof_number(system_num, jsz_var_num, 0);

    js[0] = loads_map->current_local_solution->el(dof_jx);
    js[1] = loads_map->current_local_solution->el(dof_jy);
    js[2] = loads_map->current_local_solution->el(dof_jz);

    auto &Fe = this->vectors[0][0];

    for (unsigned int k = 0; k < 3; k++)
        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            Fe(dof_i + k * n_dofs_var) *= js[k];
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

    unsigned int n_dofs_var = elem_data.elem->n_nodes();

    auto loads_object = this->context<LoadsObject>();

    auto loads_map = loads_object->get_data_map();

    Real js;

    unsigned int jsz_var_num = loads_map->variable_number("Js_z");

    dof_id_type dof_jz = elem_data.elem->dof_number(0, jsz_var_num, 0);

    js = loads_map->solution->el(dof_jz);

    auto &Fe = this->vectors[0][0];

    for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
        Fe(dof_i) *= js;
}

void ElementBuilder::build_nodal_solution_vector()
{
    // Build the nodal solution vector
    // This is a placeholder for the actual implementation
    auto etype = static_cast<ElementTypeEnum>(elem_data.element_type->get_keyopts()[0]);

    switch (etype)
    {
    case ElementTypeEnum::SOLID97:
        this->build_solid97_nodal_solution_vector();
        break;
    case ElementTypeEnum::SOLID236:
        this->build_solid236_nodal_solution_vector();
        break;
    case ElementTypeEnum::PLANE13:
        this->build_plane13_nodal_solution_vector();
        break;
    default:
        libmesh_error_msg("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
    }
}

void ElementBuilder::build_solid97_nodal_solution_vector()
{
    // ApecMSolutionSystem *solution_system = dynamic_cast<ApecMSolutionSystem *>(this->get_parent_object());
    fill_solution_vector_support(main_sys->current_local_solution.get());

    // Get Ae
    DenseVector<Number> Ae;
    for (unsigned int k = 0; k < 3; k++)
    {
        Ae.append(solution_vector_support.sol_vectors[k]);
    }
    auto &AAe = this->matrices[1][0];
    DenseVector<Number> Be_q(AAe.m());
    AAe.vector_mult(Be_q, Ae);
    auto &Me = this->matrices[1][1];
    auto &Be = this->vectors[2][0];
    Me.cholesky_solve(Be_q, Be);
    // fill_solution_vector_support(main_sys->old_local_solution);

    // // Get Ae
    // DenseVector<Number> Ae_old;
    // for (unsigned int k = 0; k < 3; k++)
    // {
    //     Ae_old.append(solution_vector_support.sol_vectors[k]);
    // }
    // Ae_dt.add(-1.0, Ae_old);
    // Ae_dt.scale(1.0 / time_step->get_active_substep()->get_dt());
}

void ElementBuilder::build_element_solution_vector()
{
    // Build the element solution vector
    // This is a placeholder for the actual implementation
    auto etype = static_cast<ElementTypeEnum>(elem_data.element_type->get_keyopts()[0]);

    switch (etype)
    {
    case ElementTypeEnum::SOLID97:
        this->build_solid97_element_solution_vector();
        break;
    case ElementTypeEnum::SOLID236:
        this->build_solid236_element_solution_vector();
        break;
    case ElementTypeEnum::PLANE13:
        this->build_plane13_element_solution_vector();
        break;
    default:
        libmesh_error_msg("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
    }
}

void ElementBuilder::build_solution_vector(SystemType sys_type)
{
    switch (sys_type)
    {
    case SystemType::NODAL_EXPLICIT:
        build_nodal_solution_vector();
        break;
    case SystemType::ELEMENT_EXPLICIT:
        build_element_solution_vector();
        break;
    default:
        libmesh_error_msg("System type not supported: " + std::to_string(static_cast<unsigned int>(sys_type)));
    }
}

void ElementBuilder::build_write_solution_vectors(SystemType sys_type)
{
    ApecMSolutionSystem *solution_system = dynamic_cast<ApecMSolutionSystem *>(this->get_parent_object());

    ApecMEquationSystems *eq_sys = dynamic_cast<ApecMEquationSystems *>(&solution_system->get_equation_systems());

    ApecMLinearSystem &system = eq_sys->get_system<ApecMLinearSystem>("ApecMLinearSystem");

    const MeshBase &mesh = system.get_mesh();

    // Get a reference to the LinearImplicitSystem we are solving
    const auto &comm = system.comm();

    // A reference to the DofMap object for this system.  The DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the DofMap
    // in future examples.
    const DofMap &dof_map = system.get_dof_map();
    const DofMap &solution_dof_map = solution_system->get_dof_map();

    auto session = this->context<Session>();

    auto time_substep = this->context<TimeStep>()->get_active_substep();

    auto sys_name = system.name();

    const auto rank_str = std::to_string(comm.rank());

    std::filesystem::path emat_folderpath = session->get_path("emat_folderpath");

    std::vector<std::string> inp_filenames = {};
    std::vector<std::string> out_filenames = {};

    unsigned int matrices_set = 0;
    unsigned int vec_set = 0;
    if (sys_name == "NodalSolutionSystem")
    {
        inp_filenames.push_back("nodal_sol_mtrx_" + rank_str + ".temp");
        out_filenames.push_back("nodal_sol_vec_" + rank_str + ".temp");
        matrices_set = 1;
        vec_set = 2;
    }
    else if (sys_name == "ElementSolutionSystem")
    {
        inp_filenames.push_back("elem_sol_mtrx_" + rank_str + ".temp");
        out_filenames.push_back("elem_sol_vec_" + rank_str + ".temp");
        matrices_set = 2;
        vec_set = 3;
    }
    auto inp_streams = session->open_file_streams(emat_folderpath, inp_filenames, FileMode::Read);

    std::ifstream *solution_mtrx_input_file = dynamic_cast<std::ifstream *>(inp_streams[0].get());

    auto out_streams = session->open_file_streams(emat_folderpath, out_filenames, FileMode::Write);

    std::ofstream *solution_vec_output_file = dynamic_cast<std::ofstream *>(out_streams[0].get());

    auto element_types = this->context<ElementTypeFactory>();

    // std::vector<unsigned int> A_num_var;
    // std::vector<unsigned int> B_num_var;
    // std::vector<unsigned int> J_num_var;
    // unsigned int v_var;

    // A_num_var.push_back(system.variable_number("A_z"));
    // B_num_var.push_back(derived_quantities_system.variable_number("B_x"));
    // B_num_var.push_back(derived_quantities_system.variable_number("B_y"));
    // J_num_var.push_back(derived_quantities_system.variable_number("J_elem_z"));

    // if (dim == 3)
    // {
    //     A_num_var.push_back(system.variable_number("A_x"));
    //     A_num_var.push_back(system.variable_number("A_y"));
    //     v_var = system.variable_number("V");

    //     B_num_var.push_back(derived_quantities_system.variable_number("B_z"));
    //     J_num_var.push_back(derived_quantities_system.variable_number("J_elem_x"));
    //     J_num_var.push_back(derived_quantities_system.variable_number("J_elem_y"));
    // }

    // auto assembler = system.get_assembler();
    // auto time_substep = this->context<TimeStep>()->get_active_substep();
    // auto transient_effect = time_substep->get_transient_effect();
    // Storage for the stress dof indices on each element
    // std::vector<std::vector<dof_id_type>>
    //     A_dof_indices_var;
    // std::vector<std::vector<dof_id_type>> B_dof_indices_var;
    // std::vector<std::vector<dof_id_type>> J_dof_indices_var;
    // std::vector<dof_id_type> A_scalar_dof_indices_var;
    // std::vector<dof_id_type> B_scalar_dof_indices_var;
    // std::vector<dof_id_type> J_scalar_dof_indices_var;
    // std::vector<dof_id_type> V_indices_var;

    // std::shared_ptr<LoadsObject> loads_object = eq_sys.parameters.get<std::shared_ptr<LoadsObject>>("loads_object");

    // DenseVector<Number> Ae;
    // DenseVector<Number> Ae_dt;
    // DenseVector<Number> Ve;
    // DenseVector<Number> Ve_dt;
    // DenseVector<Number> Be;
    // // DenseVector<Number> Fe;
    // DenseVector<Number> Be_q;
    // DenseVector<Number> J_elem;
    // DenseVector<Number> J_elem_A;
    // DenseVector<Number> J_elem_V;
    // // DenseVector<Number> B_elem;
    // // DenseVector<Number> F_elem;
    // char n_J_dofs;
    // if (dim == 2)
    // {
    //     n_J_dofs = 1;
    //     J_elem.resize(n_J_dofs);
    //     // B_elem.resize(2);
    //     // F_elem.resize(2);
    // }
    // else
    // {
    //     n_J_dofs = 3;
    //     J_elem.resize(n_J_dofs);
    //     B_elem.resize(n_J_dofs);
    //     F_elem.resize(n_J_dofs);
    // }
    // char n_A_dofs;
    // char n_V_dofs;
    // char n_B_dofs;
    // bool conductive;
    // std::set<subdomain_id_type> *conductive_subdomain = eq_sys.parameters.get<std::set<subdomain_id_type> *>("conductive_subdomain");
    std::vector<dof_id_type> dof_indices;

    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        // Start logging the shape function initialization.
        // This is done through a simple function call with
        // the name of the event to log.

        this->set_element_data(elem);

        this->read_element_matrices(*solution_mtrx_input_file, matrices_set);

        this->build_solution_vector(sys_type);

        this->write_element_vectors(*solution_vec_output_file, vec_set);
        // loads_object->local_elem_id = local_elem_idx;

        // if (elem->subdomain_id() == 1)
        //     conductive = true;
        // else
        //     conductive = false;

        // if (!trans_effects)
        //     conductive = false;

        // loads_object->compute_material_properties(conductive);

        // for (unsigned int var = 0; var < A_num_var.size(); var++)
        // {
        //     dof_map.dof_indices(elem, A_scalar_dof_indices_var, A_num_var[var]);
        //     A_dof_indices_var.push_back(A_scalar_dof_indices_var);
        // }
        // for (unsigned int var = 0; var < B_num_var.size(); var++)
        // {
        //     dof_map.dof_indices(elem, A_scalar_dof_indices_var, A_num_var[var]);
        //     A_dof_indices_var.push_back(A_scalar_dof_indices_var);
        // }
        // for (unsigned int var = 0; var < A_num_var.size(); var++)
        // {
        //     dof_map.dof_indices(elem, A_scalar_dof_indices_var, A_num_var[var]);
        //     A_dof_indices_var.push_back(A_scalar_dof_indices_var);
        // }
        // dof_map.dof_indices(elem, A_indices_var, a_var);
        // n_A_dofs = std::size(A_indices_var);
        // Ae.resize(n_A_dofs);
        // Ae_dt.resize(n_A_dofs);

        // derived_quantities_system_dof_map.dof_indices(elem, B_indices_var, b_n_var);
        // n_B_dofs = std::size(B_indices_var);
        // Be_q.resize(n_B_dofs);
        // Be.resize(n_B_dofs);

        // derived_quantities_system_elem_dof_map.dof_indices(elem, J_indices_var, j_elem_var);
        // B_elem.zero();
        // F_elem.zero();
        // if (conductive)
        // {
        //     J_elem.zero();
        //     if (dim == 3)
        //     {
        //         dof_map.dof_indices(elem, V_indices_var, v_var);
        //         n_V_dofs = std::size(V_indices_var);
        //         Ve.resize(n_V_dofs);
        //         Ve_dt.resize(n_V_dofs);
        //         J_elem_A.zero();
        //         J_elem_V.zero();
        //     }
        // }

        // for (std::size_t j = 0; j < n_A_dofs; j++)
        // {
        //     Ae(j) = system.current_solution(A_indices_var[j]);
        //     if (conductive)
        //     {
        //         Ae_dt(j) = (Ae(j) - system.old_solution(A_indices_var[j])) / dt;

        //         if (dim == 3)
        //         {
        //             if (j < n_V_dofs)
        //             {
        //                 Ve_dt(j) = (system.current_solution(V_indices_var[j]) - system.old_solution(V_indices_var[j])) / dt;
        //             }
        //         }
        //     }
        // }

        // loads_object->aae_matrix_derived_map[local_elem_idx].vector_mult(Be_q, Ae);

        // // Below we project each component of the stress tensor onto a L2_LAGRANGE discretization.
        // // Note that this gives a discontinuous stress plot on element boundaries, which is
        // // appropriate. We then also get the von Mises stress from the projected stress tensor.
        // loads_object->me_matrix_derived_map[local_elem_idx].cholesky_solve(Be_q, Be);

        // unsigned int j_qp_size;
        // if (conductive)
        // {
        //     if (dim == 2)
        //     {
        //         j_qp_size = loads_object->jae_matrix_derived_map[local_elem_idx].m();
        //         J_elem_A.resize(j_qp_size);
        //         loads_object->jae_matrix_derived_map[local_elem_idx].vector_mult(J_elem_A, Ae);
        //         for (char k = 0; k < n_J_dofs; k++)
        //             for (char qp = 0; qp < j_qp_size; qp++)
        //                 J_elem(k) -= J_elem_A(qp + k) / loads_object->rs_matrix[k] / (j_qp_size / n_J_dofs);
        //     }
        //     else
        //     {
        //         j_qp_size = loads_object->jae_matrix_derived_map[local_elem_idx].m();
        //         J_elem_A.resize(j_qp_size);
        //         J_elem_V.resize(j_qp_size);
        //         loads_object->jae_matrix_derived_map[local_elem_idx].vector_mult(J_elem_A, Ae_dt);
        //         loads_object->jve_matrix_derived_map[local_elem_idx].vector_mult(J_elem_V, Ve_dt);
        //         for (unsigned int k = 0; k < n_J_dofs; k++)
        //             for (unsigned int qp = 0; qp < j_qp_size / n_J_dofs; qp++)
        //                 J_elem(k) += (-J_elem_A(qp * 3 + k) - J_elem_V(qp * 3 + k)) / loads_object->rs_matrix[k] / (j_qp_size / n_J_dofs);
        //     }
        // }
        // else
        // {
        //     for (unsigned int k = 0; k < n_J_dofs; k++)
        //     {
        //         J_elem(k) = loads_object->compute_js(k);
        //     }
        // }
    }
    // Close file streams
    solution_mtrx_input_file->close();
    solution_vec_output_file->close();
}