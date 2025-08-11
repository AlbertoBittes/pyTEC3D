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

ElementBuilder::ElementBuilder(EquationSystems *equation_systems,
                               std::shared_ptr<ContextView> ctx_view,
                               std::shared_ptr<ContextView> old_ctx_view)
    : FEMServiceWithOldContext<Session, ElementTypesFactory, TimeStep, MaterialsObject, RealConstantsObject, LoadsObject>(ctx_view, old_ctx_view),
      FEParentObjectBase<EquationSystems>(equation_systems),
      matrices(3),
      vectors(5)
{
}

void ElementBuilder::init()
{
    MeshBase &mesh = this->get_parent_object()->get_mesh();
    mesh_data_index.etype_num_index = mesh.get_elem_integer_index("ansys_etype");
    mesh_data_index.material_num_index = mesh.get_elem_integer_index("ansys_material_type");
    mesh_data_index.real_constants_num_index = mesh.get_elem_integer_index("ansys_real_constant");
    main_sys = &this->get_parent_object()->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    main_dof_map = &main_sys->get_dof_map();
}

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
    this->matrices[matrices_set].resize(num_matrices);
    for (unsigned int i = 0; i < num_matrices; ++i)
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

void ElementBuilder::fill_solution_vector_support(const NumericVector<Number> *solution, unsigned int vectors_set)
{
    solution_vector_support.clear(vectors_set);

    auto active_vars = elem_data.element_type->get_active_variables(Globals::VariableCategory::DOFVariables);
    for (const auto &var_enum : active_vars)
    {
        unsigned int system_var_num = elem_data.element_type->get_variable_information()[var_enum].system_var_num;
        // solution_vector_support.var_nums.push_back(system_var_num);

        std::vector<dof_id_type> dof_indices;
        main_dof_map->dof_indices(elem_data.elem, dof_indices, system_var_num);
        // solution_vector_support.dof_indices_var.push_back(dof_indices);

        DenseVector<Number> local_sol(dof_indices.size());
        for (unsigned int i = 0; i < dof_indices.size(); ++i)
        {
            local_sol(i) = solution->el(dof_indices[i]);
        }
        solution_vector_support.sol_vectors[vectors_set].push_back(local_sol);
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
        libMesh::out << "Recomputing element matrices and vectors for current time step" << std::endl;
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

    libMesh::out << "Applying loads to RHS" << std::endl;
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

bool ElementBuilder::check_compute_matrix()
{
    auto eq_sys = this->get_parent_object();
    bool update_matrix;

    ApecMLinearSystem &system = eq_sys->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    PetscLinearSolver<Number> *solver = dynamic_cast<PetscLinearSolver<Number> *>(system.get_linear_solver());

    auto old_view = this->get_old_context_view();
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
        *(current_material_object->get_mat_like_map()) != *(old_material_object->get_mat_like_map()) ||
        *(current_real_constant_object->get_mat_like_map()) != *(old_real_constant_object->get_mat_like_map()))
    {
        update_matrix = true;
        solver->reuse_preconditioner(false);
    }
    else
    {
        update_matrix = false;
        solver->reuse_preconditioner(true);
    }
    if (update_matrix)
        libMesh::out << "Element matrices will be recomputed for current time step" << std::endl;
    else
        libMesh::out << "Element matrices will be reused from previous time step" << std::endl;
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

    elem_builder_misc_data.materials_object->fill_mu_matrix(this->elem_data.material_num);
    auto &mu_matrix = elem_builder_misc_data.materials_object->get_mu_matrix();

    auto &Ke = this->matrices[0][0];
    auto &Ce = this->matrices[0][1];
    auto &Me = this->matrices[1][0];
    auto &AAe = this->matrices[1][1];
    auto &JAe = this->matrices[2][0];
    auto &JVe = this->matrices[2][1];
    auto &Fe = this->vectors[0][0];

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

    switch (elem_data.element_type->get_physical_region())
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int k = 0; k < 3; k++)
                    Fe(dof_i + k * n_dofs_var) += JxW[qp] * phi[dof_i][qp]; //* js_loads[elem->id()][k];
        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    {
        Real CAV_value;
        elem_builder_misc_data.materials_object->fill_rs_matrix(this->elem_data.material_num);
        auto &rs_matrix = elem_builder_misc_data.materials_object->get_rs_matrix();
        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            {

                for (unsigned int k = 0; k < 3; k++)
                {
                    JVe(qp + k * n_dofs_var, dof_i) += dphi[dof_i][qp](k);
                    JAe(qp + k * n_dofs_var, dof_i + k * n_dofs_var) += phi[dof_i][qp];
                }

                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    for (unsigned int l = 0; l < 3; l++)
                    {
                        Ce(dof_i + l * n_dofs_var, dof_j + l * n_dofs_var) += JxW[qp] / rs_matrix[l] * (phi[dof_i][qp] * phi[dof_j][qp]);
                        CAV_value = JxW[qp] / rs_matrix[l] * phi[dof_i][qp] * dphi[dof_j][qp](l);
                        Ce(dof_i + l * n_dofs_var, n_dofs_A + dof_j) += CAV_value;
                        Ce(n_dofs_A + dof_j, dof_i + l * n_dofs_var) += CAV_value;
                        Ce(n_dofs_A + dof_i, n_dofs_A + dof_j) += JxW[qp] * dphi[dof_i][qp](l) / rs_matrix[l] * dphi[dof_j][qp](l);
                    }
            }

        // Here when you will implement loads with F
        // for (unsigned int qp = 0; qp < nq; qp++)
        //     for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
        //         for (unsigned int k = 0; k < 3; k++)
        //             Fe(dof_i + k * n_dofs_var) += JxW[qp] * phi[dof_i][qp]; //* js_loads[elem->id()][k];
        break;
    }
    default:
        libmesh_error_msg("Physical region not implemented");
    }
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

    const std::vector<Point> &xyz = fe->get_xyz();

    fe->reinit(elem_data.elem);

    this->prepare_element_matrix_building();

    // Get the scalar shape functions and vector shape functions
    //   const unsigned int n_sf =
    // FEInterface::n_dofs(cfe->get_fe_type(), elem);
    unsigned int n_dofs_var = fe->n_shape_functions();
    unsigned int nq = fe->n_quadrature_points();

    elem_builder_misc_data.materials_object->fill_mu_matrix(this->elem_data.material_num);
    auto &mu_matrix = elem_builder_misc_data.materials_object->get_mu_matrix();

    auto &Ke = this->matrices[0][0];
    auto &Ce = this->matrices[0][1];
    auto &Me = this->matrices[1][0];
    auto &AAe = this->matrices[1][1];
    auto &JAe = this->matrices[2][0];
    auto &Fe = this->vectors[0][0];

    unsigned int axisymm = elem_data.element_type->keyopts[3];

    for (unsigned int qp = 0; qp < nq; qp++)
    {
        Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
            {
                AAe(dof_i, dof_j) += JxW[qp] * (phi[dof_i][qp] * dphi[dof_j][qp](1));
                AAe(dof_i + n_dofs_var, dof_j) -= JxW[qp] * (phi[dof_i][qp] * dphi[dof_j][qp](0));
                Me(dof_i, dof_j) += JxW[qp] * (phi[dof_i][qp] * phi[dof_j][qp]);
                Me(dof_i + n_dofs_var, dof_j + n_dofs_var) += JxW[qp] * (phi[dof_i][qp] * phi[dof_j][qp]);
            }

        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                Ke(dof_i, dof_j) += JxW[qp] * (dphi[dof_i][qp] * dphi[dof_j][qp] + axisymm * phi[dof_i][qp] * phi[dof_j][qp] / pow(r, 2)) * pow(r, axisymm);
    }
    Ke.scale(1.0 / mu_matrix[0]);
    PhysicalRegionEnum physical_region = elem_data.element_type->get_physical_region();

    switch (physical_region)
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        for (unsigned int qp = 0; qp < nq; qp++)
        {
            Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                Fe(dof_i) += JxW[qp] * phi[dof_i][qp] * pow(r, axisymm); // * js_loads[elem->id()][2];
        }
        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    {
        elem_builder_misc_data.materials_object->fill_rs_matrix(this->elem_data.material_num);
        auto &rs_matrix = elem_builder_misc_data.materials_object->get_rs_matrix();
        for (unsigned int qp = 0; qp < nq; qp++)
        {
            Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            {
                // matrices[4](0, i) -= scalar_shape_functions->JxW[qp] * scalar_shape_functions->phi[i][qp] / scalar_shape_functions->phi[0].size();
                JAe(qp, dof_i) += phi[dof_i][qp];

                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    Ce(dof_i, dof_j) += JxW[qp] * (phi[dof_i][qp] * phi[dof_j][qp] * pow(r, axisymm));
            }
        }
        Ce.scale(1.0 / rs_matrix[0]);
        break;
        // Here when you will implement loads with F
        // for (unsigned int qp = 0; qp < nq; qp++)
        //     for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
        //         Fe(dof_i) += JxW[qp] * phi[dof_i][qp]; // * js_loads[elem->id()][2];
    }
    default:
        libmesh_error_msg("Physical region not implemented");
    }
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

    switch (elem_data.element_type->get_physical_region())
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        auto &Fe = this->vectors[0][0];

        unsigned int n_dofs_var = elem_data.elem->n_nodes();

        Real js[3];

        dof_id_type dof_jx = elem_data.elem->dof_number(elem_builder_misc_data.loads_sys_number, elem_builder_misc_data.jsx_var_num, 0);
        dof_id_type dof_jy = elem_data.elem->dof_number(elem_builder_misc_data.loads_sys_number, elem_builder_misc_data.jsy_var_num, 0);
        dof_id_type dof_jz = elem_data.elem->dof_number(elem_builder_misc_data.loads_sys_number, elem_builder_misc_data.jsz_var_num, 0);

        js[0] = elem_builder_misc_data.loads_map->current_local_solution->el(dof_jx);
        js[1] = elem_builder_misc_data.loads_map->current_local_solution->el(dof_jy);
        js[2] = elem_builder_misc_data.loads_map->current_local_solution->el(dof_jz);

        for (unsigned int k = 0; k < 3; k++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                Fe(dof_i + k * n_dofs_var) *= js[k];
        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    {
        // Here you will need to implement loads with F
        break;
    }
    default:
        libmesh_error_msg("Physical region not implemented");
    }
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

    switch (elem_data.element_type->get_physical_region())
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        auto &Fe = this->vectors[0][0];

        unsigned int n_dofs_var = elem_data.elem->n_nodes();

        Real js;

        dof_id_type dof_jz = elem_data.elem->dof_number(elem_builder_misc_data.loads_sys_number, elem_builder_misc_data.jsz_var_num, 0);

        js = elem_builder_misc_data.loads_map->current_local_solution->el(dof_jz);

        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            Fe(dof_i) *= js;
        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    {
        // Here you will need to implement loads with F
        break;
    }
    default:
        libmesh_error_msg("Physical region not implemented");
    }
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

    // Get Ae
    DenseVector<Number> Ae;
    for (unsigned int k = 0; k < 3; k++)
    {
        Ae.append(solution_vector_support.sol_vectors[0][k]);
    }

    auto &AAe = this->matrices[1][1];
    DenseVector<Number> Be_q(AAe.m());
    AAe.vector_mult(Be_q, Ae);
    auto &Me = this->matrices[1][0];
    auto &Be = this->vectors[2][0];
    Be.resize(Me.m());
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

void ElementBuilder::build_solid236_nodal_solution_vector()
{
    libmesh_error_msg("build_solid236_nodal_solution_vector is not implemented yet.");
}

void ElementBuilder::build_plane13_nodal_solution_vector()
{
    // ApecMSolutionSystem *solution_system = dynamic_cast<ApecMSolutionSystem *>(this->get_parent_object());

    auto &AAe = this->matrices[1][1];
    DenseVector<Number> Be_q(AAe.m());
    AAe.vector_mult(Be_q, solution_vector_support.sol_vectors[0][0]);
    auto &Me = this->matrices[1][0];
    auto &Be = this->vectors[2][0];
    Be.resize(Me.m());
    Me.cholesky_solve(Be_q, Be);
}

void ElementBuilder::build_solid97_element_solution_vector()
{
    // Materials could be retrieved in another point

    switch (elem_data.element_type->get_physical_region())
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        auto &J_elem = this->vectors[3][0];
        dof_id_type dof_jx = elem_data.elem->dof_number(elem_builder_misc_data.loads_sys_number, elem_builder_misc_data.jsx_var_num, 0);
        dof_id_type dof_jy = elem_data.elem->dof_number(elem_builder_misc_data.loads_sys_number, elem_builder_misc_data.jsy_var_num, 0);
        dof_id_type dof_jz = elem_data.elem->dof_number(elem_builder_misc_data.loads_sys_number, elem_builder_misc_data.jsz_var_num, 0);

        J_elem(0) = elem_builder_misc_data.loads_map->current_local_solution->el(dof_jx);
        J_elem(1) = elem_builder_misc_data.loads_map->current_local_solution->el(dof_jy);
        J_elem(2) = elem_builder_misc_data.loads_map->current_local_solution->el(dof_jz);
        if (!this->elem_builder_misc_data.ramped)
        {
            break;
        }
        else
        {
            if (this->elem_builder_misc_data.old_loads_object != nullptr)
            {
                auto &J_elem_old = this->vectors[4][0];
                // Apply ramping to the current load
                J_elem_old(0) = this->elem_builder_misc_data.loads_map->old_local_solution->el(dof_jx);
                J_elem_old(1) = this->elem_builder_misc_data.loads_map->old_local_solution->el(dof_jy);
                J_elem_old(2) = this->elem_builder_misc_data.loads_map->old_local_solution->el(dof_jz);
            }
            this->ramp_vectors(3, 4, 0, 0);
        }
        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    {
        if (elem_builder_misc_data.transient_effect)
        {
            auto &J_elem = this->vectors[3][0];

            elem_builder_misc_data.materials_object->fill_rs_matrix(this->elem_data.material_num);

            auto &rs_matrix = elem_builder_misc_data.materials_object->get_rs_matrix();

            unsigned int solv_var_vector_size = solution_vector_support.sol_vectors[0].size();

            // ApecMSolutionSystem *solution_system = dynamic_cast<ApecMSolutionSystem *>(this->get_parent_object());

            fill_solution_vector_support(main_sys->older_local_solution, 0); // solution_derivative

            DenseVector<Number> Ae_dt;
            for (unsigned int k = 0; k < 3; k++)
            {
                Ae_dt.append(solution_vector_support.sol_vectors[0][k]);
            }

            auto &JAe = this->matrices[2][0];
            unsigned int j_qp_size = JAe.m();
            DenseVector<Number> J_elem_A(j_qp_size);
            JAe.vector_mult(J_elem_A, Ae_dt);

            DenseVector<Number> J_elem_V(j_qp_size);

            auto &JVe = this->matrices[2][1];
            unsigned int volt_var_num = elem_data.element_type->get_variable_information()[Globals::GlobalVariableName::VOLT].system_var_num;
            auto &Ve_dt = solution_vector_support.sol_vectors[0][volt_var_num];
            JVe.vector_mult(J_elem_V, Ve_dt);
            J_elem_A += J_elem_V;

            unsigned int n_qp = j_qp_size / 3;

            for (unsigned int k = 0; k < 3; k++)
            {
                for (unsigned int qp = 0; qp < n_qp; qp++)
                    J_elem(k) -= J_elem_A(qp + k * n_qp);
                J_elem(k) = J_elem(k) / rs_matrix[k] / n_qp;
            }
        }
        break;
    }
    default:
        libmesh_error_msg("Physical region not implemented");
    }
}

void ElementBuilder::build_solid236_element_solution_vector()
{
    // This is a placeholder for the actual implementation
    libmesh_error_msg("build_solid236_element_solution_vector is not implemented yet.");
}

void ElementBuilder::build_plane13_element_solution_vector()
{
    // Materials could be retrieved in another point

    switch (elem_data.element_type->get_physical_region())
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {

        auto &J_elem = this->vectors[3][0];
        dof_id_type dof_jz = elem_data.elem->dof_number(elem_builder_misc_data.loads_sys_number, elem_builder_misc_data.jsz_var_num, 0);

        J_elem(0) = elem_builder_misc_data.loads_map->current_local_solution->el(dof_jz);
        if (!this->elem_builder_misc_data.ramped)
        {
            break;
        }
        else
        {
            if (this->elem_builder_misc_data.old_loads_object != nullptr)
            {
                auto &J_elem_old = this->vectors[4][0];
                // Apply ramping to the current load
                J_elem_old(0) = this->elem_builder_misc_data.loads_map->old_local_solution->el(dof_jz);
            }
            // else not needed because resizing already zeros the vector
            this->ramp_vectors(3, 4, 0, 0);
        }

        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    {
        if (elem_builder_misc_data.transient_effect)
        {
            auto &J_elem = this->vectors[3][0];

            elem_builder_misc_data.materials_object->fill_rs_matrix(this->elem_data.material_num);

            auto &rs_matrix = elem_builder_misc_data.materials_object->get_rs_matrix();
            // ApecMSolutionSystem *solution_system = dynamic_cast<ApecMSolutionSystem *>(this->get_parent_object());

            fill_solution_vector_support(main_sys->request_vector(0), 0); // solution_derivative

            auto &JAe = this->matrices[2][0];

            unsigned int j_qp_size = JAe.m();

            DenseVector<Number> J_elem_A(j_qp_size);
            JAe.vector_mult(J_elem_A, solution_vector_support.sol_vectors[0][0]);

            for (unsigned int qp = 0; qp < j_qp_size; qp++)
                J_elem(0) -= J_elem_A(qp);

            J_elem(0) = J_elem(0) / rs_matrix[0] / j_qp_size;
        }
        break;
    }
    default:
        libmesh_error_msg("Physical region not implemented");
    }
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
    {
        this->fill_solution_vector_support(main_sys->current_local_solution.get(), 0);
        build_nodal_solution_vector();
        break;
    }
    case SystemType::ELEMENT_EXPLICIT:
        build_element_solution_vector();
        break;
    default:
        libmesh_error_msg("System type not supported: " + std::to_string(static_cast<unsigned int>(sys_type)));
    }
}

void ElementBuilder::build_write_solution_vectors(const std::string sys_name)
{

    auto eq_sys = this->get_parent_object();

    ApecMSolutionSystem &solution_system = eq_sys->get_system<ApecMSolutionSystem>(sys_name);
    ApecMLinearSystem &system = eq_sys->get_system<ApecMLinearSystem>("ApecMLinearSystem");

    const MeshBase &mesh = system.get_mesh();

    // Get a reference to the LinearImplicitSystem we are solving
    const auto &comm = system.comm();

    // A reference to the DofMap object for this system.  The DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the DofMap
    // in future examples.
    // const DofMap &dof_map = system.get_dof_map();
    // const DofMap &solution_dof_map = solution_system->get_dof_map();

    auto session = this->context<Session>();

    // auto time_substep = this->context<TimeStep>()->get_active_substep();

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
        libMesh::out << "Computing nodal solutions" << std::endl;
    }
    else if (sys_name == "ElementSolutionSystem")
    {
        inp_filenames.push_back("elem_sol_mtrx_" + rank_str + ".temp");
        out_filenames.push_back("elem_sol_vec_" + rank_str + ".temp");
        matrices_set = 2;
        vec_set = 3;
        libMesh::out << "Computing element solutions" << std::endl;
    }
    auto sys_type = get_system_type_from_name(sys_name);

    auto inp_streams = session->open_file_streams(emat_folderpath, inp_filenames, FileMode::Read);

    std::ifstream *solution_mtrx_input_file = dynamic_cast<std::ifstream *>(inp_streams[0].get());

    auto out_streams = session->open_file_streams(emat_folderpath, out_filenames, FileMode::Write);

    std::ofstream *solution_vec_output_file = dynamic_cast<std::ofstream *>(out_streams[0].get());

    auto element_types = this->context<ElementTypesFactory>();

    std::vector<dof_id_type> dof_indices;

    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        // Start logging the shape function initialization.
        // This is done through a simple function call with
        // the name of the event to log.

        this->set_element_data(elem);

        this->read_element_matrices(*solution_mtrx_input_file, matrices_set);

        this->prepare_solution_vectors();

        this->build_solution_vector(sys_type);

        this->write_element_vectors(*solution_vec_output_file, vec_set);
    }
    // Close file streams
    solution_mtrx_input_file->close();
    solution_vec_output_file->close();
}