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
#include "apec_m_system_objects.h"
#include "element_types_factory.h"
#include "globals.h"

//
//
// constructor and init
//
//

ElementBuilder::ElementBuilder(EquationSystems *equation_systems,
                               std::shared_ptr<ContextView> ctx_view,
                               std::shared_ptr<ContextView> old_ctx_view)
    : FEMServiceWithOldContext<Session, ElementTypesFactory, TimeStep, MaterialsObject, RealConstantsObject>(ctx_view, old_ctx_view),
      FEParentObjectBase<EquationSystems>(equation_systems),
      matrices(MatricesSet::NUM_MATRICES_SETS),
      vectors(VectorsSet::NUM_VECTORS_SETS)
{
}

void ElementBuilder::init()
{
    MeshBase &mesh = this->get_parent_object()->get_mesh();
    mesh_data_index.etype_num_index = mesh.get_elem_integer_index("ansys_etype");
    mesh_data_index.material_num_index = mesh.get_elem_integer_index("ansys_material_type");
    mesh_data_index.real_constants_num_index = mesh.get_elem_integer_index("ansys_real_constant");

    populate_systems_and_dofmaps();
}

//
//
// read/write matrices/vector from/to file
//
//

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

//
//
// Methods to size element matrices and vectors
//
//

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

void ElementBuilder::prepare_element_matrix_building()
{
    // reinit to be done inside each element
    // as well as the reception of shape functions
    // this->elem_data.element_type->reinit_fe(this->elem);
    unsigned int n_dofs = dof_indices.size();

    std::vector<std::vector<unsigned int>> matrix_sizes = {
        {n_dofs, n_dofs},
        {n_dofs, n_dofs}};

    // One vector of size n_dofs
    std::vector<std::vector<unsigned int>> vector_sizes = {
        {n_dofs}};

    // Usage:
    this->build_matrices(&matrix_sizes, MatricesSet::DOF_SOLUTION_MATRICES);
    this->build_vectors(&vector_sizes, VectorsSet::RHS_VECTORS);
    this->build_matrices(elem_data.element_type->get_solution_matrix_idx_size_map().get(), MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES);
    // this->build_vectors(vector_idx_size_map, 1);
}

//
//
// Setters
//
//

void ElementBuilder::set_element_data(Elem *elem)
{
    auto element_types = this->context<ElementTypesFactory>();

    elem_data.element_type = element_types->get_element_type(elem->get_extra_integer(mesh_data_index.etype_num_index));
    elem_data.material_num = elem->get_extra_integer(mesh_data_index.material_num_index);
    elem_data.real_constant_num = elem->get_extra_integer(mesh_data_index.real_constants_num_index);

    elem_data.elem = elem;
}

//
//
// Methods to build element matrices and rhs
//
//

void ElementBuilder::build_write_element_matrices_rhs()
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

    auto main_dof_map = this->dof_maps[Globals::DOF_SOLUTION_SYSTEM];
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
    // Print check results
    if (update_matrix)
    {
        libMesh::out << "Recomputing element matrices and vectors for current time step" << std::endl;
        std::vector<std::string> filenames = {
            "mtrx_" + rank_str + ".temp",
            "sol_mtrx_" + rank_str + ".temp",
            "base_rhs_" + rank_str + ".temp"};

        // Open all files using the session method
        auto streams = session->open_file_streams(emat_folderpath, filenames, FileMode::Write);

        // Downcast to use as ofstream
        std::ofstream *mtrx_output_file = dynamic_cast<std::ofstream *>(streams[0].get());
        std::ofstream *sol_mtrx_output_file = dynamic_cast<std::ofstream *>(streams[1].get());
        std::ofstream *base_rhs_output_file = dynamic_cast<std::ofstream *>(streams[2].get());

        for (const auto &elem : mesh.active_local_element_ptr_range())
        {
            // this can be moved
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
            this->write_element_matrices(*mtrx_output_file, MatricesSet::DOF_SOLUTION_MATRICES);

            this->write_element_matrices(*sol_mtrx_output_file, MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES);

            this->write_element_vectors(*base_rhs_output_file, VectorsSet::RHS_VECTORS);
        }

        session->close_file_streams(streams);
    }
    else
        libMesh::out << "Element matrices will be reused from previous time step" << std::endl;

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

        this->read_element_vectors(*base_rhs_input_file, VectorsSet::RHS_VECTORS);

        this->build_element_rhs();

        this->write_element_vectors(*rhs_output_file, VectorsSet::RHS_VECTORS);
    }
    session->close_file_streams(streams_in);
    session->close_file_streams(streams_out);
}

void ElementBuilder::check_compute_matrix()
{
    auto eq_sys = this->get_parent_object();

    ApecMLinearSystem &system = eq_sys->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    PetscLinearSolver<Number> *solver = dynamic_cast<PetscLinearSolver<Number> *>(system.get_linear_solver());

    auto old_view = this->get_old_context_view();
    // Get current and old materials/real constants from the views
    auto current_material_object = this->context<MaterialsObject>();
    auto old_material_object = old_view ? old_view->template get<MaterialsObject>() : nullptr;
    auto current_real_constant_object = this->context<RealConstantsObject>();
    auto old_real_constant_object = old_view ? old_view->template get<RealConstantsObject>() : nullptr;

    // If not the first substep, we can reuse the preconditioner (in linear system assuming linear mat props)
    if (this->context<TimeStep>()->get_active_substep()->get_index())
    {
        update_matrix = false;
        solver->reuse_preconditioner(true);
        return;
    }

    // If first time substep of first time step, generate the matrix
    if (old_view == nullptr)
    {
        update_matrix = true;
        solver->reuse_preconditioner(false);
        return;
    }
    else if (!old_material_object || !old_real_constant_object)
    {
        update_matrix = true;
        solver->reuse_preconditioner(false);
        return;
    }

    // Check if material properties/real constants have changed
    else if (
        *(current_material_object->get_dataframe_map()) != *(old_material_object->get_dataframe_map()) ||
        *(current_real_constant_object->get_dataframe_map()) != *(old_real_constant_object->get_dataframe_map()))
    {
        update_matrix = true;
        solver->reuse_preconditioner(false);
        return;
    }
    else
    {
        update_matrix = false;
        solver->reuse_preconditioner(true);
        return;
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
    case ElementTypeEnum::INFIN110:
        this->build_infin110();
        break;
    case ElementTypeEnum::INFIN111:
        this->build_infin111();
        break;
    case ElementTypeEnum::INFIN112:
        this->build_infin112();
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
    case ElementTypeEnum::INFIN110:
        this->build_infin110_rhs();
        break;
    case ElementTypeEnum::INFIN111:
        this->build_infin111_rhs();
        break;
    case ElementTypeEnum::INFIN112:
        this->build_infin112_rhs();
        break;
    default:
        libmesh_error_msg("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
    }
}

//
//
// solution calculation methods
//
//

void ElementBuilder::build_write_solution_vectors()
{
    auto eq_sys = this->get_parent_object();

    const MeshBase &mesh = eq_sys->get_mesh();

    // Get a reference to the LinearImplicitSystem we are solving
    const auto &comm = mesh.comm();

    // A reference to the DofMap object for this system.  The DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the DofMap
    // in future examples.

    auto session = this->context<Session>();

    // auto time_substep = this->context<TimeStep>()->get_active_substep();

    const auto rank_str = std::to_string(comm.rank());

    std::filesystem::path emat_folderpath = session->get_path("emat_folderpath");

    std::vector<std::string> inp_filenames = {"sol_mtrx_" + rank_str + ".temp"};

    // matrices_set = 1;
    // vec_set = 2;
    libMesh::out << "Computing nodal and elemental solutions" << std::endl;

    // inp_filenames.push_back("elem_sol_mtrx_" + rank_str + ".temp");
    // out_filenames.push_back("elem_sol_vec_" + rank_str + ".temp");
    // matrices_set = 2;
    // vec_set = 3;
    // libMesh::out << "Computing element solutions" << std::endl;

    // auto sys_type = get_key_by_predicate_one(
    //     Globals::SYSTEM_INFO_MAP,
    //     [sys_name](const Globals::SystemInfo &info)
    //     {
    //         return info.name == sys_name;
    //     });

    auto inp_streams = session->open_file_streams(emat_folderpath, inp_filenames, FileMode::Read);

    std::ifstream *solution_mtrx_input_file = dynamic_cast<std::ifstream *>(inp_streams[0].get());

    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        // Start logging the shape function initialization.
        // This is done through a simple function call with
        // the name of the event to log.

        this->set_element_data(elem);

        this->read_element_matrices(*solution_mtrx_input_file, MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES);

        this->build_solution_vectors();
    }
    // Close file streams
    solution_mtrx_input_file->close();
    // Close systems solutions
    systems[Globals::NODAL_SOLUTION]->solution->close();
    systems[Globals::ELEMENT_SOLUTION]->solution->close();
    systems[Globals::NODAL_SOLUTION]->update();
    systems[Globals::ELEMENT_SOLUTION]->update();
}

void ElementBuilder::build_solution_vectors()
{
    // Build the nodal solution vector
    // This is a placeholder for the actual implementation
    auto etype = static_cast<ElementTypeEnum>(elem_data.element_type->get_keyopts()[0]);

    switch (etype)
    {
    case ElementTypeEnum::SOLID97:
        this->build_solid97_solution_vectors();
        break;
    case ElementTypeEnum::SOLID236:
        this->build_solid236_solution_vectors();
        break;
    case ElementTypeEnum::PLANE13:
        this->build_plane13_solution_vectors();
        break;
    case ElementTypeEnum::INFIN110:
        this->build_infin110_solution_vectors();
        break;
    case ElementTypeEnum::INFIN111:
        this->build_infin111_solution_vectors();
        break;
    case ElementTypeEnum::INFIN112:
        this->build_infin112_solution_vectors();
        break;
    default:
        libmesh_error_msg("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
    }
}

void ElementBuilder::fill_solution_vector_support(const NumericVector<Number> *solution, int dof_map_index)
{
    solution_vector_support.clear();

    // Check index validity
    if (dof_map_index < 0 || dof_map_index >= static_cast<int>(dof_maps.size()))
        return;

    DofMap *dof_map = dof_maps[dof_map_index];
    // Cast dof_map_index to SystemType for SYSTEM_INFO_MAP lookup
    Globals::SystemType sys_type = static_cast<Globals::SystemType>(dof_map_index);
    auto sys_info = Globals::get_info_by_key(Globals::SYSTEM_INFO_MAP, sys_type);
    auto var_category = sys_info.category;
    auto active_vars = elem_data.element_type->get_active_variables(var_category);

    for (const auto &var_enum : active_vars)
    {
        unsigned int system_var_num = elem_data.element_type->get_variable_information()[var_enum].system_var_num;
        std::vector<dof_id_type> dof_indices_variable;
        if (dof_map)
        {
            dof_map->dof_indices(elem_data.elem, dof_indices_variable, system_var_num);
        }
        solution_vector_support.dof_indices_var.emplace_back(std::move(dof_indices_variable));

        if (solution && dof_map)
        {
            const auto &dof_indices_ref = solution_vector_support.dof_indices_var.back();
            DenseVector<Number> local_sol(dof_indices_ref.size());
            for (unsigned int i = 0; i < dof_indices_ref.size(); ++i)
            {
                local_sol(i) = solution->el(dof_indices_ref[i]);
            }
            solution_vector_support.sol_vectors.emplace_back(std::move(local_sol));
        }
    }
}

void ElementBuilder::populate_systems_and_dofmaps()
{
    // Clear vectors first
    systems.clear();
    dof_maps.clear();

    EquationSystems *eq_sys = this->get_parent_object();

    // Loop through all system info entries
    for (const auto &entry : Globals::SYSTEM_INFO_MAP)
    {
        const std::string &sys_name = entry.second.name;
        System *base_sys = &eq_sys->get_system(sys_name);
        if (base_sys)
        {
            systems.push_back(base_sys);
            dof_maps.push_back(&base_sys->get_dof_map());
        }
    }
}