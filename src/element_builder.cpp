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
#include "element_types_factory.h"

//
//
// constructor and init
//
//

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

    this->build_matrices(elem_data.element_type->get_dof_matrix_idx_size_map().get(), 0);

    this->build_matrices(elem_data.element_type->get_nodal_solution_matrix_idx_size_map().get(), 1);
    this->build_matrices(elem_data.element_type->get_element_solution_matrix_idx_size_map().get(), 2);

    this->build_vectors(elem_data.element_type->get_vector_idx_size_map().get(), 0);
    // this->build_vectors(vector_idx_size_map, 1);
}

void ElementBuilder::prepare_solution_vectors()
{
    // reinit to be done inside each element
    // as well as the reception of shape functions
    // this->elem_data.element_type->reinit_fe(this->elem);
    this->build_vectors(elem_data.element_type->get_nodal_solution_vector_idx_size_map().get(), 2);
    this->build_vectors(elem_data.element_type->get_element_solution_vector_idx_size_map().get(), 3);
    this->build_vectors(elem_data.element_type->get_element_solution_vector_idx_size_map().get(), 4);
    // this->build_vectors(elem_data.element_type->get_vector_idx_size_map().get(), 0);
    // this->build_vectors(vector_idx_size_map, 1);
}

//
//
// Ramping method
//
//

void ElementBuilder::ramp_vectors(unsigned int set_vector_1, unsigned int set_vector_2, unsigned int idx_vector_1, unsigned int idx_vector_2)
{
    Real y1;
    Real y2;
    Real dt = elem_builder_misc_data.time_substep->get_parent_object()->get_delta_t();
    Real time = elem_builder_misc_data.time_substep->get_time();
    Real time_beginning = elem_builder_misc_data.time_substep->get_parent_object()->get_time_beginning();
    Real deltim_subst = time - time_beginning;
    unsigned int vec_size = vectors[set_vector_1][idx_vector_1].size();
    for (unsigned int i = 0; i < vec_size; i++)
    {
        y2 = vectors[set_vector_1][idx_vector_1](i);
        y1 = vectors[set_vector_2][idx_vector_2](i);

        vectors[set_vector_1][idx_vector_1](i) = y1 + (y2 - y1) / dt * deltim_subst;
    }
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

void ElementBuilder::set_elem_builder_misc_data()
{
    auto curr_time_step = this->context<TimeStep>();
    auto old_time_step = this->old_context<TimeStep>();
    elem_builder_misc_data.time_substep = curr_time_step ? curr_time_step->get_active_substep() : nullptr;
    elem_builder_misc_data.dt = elem_builder_misc_data.time_substep ? elem_builder_misc_data.time_substep->get_dt() : 0.0;
    elem_builder_misc_data.ramped = elem_builder_misc_data.time_substep ? elem_builder_misc_data.time_substep->get_ramped_load() : true;
    elem_builder_misc_data.transient_effect = elem_builder_misc_data.time_substep ? elem_builder_misc_data.time_substep->get_transient_effect() : false;

    // Use old and current LoadsObject for ramping, etc.
    if (elem_builder_misc_data.time_substep->get_index() == 0)
    {

        elem_builder_misc_data.old_loads_object = this->old_context<LoadsObject>();
        // elem_builder_misc_data.old_loads_map = elem_builder_misc_data.old_loads_object ? elem_builder_misc_data.old_loads_object->get_data_map() : nullptr;
        elem_builder_misc_data.loads_object = this->context<LoadsObject>();
        elem_builder_misc_data.loads_map = elem_builder_misc_data.loads_object ? elem_builder_misc_data.loads_object->get_data_map() : nullptr;

        elem_builder_misc_data.materials_object = this->context<MaterialsObject>();

        elem_builder_misc_data.boundary_info = &(elem_builder_misc_data.loads_map->get_mesh().get_boundary_info());

        if (elem_builder_misc_data.loads_map)
        {
            elem_builder_misc_data.loads_sys_number = elem_builder_misc_data.loads_map->number();
            if (elem_builder_misc_data.loads_map->has_variable("Js_x"))
                elem_builder_misc_data.jsx_var_num = elem_builder_misc_data.loads_map->variable_number("Js_x");
            if (elem_builder_misc_data.loads_map->has_variable("Js_y"))
                elem_builder_misc_data.jsy_var_num = elem_builder_misc_data.loads_map->variable_number("Js_y");
            if (elem_builder_misc_data.loads_map->has_variable("Js_z"))
                elem_builder_misc_data.jsz_var_num = elem_builder_misc_data.loads_map->variable_number("Js_z");
        }
    }
}

//
//
// Methods to build element matrices and rhs
//
//

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
        *(current_material_object->get_dataframe_map()) != *(old_material_object->get_dataframe_map()) ||
        *(current_real_constant_object->get_dataframe_map()) != *(old_real_constant_object->get_dataframe_map()))
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
    case ElementTypeEnum::INFIN110:
        this->build_infin110_nodal_solution_vector();
        break;
    case ElementTypeEnum::INFIN111:
        this->build_infin111_nodal_solution_vector();
        break;
    case ElementTypeEnum::INFIN112:
        this->build_infin112_nodal_solution_vector();
        break;
    default:
        libmesh_error_msg("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
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
    case ElementTypeEnum::INFIN110:
        this->build_infin110_element_solution_vector();
        break;
    case ElementTypeEnum::INFIN111:
        this->build_infin111_element_solution_vector();
        break;
    case ElementTypeEnum::INFIN112:
        this->build_infin112_element_solution_vector();
        break;
    default:
        libmesh_error_msg("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
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