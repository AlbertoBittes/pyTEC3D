#include "../include/load_step.h"

void LoadStep::read_from_file(const std::string &filename)
{
    // Placeholder for actual file reading logic
    std::cout << "Reading LoadStep data from file: " << filename << std::endl;
    // Open the binary file
    std::ifstream file(filename, std::ios::binary);
    if (!file)
    {
        libmesh_error_msg("Failed to open file: " + filename);
    }

    TransientLinearImplicitSystem &system = equation_systems->get_system<TransientLinearImplicitSystem>("ApecMLinearSystem");

    time_step = std::make_shared<TimeStep>(this); // Example initialization
    time_step->read_data(file);
    time_step->populate_time_step(load_stepper->get_old_load_step());

    solver_config = std::make_shared<PetscSolverConfiguration>(this);
    solver_config->read_data(file);
    PetscLinearSolver<Number> *solver = dynamic_cast<PetscLinearSolver<Number> *>(system.get_linear_solver());
    // SolverType solver_type = Utility::string_to_enum<SolverType>("GMRES");
    // solver->set_solver_type(GMRES);
    PetscSolverConfiguration solver_configuration(*solver);
    solver_configuration.configure_solver();
    solver->set_solver_configuration(solver_configuration);

    materials = std::make_shared<MaterialsObject>(this);
    materials->read_data(file);

    real_constants = std::make_shared<RealConstantsObject>(this);
    real_constants->read_data(file);

    const MeshBase &mesh = equation_systems->get_mesh();
    loads = std::make_shared<LoadsObject>(mesh, this);
    loads->read_data(file);

    node_bc = std::make_shared<NodeBCObject>(this);
    node_bc->read_data(file);
    NodeBCObject *old_object = dynamic_cast<NodeBCObject *>(&system.get_constraint_object());

    if (old_object)
    {
        // Explicitly call the destructor on the old NodeBCObject
        old_object->~NodeBCObject();
        libMesh::out << "Old NodeBCObject destructor called." << std::endl;
    }
    else
    {
        libMesh::out << "Warning: Constraint object is not a NodeBCObject." << std::endl;
    }
    system.attach_constraint_object(*(node_bc));

    file.close();
}