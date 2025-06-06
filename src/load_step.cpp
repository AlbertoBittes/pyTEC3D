#include "load_step.h"
#include "fe_data_object.h"
#include "apec_m_linear_system.h"
#include "libmesh/explicit_system.h"

using namespace libMesh;
void LoadStep::read_from_file(const std::string &filename)
{
    std::cout << "Reading LoadStep data from file: " << filename << std::endl;
    std::ifstream file(filename, std::ios::binary);
    if (!file)
    {
        libmesh_error_msg("Failed to open file: " + filename);
    }

    auto *load_stepper = this->get_parent_object();
    ApecMLinearSystem &sys = load_stepper->get_parent_object()->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    MeshBase &mesh = sys.get_mesh();
    PetscLinearSolver<Number> *petsc_linear_solver =
        cast_ptr<PetscLinearSolver<Number> *>(sys.get_linear_solver());
    std::shared_ptr<LoadStep> prev_load_step = load_stepper->get_old_load_step();

    // Clear any existing data objects
    data_objects.clear();

    // Example: read the number of data objects
    uint32_t num_objects = 0;
    file.read(reinterpret_cast<char *>(&num_objects), sizeof(uint32_t));

    // For each data object, read a type identifier and construct the appropriate object
    for (size_t i = 0; i < num_objects; ++i)
    {
        // Example: read a type id (assume uint32_t)
        uint32_t type_id = 0;
        file.read(reinterpret_cast<char *>(&type_id), sizeof(uint32_t));

        auto obj_type = static_cast<LoadStepObjectType>(type_id);

        // Based on type_id, create the correct FEDataObjectBase<LoadStepper> derived object
        std::shared_ptr<FEDataObjectAbstract> obj;

        // Example switch (replace with your actual type ids and classes)
        switch (obj_type)
        {
        case LoadStepObjectType::TIME_STEP:
            set_object(std::make_shared<TimeStep>(file, prev_load_step));
            break;
        case LoadStepObjectType::PETSC_SOLVER_CONFIG:
            set_object(std::make_shared<PetscSolverConfiguration>(file, petsc_linear_solver));
            break;
        case LoadStepObjectType::MATERIALS_OBJECT:
            set_object(std::make_shared<MaterialsObject>(file));
            break;
        case LoadStepObjectType::REAL_CONSTANTS_OBJECT:
            set_object(std::make_shared<RealConstantsObject>(file));
            break;
        case LoadStepObjectType::LOADS_OBJECT:
            set_object(std::make_shared<LoadsObject>(file, mesh));
            break;
        case LoadStepObjectType::NODE_BC_OBJECT:
            set_object(std::make_shared<NodeBCObject>(file));
            break;
        default:
            libmesh_error_msg("Invalid load step binary file: unknown type_id " + std::to_string(type_id));
            continue;
        }
    }

    file.close();
}