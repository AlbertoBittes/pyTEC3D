#include <libmesh/explicit_system.h>

#include "load_step.h"
#include "load_stepper.h"
#include "fe_data_object.h"
#include "apec_m_system_objects.h"

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
    auto eq_sys = load_stepper->get_parent_object();
    ApecMLinearSystem &sys = eq_sys->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    ApecMLoadsSystem &loads_sys = eq_sys->get_system<ApecMLoadsSystem>("ElementLoadsSystem");
    auto shared_loads_sys = std::shared_ptr<ApecMLoadsSystem>(&loads_sys, [](ApecMLoadsSystem *) {});
    PetscLinearSolver<Number> *petsc_linear_solver =
        cast_ptr<PetscLinearSolver<Number> *>(sys.get_linear_solver());
    std::shared_ptr<LoadStep> prev_load_step = load_stepper->old_context<LoadStep>();

    // Clear any existing data objects
    data_objects.clear();

    // Example: read the number of data objects
    dof_id_type num_objects = 0;
    file.read(reinterpret_cast<char *>(&num_objects), sizeof(num_objects));

    // For each data object, read a type identifier and construct the appropriate object
    for (size_t i = 0; i < num_objects; ++i)
    {
        // Example: read a type id (assume uint32_t)
        dof_id_type type_id = 0;
        file.read(reinterpret_cast<char *>(&type_id), sizeof(type_id));

        auto obj_type = static_cast<LoadStepObjectType>(type_id);

        // Based on type_id, create the correct FEParentObjectBase<LoadStepper> derived object
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
            set_object(std::make_shared<LoadsObject>(file, shared_loads_sys));
            break;
        case LoadStepObjectType::NODE_BC_OBJECT:
            set_object(std::make_shared<NodeDOFConstraints>(file));
            break;
        default:
            libmesh_error_msg("Invalid load step binary file: unknown type_id " + std::to_string(type_id));
            continue;
        }
    }

    file.close();
}