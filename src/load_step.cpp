#include <libmesh/explicit_system.h>

#include "load_step.h"
#include "load_stepper.h"
#include "fe_data_object.h"
#include "apec_m_system_objects.h"
#include "materials_object.h"
#include "real_constants_object.h"
#include "node_loads_dof_constraints.h"
#include "element_loads_object.h"
#include "node_loads_dof_constraints.h"
#include "time_step.h"
#include "solver_configuration.h"

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
    auto raw_eq_sys = static_cast<ApecMEquationSystems *>(load_stepper->get_parent_object());
    std::shared_ptr<ApecMEquationSystems> eq_sys(raw_eq_sys, [](ApecMEquationSystems *) {}); // aliasing shared_ptr, no ownership

    ApecMLinearSystem &sys = eq_sys->get_system<ApecMLinearSystem>("ApecMLinearSystem");

    // ApecMLoadsSystem &node_loads_sys = eq_sys->get_system<ApecMLoadsSystem>("NodeLoadsSystem");
    // auto shared_node_loads_sys = std::shared_ptr<ApecMLoadsSystem>(&node_loads_sys, [](ApecMLoadsSystem *) {});

    ApecMLoadsSystem &element_loads_sys = eq_sys->get_system<ApecMLoadsSystem>("ElementLoadsSystem");
    auto shared_element_loads_sys = std::shared_ptr<ApecMLoadsSystem>(&element_loads_sys, [](ApecMLoadsSystem *) {});

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
        case LoadStepObjectType::NODE_LOADS_DOF_CONSTRAINTS_OBJECT:
            set_object(std::make_shared<NodeLoadsDOFConstraints>(file, eq_sys));
            break;
        case LoadStepObjectType::ELEMENT_LOADS_OBJECT:
            set_object(std::make_shared<ElementLoadsObject>(file, shared_element_loads_sys));
            break;
        default:
            libmesh_error_msg("Invalid load step binary file: unknown type_id " + std::to_string(type_id));
            continue;
        }
    }

    file.close();
}