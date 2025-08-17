#ifndef LOAD_STEP_H
#define LOAD_STEP_H

#include <memory>
#include "fe_data_object.h"
#include "fe_parent_object.h"
#include "fem_data_repository.h"

class LoadStepper;

using namespace libMesh;

// LoadStep is a repository that holds all information regarding the LoadStep to be solved.
// The various services (e.g. ElementBuilder, Assembler) will use this data to perform their tasks
class LoadStep : public FEParentObjectBase<LoadStepper>, public FEMDataRepository, public FEDataObjectAbstract
{
public:
    enum class LoadStepObjectType : uint32_t
    {
        TIME_STEP = 0,
        PETSC_SOLVER_CONFIG,
        MATERIALS_OBJECT,
        REAL_CONSTANTS_OBJECT,
        NODE_LOADS_DOF_CONSTRAINTS_OBJECT,
        ELEMENT_LOADS_OBJECT
    };
    // Default constructor
    LoadStep(LoadStepper *load_stepper = nullptr)
        : FEParentObjectBase<LoadStepper>(load_stepper)
    {
    }

    // Implement abstract destructor
    virtual ~LoadStep() override = default;

    void read_from_file(const std::string &filename);
    // // Add a FEDataObject to the array
    // void add_data_object(const std::shared_ptr<FEParentObjectBase<LoadStepper>> &obj)
    // {
    //     data_objects.push_back(obj);
    // }
};

#endif // LOAD_STEP_H