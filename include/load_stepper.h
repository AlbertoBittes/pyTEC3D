#ifndef LOAD_STEPPER_H
#define LOAD_STEPPER_H

#include "load_step.h"
#include "fe_parent_object.h"
#include "fem_service_abstract.h"
// Now LoadStepper is a FEMServiceWithOldContext, using LoadStep as context
class LoadStepper : public FEParentObjectBase<EquationSystems>,
                    public FEMServiceWithOldContext<LoadStep>
{
public:
    // Constructor
    LoadStepper(EquationSystems *eq_sys = nullptr)
        : FEParentObjectBase<EquationSystems>(eq_sys)
    {
        // Initialize context views as nullptr
        this->ctx_view = nullptr;
        this->old_ctx_view = nullptr;
    }

    // Getters for index
    unsigned int get_load_step_index() const { return load_step_index; }
    unsigned int get_cumulative_substep_index() const { return cumulative_substep_index; }

    // Switch pointers and read a new LoadStep, updating context views
    void read_new_load_step(const std::string &filename)
    {
        // Create and load the new LoadStep
        auto new_load_step = std::make_shared<LoadStep>(this);
        new_load_step->read_from_file(filename);

        // Create a temporary repository and set the new LoadStep
        FEMDataRepository temp_repo;
        temp_repo.set_object<LoadStep>(new_load_step);

        // Slice a view from the temporary repository
        auto new_view = std::make_shared<ContextView>(ContextView::from(temp_repo));

        // Update context views
        this->update_context_views(new_view);

        // Update the load step index
        load_step_index++;
    }
    void update_cumulative_substep_index()
    {
        cumulative_substep_index++;
    }

private:
    unsigned int load_step_index = 0;
    unsigned int cumulative_substep_index = 0;
};

#endif // LOAD_STEPPER_H