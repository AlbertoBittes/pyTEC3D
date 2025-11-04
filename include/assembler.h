#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include "fem_service_abstract.h"
#include "fem_data_repository.h"
#include "fe_parent_object.h"
#include "time_step.h"
#include "session.h"
#include "element_types.h"
#include "element_types_factory.h"
#include <libmesh/system.h>

using namespace libMesh;

class ElementBuilder;

// The AssemblerLinear:
// Reads the element matrices/vectors written to file by an element builder
// Builds the main system (DOFS: Ax, Ay, etc.) matrix and rhs according to the chosen time integration scheme/ramping/etc.
// Once the system matrix is built, the libMesh solver will take care of solving it
class AssemblerLinear
    : public System::Assembly,
      public FEMServiceAbstract<TimeStep, Session>,
      public FEParentObjectBase<System>
{
public:
    using ContextView = TypedContextView<TimeStep, Session>;
    using FEParentObjectBase<System>::FEParentObjectBase;

    // Constructor accepting a shared_ptr to a context view
    AssemblerLinear(std::shared_ptr<ContextView> ctx_view = nullptr)
        : FEMServiceAbstract<TimeStep, Session>(ctx_view)
    {
    }

    void assemble() override;

    bool check_update_matrix(Substep *time_substep);

    void write_element_matrices_and_rhs(ElementBuilder &builder);
    void write_solution_vectors(ElementBuilder &builder);

private:
    Substep *prev_substep = nullptr; // Pointer to the previous substep
};

#endif // ASSEMBLER_H
