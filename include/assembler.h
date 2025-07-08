#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include "fem_service_abstract.h"
#include "fem_data_repository.h"
#include "fe_parent_object.h"
#include "time_step.h"
#include "session.h"
#include "element_types.h"
#include <libmesh/system.h>

using namespace libMesh;

// The assembler expects a context view with TimeStep and Session (shared_ptr-based)
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

private:
    Substep *prev_substep = nullptr; // Pointer to the previous substep
};

class AssemblerSolution
    : public System::Assembly,
      public FEMServiceAbstract<ElementTypeFactory, TimeStep, Session>,
      public FEParentObjectBase<System>
{
public:
    using ContextView = TypedContextView<ElementTypeFactory, TimeStep, Session>;
    using FEParentObjectBase<System>::FEParentObjectBase;

    // Constructor accepting a shared_ptr to a context view
    AssemblerSolution(std::shared_ptr<ContextView> ctx_view = nullptr)
        : FEMServiceAbstract<ElementTypeFactory, TimeStep, Session>(ctx_view)
    {
    }

    void assemble() override {};
    void compute_derived_quantities();
};

#endif // ASSEMBLER_H