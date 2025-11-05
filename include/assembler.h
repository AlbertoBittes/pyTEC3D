#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include <memory>

#include "fem_service_abstract.h"
#include "fem_data_repository.h"
#include "fe_parent_object.h"
#include "time_step.h"
#include "session.h"
#include <libmesh/system.h>

using namespace libMesh;

class ElementBuilder;
class ApecMEquationSystems;
class ApecMLinearSystem;

// The AssemblerLinear:
// Reads the element matrices/vectors written to file by an element builder
// Builds the main system (DOFS: Ax, Ay, etc.) matrix and rhs according to the chosen time integration scheme/ramping/etc.
// Once the system matrix is built, the libMesh solver will take care of solving it
class AssemblerLinear
    : public System::Assembly,
      public FEMServiceWithOldContext<TimeStep, Session>,
      public FEParentObjectBase<ApecMEquationSystems>
{
public:
    using ContextView = TypedContextView<TimeStep, Session>;
    using FEParentObjectBase<ApecMEquationSystems>::FEParentObjectBase;

    AssemblerLinear(ApecMEquationSystems *equation_systems = nullptr,
                    std::shared_ptr<ContextView> ctx_view = nullptr);

    void assemble() override;

    bool check_update_system_matrix();
    bool check_update_element_matrices();

    void write_element_matrices_and_rhs();
    void write_solution_vectors();

private:
    std::unique_ptr<ElementBuilder> make_element_builder() const;
};

#endif // ASSEMBLER_H
