#ifndef SYSTEM_BUILDER_H
#define SYSTEM_BUILDER_H

#include <string>
#include "fe_data_object.h"
#include "fem_service_abstract.h"
#include "fe_parent_object.h"
#include <libmesh/system.h>
#include "globals.h"
#include "element_types.h" // for VariableCategory

using namespace libMesh;

// SystemBuilder takes care of assigning the proper variables to the systems, depending on the element formulation
// chosen by the user. It also sets proper initial settings
class SystemBuilder
    : public System::Initialization,
      public FEMServiceAbstract<ElementTypesFactory>,
      public FEParentObjectBase<System>
{
public:
    using ContextView = TypedContextView<ElementTypesFactory>;

    // Constructor: pass context view as shared_ptr to FEMServiceAbstract
    SystemBuilder(System *sys, std::shared_ptr<ContextView> ctx_view = nullptr)
        : FEParentObjectBase<System>(sys), FEMServiceAbstract<ElementTypesFactory>(ctx_view) {}

    // Override the initialize method
    void initialize() override;

    void init_variables();
    // Example: access the ElementTypesFactory from the context view
    std::shared_ptr<ElementTypesFactory> element_type_factory() { return this->context<ElementTypesFactory>(); }
};

#endif // SYSTEM_BUILDER_H