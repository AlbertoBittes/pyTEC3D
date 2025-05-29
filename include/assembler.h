#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include "libmesh/equation_systems.h"
#include "element_types.h" // Include the ElementTypesMap header file
#include "fem_service_abstract.h"
#include "fe_data_object.h" // Include the FEDataObjectBase header file
#include "time_step.h"

using namespace libMesh;

class AssemblerLinear : public System::Assembly, public FEMServiceAbstract, public FEDataObjectBase<System>
{
public:
    using FEDataObjectBase<System>::FEDataObjectBase;

    void assemble() override;

    // void build_write_element_matrices();

    bool check_update_matrix(Substep *time_substep);

    // bool check_compute_matrix();

    // void compute_derived_quantities();

    // void update_rhs();

private:
    Substep *prev_substep = nullptr; // Pointer to the previous substep
};

#endif // ASSEMBLER_H