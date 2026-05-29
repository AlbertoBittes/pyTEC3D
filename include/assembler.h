#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include <memory>

#include <libmesh/system.h>
#include <pybind11/pybind11.h>

#include "fe_parent_object.h"

using namespace libMesh;

class ElementBuilder;
class ApecMEquationSystems;
class ApecMLinearSystem;

#if defined(__GNUC__) || defined(__clang__)
#define APEC_M_ASSEMBLER_CLASS_HIDDEN __attribute__((visibility("hidden")))
#else
#define APEC_M_ASSEMBLER_CLASS_HIDDEN
#endif

// The AssemblerLinear:
// Reads the element matrices/vectors written to file by an element builder
// Builds the main system (DOFS: Ax, Ay, etc.) matrix and rhs according to the chosen time integration scheme/ramping/etc.
// Once the system matrix is built, the libMesh solver will take care of solving it
class APEC_M_ASSEMBLER_CLASS_HIDDEN AssemblerLinear
    : public System::Assembly,
      public FEParentObjectBase<ApecMEquationSystems>
{
public:
    using FEParentObjectBase<ApecMEquationSystems>::FEParentObjectBase;

    AssemblerLinear(ApecMEquationSystems *equation_systems = nullptr);

    void assemble() override;

    bool check_update_system_matrix();
    bool check_update_element_matrices();

    void write_element_matrices_and_rhs();
    void write_nodal_solution();
    void write_element_solution();
};

#undef APEC_M_ASSEMBLER_CLASS_HIDDEN

#endif // ASSEMBLER_H
