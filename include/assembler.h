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

// The AssemblerLinear:
// Reads the element matrices/vectors written to file by an element builder
// Builds the main system (DOFS: Ax, Ay, etc.) matrix and rhs according to the chosen time integration scheme/ramping/etc.
// Once the system matrix is built, the libMesh solver will take care of solving it
class AssemblerLinear
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

#endif // ASSEMBLER_H
