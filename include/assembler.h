#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include "libmesh/equation_systems.h"

using namespace libMesh;

class AssemblerLinear : public System::Assembly
{
public:
    AssemblerLinear(EquationSystems &eq_sys)
        : eq_sys(eq_sys)
    {
    }

    void assemble() override;

    bool check_update_matrix();

    void compute_derived_quantities();

private:
    EquationSystems &eq_sys;
};

#endif // ASSEMBLER_H