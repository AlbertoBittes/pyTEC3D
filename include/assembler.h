#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include "libmesh/equation_systems.h"
#include "../include/time_stepper.h"

using namespace libMesh;

class AssemblerLinear : public System::Assembly
{
public:
    AssemblerLinear(EquationSystems &eq_sys)
        : updated_bcs(true), eq_sys(eq_sys)
    {
    }

    void assemble() override;

    bool check_update_matrix();

    void compute_derived_quantities();

private:
    EquationSystems &eq_sys;
    bool updated_bcs;
};

#endif // ASSEMBLER_H