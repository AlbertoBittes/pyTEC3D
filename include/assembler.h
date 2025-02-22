#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include "libmesh/equation_systems.h"
#include "../include/materials_properties.h"
#include "../include/json.hpp"
#include <vtkCellData.h>

using namespace libMesh;

class AssemblerLinear : public System::Assembly
{
public:
    AssemblerLinear(EquationSystems &eq_sys)
        : eq_sys(eq_sys), updated_bcs(false), update_matrix(true)
    {
    }

    void assemble() override;

    void update_rhs();

    void check_update_matrix();

    bool updated_bcs;
    bool update_matrix;

protected:
    EquationSystems &eq_sys;
};

#endif // ASSEMBLER_H