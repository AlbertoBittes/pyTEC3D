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
        : eq_sys(eq_sys), time_substep(nullptr), old_time_substep(nullptr), updated_bcs(true)
    {
    }

    void assemble() override;

    void set_time_substep(Substep *time_substep)
    {
        this->old_time_substep = this->time_substep;
        this->time_substep = time_substep;
    };
    Substep *get_time_substep() const { return time_substep; };
    Substep *get_old_time_substep() const { return old_time_substep; };

    bool check_update_matrix();

    void compute_derived_quantities();

protected:
    EquationSystems &eq_sys;
    Substep *time_substep;
    Substep *old_time_substep;
    bool updated_bcs;
};

#endif // ASSEMBLER_H