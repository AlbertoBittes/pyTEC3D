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
    AssemblerLinear(EquationSystems &eq_sys, nlohmann::json config)
        : eq_sys(eq_sys), config(config)
    {
        transient = config["Transient Analysis Settings"]["transient_analysis"].get<bool>();
        ferromagnetic = config["Finite Elements Settings"]["nonlinear_analysis_ferrom"].get<bool>();
        fe_family = config["Finite Elements Settings"]["fe_family"].get<std::string>();
        formulation = config["Finite Elements Settings"]["formulation"].get<std::string>();
        order = config["Finite Elements Settings"]["approx_order"].get<std::string>();
        dim = eq_sys.get_mesh().mesh_dimension();
    }

    void assemble() override;

private:
    EquationSystems &eq_sys;
    bool transient;
    bool ferromagnetic;
    std::string fe_family;
    std::string formulation;
    std::string order;
    int dim;
    nlohmann::json config;
};

#endif // ASSEMBLER_H