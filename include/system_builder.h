#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "../include/json.hpp"
#include <map> // Include for std::map
#include "vtkCellData.h"

using namespace libMesh;

class AssemblerLinear : public System::Assembly
{
public:
    AssemblerLinear(EquationSystems &eq_sys, vtkCellData *cellData, bool transient, nlohmann.json config)
        : eq_sys(eq_sys), cell_data(cellData), config(config)
    {
        transient = config["Transient Analysis Settings"]["transient_analysis"].get<bool>();
        ferromagnetic = config["Finite Elements Settings"]["nonlinear_analysis_ferrom"].get<bool>();
        fe_family = config["Finite Elements Settings"]["fe_family"].get<std::string>();
        formulation = config["Finite Elements Settings"]["formulation"].get<std::string>();
        order = config["Finite Elements Settings"]["approx_order"].get<std::string>();
        dim = equation_systems.get_mesh().mesh_dimension();
    }

    void assemble() override;

private:
    EquationSystems &eq_sys;
    vtkCellData *cell_data;
    bool transient;
    bool ferromagnetic;
    std::string fe_family;
    std::string formulation;
    std::string order;
    int dim;
};

class SystemBuilder
{
public:
    EquationSystems &equation_systems;
    bool transient;
    bool ferromagnetic;
    // Assuming var1, var2, var3 are placeholders; replace with actual types and names
    std::string fe_family;
    int dim;
    std::string formulation;
    int type_num;
    std::string system_name;
    std::string order;
    const nlohmann::json &config;
    vtkCellData *cell_data;
    std::shared_ptr<System::Assembly> assembler; // Add AssemblerBase attribute

    // Modified constructor to accept additional parameters
    SystemBuilder(EquationSystems &eq_systems, const nlohmann::json &config, vtkCellData *cellData)
        : equation_systems(eq_systems), config(config), cell_data(cellData)
    {
        createSystem(); // Call createSystem at the end of the constructor
        setVariables(); // Call setVariables at the end of the constructor
    }

    void createSystem();

    void setVariables();

    // System &getSystemByType();
};