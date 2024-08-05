#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "../include/json.hpp"
#include <map> // Include for std::map
#include "vtkCellData.h"

using namespace libMesh;

class AssemblerStatic2D : public System::Assembly
{
public:
    AssemblerStatic2D(EquationSystems &eq_sys, vtkCellData *cellData)
        : eq_sys(eq_sys), cell_data(cellData) {}

    void assemble() override;

private:
    EquationSystems &eq_sys;
    vtkCellData *cell_data;
};

class AssemblerStatic3D : public System::Assembly
{
public:
    AssemblerStatic3D(EquationSystems &eq_sys, vtkCellData *cellData)
        : eq_sys(eq_sys), cell_data(cellData) {}

    void assemble() override;

private:
    EquationSystems &eq_sys;
    vtkCellData *cell_data;
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