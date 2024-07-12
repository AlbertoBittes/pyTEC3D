#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/transient_system.h"
#include <map> // Include for std::map

using namespace libMesh;

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
    // Modified constructor to accept additional parameters
    SystemBuilder(EquationSystems &eq_systems, bool transient, bool ferromagnetic, std::string fe_family, int dim, std::string formulation)
        : equation_systems(eq_systems), transient(transient), ferromagnetic(ferromagnetic), fe_family(fe_family), dim(dim), formulation(formulation)
    {
        createSystem(); // Call createSystem at the end of the constructor
    }

    void createSystem();

    void setVariables();

    template <typename SystemType>
    SystemType &getSystemByType();
};