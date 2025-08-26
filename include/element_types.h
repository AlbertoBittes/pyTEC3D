#ifndef ELEMENT_TYPES_H
#define ELEMENT_TYPES_H

#include <unordered_map>
#include <set>
#include <string>
#include <memory>
#include <vector>
#include <libmesh/libmesh.h> // Include the libMesh header
#include <libmesh/quadrature_gauss.h>
#include <libmesh/string_to_enum.h>
#include <libmesh/fe.h>

#include "fe_data_object.h"
#include "globals.h" // Include the header for global constants

using namespace libMesh;

class ElementType;

enum ElementTypeEnum : unsigned int
{
    SOLID97 = 97,
    SOLID236 = 236,
    PLANE233 = 233,
    INFIN110 = 110,
    INFIN111 = 111,
    INFIN112 = 112
};

enum PhysicalRegionEnum : unsigned int
{
    NON_MAGNETIC,
    SOFT_MAGNETIC,
    HARD_MAGNETIC,
    STRANDED_CONDUCTOR,
    SOLID_CONDUCTOR,
    SOLID_CONDUCTOR_CIRCUIT_COUPLED,
    STRANDED_COIL,
    STRANDED_COIL_VOLTAGE_FED,
    STRANDED_COIL_CIRCUIT_COUPLED,
    UNKNOWN
};

struct FEDataShapeFunctions
{
    Order order;
    FEFamily family;
    unsigned int dimension;
    Order quadrature_order;
    unsigned int quadrature_dim; // If 0, use default quadrature
    bool is_infinite;
    bool is_vector; // NEW: true for FEVectorBase, false for FEBase
};
typedef std::unordered_map<unsigned int, std::unique_ptr<ElementType>> ElementTypesMap;
typedef std::map<Globals::GlobalVariableName, std::set<subdomain_id_type>> VariableToElementIdsMapType;
typedef std::map<Globals::GlobalVariableName, Globals::VariableInfo> VariableInformation;

// This is probably one of the most important classes of the program (ETYPE in APDL).
// It holds the ElementTypes ID numbers (formulations) and all the information relevant for the formulation of the problem:
// The physical region represented by the element (KEYOPT1)
//
// Axisymmetric options (KEYOPT3)
// Last, it contains an array of libMesh FE objects, which hold all the shape functions, quadrature points, etc.
// These infomration will be used to build by the ElementBuilder e.g. to build the element matrices
// Also, each element has associated DOF variables (e.g. SOLI97, see APDL, can have A_x, A_y, A_z, VOLT...)
// given a mesh containing the element type IDs for each element, the SystemBuilder will distribute properly the variables
// in the proper subdomains of the mesh. The ElementType contains info about name, number of variables, if they are active, etc.
class ElementType
{
public:
    ElementType(const std::vector<unsigned int> &keyopts,
                std::function<void(std::vector<Globals::GlobalVariableName> &)> variable_collector,
                std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);

    void build_variable_information(std::function<void(std::vector<Globals::GlobalVariableName> &)> variable_collector);

    void set_variables_order();
    // Virtual methods to be implemented by derived classes
    inline std::vector<Globals::GlobalVariableName> get_active_variables(Globals::SystemType variable_system_type) const;

    void ElementType::build_fe_array(std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);
    // Getters
    inline const std::vector<unsigned int> &get_keyopts() const { return keyopts; }
    // Getter for ElementTypeData
    inline const VariableInformation &get_variable_information() const { return variable_info; }

    inline VariableInformation &get_variable_information() { return variable_info; }

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_solution_matrix_idx_size_map() = 0;

    std::vector<std::unique_ptr<FEAbstract>> fe_array; // Finite element objects

protected:
    VariableInformation variable_info;
    const std::vector<unsigned int> keyopts; // Element type keyopts

    // PhysicalRegionEnum physical_region; // Flag to indicate the physical region of the element type
    // Order order;
};

class AVElementType : public ElementType
{
public:
    // Constructor
    AVElementType::AVElementType(const std::vector<unsigned int> &keyopts,
                                 std::function<void(std::vector<Globals::GlobalVariableName> &)> variable_collector,
                                 std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);
};

class AV3DElementType : public AVElementType
{
public:
    // Constructor
    AV3DElementType::AV3DElementType(const std::vector<unsigned int> &keyopts,
                                     std::function<void(std::vector<Globals::GlobalVariableName> &)> variable_collector,
                                     std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);
};

class AV2DElementType : public AVElementType
{
public:
    // Constructor
    AV2DElementType(const std::vector<unsigned int> &keyopts,
                    std::function<void(std::vector<Globals::GlobalVariableName> &)> variable_collector,
                    std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector)
        : AVElementType(keyopts, variable_collector, fe_collector)
    {
        // Optionally call set_fe_array() or other initialization here
    }
};

class SOLID97Type : public AV3DElementType
{
public:
    // Direct constructor
    SOLID97Type(const std::vector<unsigned int> &keyopts);

    // Chain-prolonging constructor for derived types
    SOLID97Type(const std::vector<unsigned int> &keyopts,
                std::function<void(std::vector<Globals::GlobalVariableName> &)> variable_collector,
                std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);
    // Helper to add SOLID97 variables
    static void add_solid97_variables(const PhysicalRegionEnum physical_region,
                                      std::vector<Globals::GlobalVariableName> &vars);

    // Helper to add two FEDataShapeFunctions structs
    static void add_solid97_fe_structs(std::vector<FEDataShapeFunctions> &fe_structs,
                                       const std::vector<unsigned int> &keyopts);

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_solution_matrix_idx_size_map() override;
};

class SOLID236Type : public AV3DElementType
{
public:
    // Direct constructor
    SOLID236Type(const std::vector<unsigned int> &keyopts);

    // Chain-prolonging constructor for derived types
    SOLID236Type(const std::vector<unsigned int> &keyopts,
                 std::function<void(std::vector<Globals::GlobalVariableName> &)> variable_collector,
                 std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);

    // Helper to add SOLID236 variables
    static void add_solid236_variables(const PhysicalRegionEnum physical_region,
                                       std::vector<Globals::GlobalVariableName> &vars);

    // Helper to add two FEDataShapeFunctions structs
    static void add_solid236_fe_structs(std::vector<FEDataShapeFunctions> &fe_structs,
                                        const std::vector<unsigned int> &keyopts);

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_solution_matrix_idx_size_map() override;
};

class PLANE233Type : public AV2DElementType
{
public:
    // Direct constructor
    PLANE233Type(const std::vector<unsigned int> &keyopts);

    // Chain-prolonging constructor for derived types
    PLANE233Type(const std::vector<unsigned int> &keyopts,
                 std::function<void(std::vector<Globals::GlobalVariableName> &)> variable_collector,
                 std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);

    // Helper to add PLANE233 variables
    static void add_plane233_variables(const PhysicalRegionEnum physical_region,
                                       std::vector<Globals::GlobalVariableName> &vars);

    // Helper to add two FEDataShapeFunctions structs
    static void add_plane233_fe_structs(std::vector<FEDataShapeFunctions> &fe_structs,
                                        const std::vector<unsigned int> &keyopts);

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_solution_matrix_idx_size_map() override;
};

class INFIN110Type : public PLANE233Type
{
public:
    INFIN110Type(const std::vector<unsigned int> &keyopts);
};

class INFIN111Type : public SOLID97Type
{
public:
    INFIN111Type(const std::vector<unsigned int> &keyopts);
};

class INFIN112Type : public SOLID236Type
{
public:
    INFIN112Type(const std::vector<unsigned int> &keyopts);
};
#endif // ELEMENT_TYPES_H