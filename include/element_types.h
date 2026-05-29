#ifndef ELEMENT_TYPES_H
#define ELEMENT_TYPES_H

#include <unordered_map>
#include <map>
#include <set>
#include <string>
#include <memory>
#include <vector>
#include <utility>
#include <libmesh/libmesh.h> // Include the libMesh header
#include <libmesh/quadrature_gauss.h>
#include <libmesh/string_to_enum.h>
#include <libmesh/fe.h>

using namespace libMesh;

class ElementType;
enum class GlobalVariableName;

enum ElementTypeEnum : unsigned int
{
    AV2DNODAL,
    AV3DNODAL,
    AV3DEDGE,
    UNKNOWN
};

enum PhysicalRegionEnum : unsigned int
{
    NON_MAGNETIC,
    SOFT_MAGNETIC,
    HARD_MAGNETIC,
    STRANDED_CONDUCTOR,
    SHORT_CIRCUIT_CONDUCTOR,
    SOLID_CONDUCTOR,
    STRANDED_COIL,
    UNKNOWN_REGION
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
typedef std::map<GlobalVariableName, std::set<subdomain_id_type>> VariableToElementIdsMapType;

// This is probably one of the most important classes of the program (ETYPE in APDL).
// It holds the ElementTypes ID numbers (formulations) and all the information relevant for the formulation of the problem:
// The physical region represented by the element (OPTION1)
//
// Axisymmetric options (OPTION3)
// Last, it contains an array of libMesh FE objects, which hold all the shape functions, quadrature points, etc.
// These infomration will be used to build by the ElementBuilder e.g. to build the element matrices
// Variable ownership and mapping live in VariablesMap; ElementType only stores formulation and FE data.
class ElementType
{
public:
    ElementType(const std::vector<unsigned int> &OPTIONs,
                std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);

    void build_fe_array(std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);
    // Getters
    inline const std::vector<unsigned int> &get_OPTIONs() const { return OPTIONs; }

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_nodal_solution_matrix_idx_size_map() = 0;
    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_element_solution_matrix_idx_size_map() = 0;

    // FE objects keep non-owning qrule pointers; retain ownership here.
    std::vector<std::unique_ptr<QBase>> quadrature_rules;
    std::vector<std::unique_ptr<FEAbstract>> fe_array; // Finite element objects

protected:
    const std::vector<unsigned int> OPTIONs; // Element type OPTIONs
};

class AVElementType : public ElementType
{
public:
    // Constructor
    AVElementType(const std::vector<unsigned int> &OPTIONs,
                  std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);
};

class AV3DElementType : public AVElementType
{
public:
    // Constructor
    AV3DElementType(const std::vector<unsigned int> &OPTIONs,
                    std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);
};

class AV2DElementType : public AVElementType
{
public:
    // Constructor
    AV2DElementType(const std::vector<unsigned int> &OPTIONs,
                    std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector)
        : AVElementType(OPTIONs, fe_collector)
    {
        // Optionally call set_fe_array() or other initialization here
    }
};

class AV3DNODALType : public AV3DElementType
{
public:
    // Direct constructor
    AV3DNODALType(const std::vector<unsigned int> &OPTIONs);

    // Chain-prolonging constructor for derived types
    AV3DNODALType(const std::vector<unsigned int> &OPTIONs,
                  std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);
    // Helper to add two FEDataShapeFunctions structs
    static void add_av3dnodal_fe_structs(std::vector<FEDataShapeFunctions> &fe_structs,
                                         const std::vector<unsigned int> &OPTIONs);

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_nodal_solution_matrix_idx_size_map() override;
    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_element_solution_matrix_idx_size_map() override;
};

class AV3DEDGEType : public AV3DElementType
{
public:
    // Direct constructor
    AV3DEDGEType(const std::vector<unsigned int> &OPTIONs);

    // Chain-prolonging constructor for derived types
    AV3DEDGEType(const std::vector<unsigned int> &OPTIONs,
                 std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);

    // Helper to add two FEDataShapeFunctions structs
    static void add_av3dedge_fe_structs(std::vector<FEDataShapeFunctions> &fe_structs,
                                        const std::vector<unsigned int> &OPTIONs);

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_nodal_solution_matrix_idx_size_map() override;
    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_element_solution_matrix_idx_size_map() override;
};

class AV2DNODALType : public AV2DElementType
{
public:
    // Direct constructor
    AV2DNODALType(const std::vector<unsigned int> &OPTIONs);

    // Chain-prolonging constructor for derived types
    AV2DNODALType(const std::vector<unsigned int> &OPTIONs,
                  std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);

    // Helper to add two FEDataShapeFunctions structs
    static void add_av2dnodal_fe_structs(std::vector<FEDataShapeFunctions> &fe_structs,
                                         const std::vector<unsigned int> &OPTIONs);

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_nodal_solution_matrix_idx_size_map() override;
    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_element_solution_matrix_idx_size_map() override;
};

#endif // ELEMENT_TYPES_H
