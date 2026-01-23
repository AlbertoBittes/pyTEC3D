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

#include "fe_data_object.h"
using namespace libMesh;

class ElementType;
enum class GlobalVariableName;

enum ElementTypeEnum : unsigned int
{
    SOLID97 = 97,
    PLANE233 = 233,
    SOLID236 = 236
};

enum PhysicalRegionEnum : unsigned int
{
    NON_MAGNETIC,
    SOFT_MAGNETIC,
    HARD_MAGNETIC,
    STRANDED_CONDUCTOR,
    SHORT_CIRCUIT_CONDUCTOR,
    SOLID_CONDUCTOR,
    SOLID_CONDUCTOR_SOLENOIDAL,
    STRANDED_COIL,
    STRANDED_COIL_SOLENOIDAL,
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
typedef std::map<GlobalVariableName, std::set<subdomain_id_type>> VariableToElementIdsMapType;

// This is probably one of the most important classes of the program (ETYPE in APDL).
// It holds the ElementTypes ID numbers (formulations) and all the information relevant for the formulation of the problem:
// The physical region represented by the element (KEYOPT1)
//
// Axisymmetric options (KEYOPT3)
// Last, it contains an array of libMesh FE objects, which hold all the shape functions, quadrature points, etc.
// These infomration will be used to build by the ElementBuilder e.g. to build the element matrices
// Variable ownership and mapping live in VariablesMap; ElementType only stores formulation and FE data.
class ElementType
{
public:
    ElementType(const std::vector<unsigned int> &keyopts,
                std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);

    void build_fe_array(std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);
    // Getters
    inline const std::vector<unsigned int> &get_keyopts() const { return keyopts; }

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_solution_matrix_idx_size_map(bool nodal_solution) = 0;

    std::vector<std::unique_ptr<FEAbstract>> fe_array; // Finite element objects

protected:
    const std::vector<unsigned int> keyopts; // Element type keyopts
};

class AVElementType : public ElementType
{
public:
    // Constructor
    AVElementType(const std::vector<unsigned int> &keyopts,
                  std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);
};

class AV3DElementType : public AVElementType
{
public:
    // Constructor
    AV3DElementType(const std::vector<unsigned int> &keyopts,
                    std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);
};

class AV2DElementType : public AVElementType
{
public:
    // Constructor
    AV2DElementType(const std::vector<unsigned int> &keyopts,
                    std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector)
        : AVElementType(keyopts, fe_collector)
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
                std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);
    // Helper to add two FEDataShapeFunctions structs
    static void add_solid97_fe_structs(std::vector<FEDataShapeFunctions> &fe_structs,
                                       const std::vector<unsigned int> &keyopts);

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_solution_matrix_idx_size_map(bool nodal_solution) override;
};

class SOLID236Type : public AV3DElementType
{
public:
    // Direct constructor
    SOLID236Type(const std::vector<unsigned int> &keyopts);

    // Chain-prolonging constructor for derived types
    SOLID236Type(const std::vector<unsigned int> &keyopts,
                 std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);

    // Helper to add two FEDataShapeFunctions structs
    static void add_solid236_fe_structs(std::vector<FEDataShapeFunctions> &fe_structs,
                                        const std::vector<unsigned int> &keyopts);

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_solution_matrix_idx_size_map(bool nodal_solution) override;
};

class PLANE233Type : public AV2DElementType
{
public:
    // Direct constructor
    PLANE233Type(const std::vector<unsigned int> &keyopts);

    // Chain-prolonging constructor for derived types
    PLANE233Type(const std::vector<unsigned int> &keyopts,
                 std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);

    // Helper to add two FEDataShapeFunctions structs
    static void add_plane233_fe_structs(std::vector<FEDataShapeFunctions> &fe_structs,
                                        const std::vector<unsigned int> &keyopts);

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_solution_matrix_idx_size_map(bool nodal_solution) override;
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

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_solution_matrix_idx_size_map(bool nodal_solution) override;
};

class INFIN112Type : public SOLID236Type
{
public:
    INFIN112Type(const std::vector<unsigned int> &keyopts);
};
#endif // ELEMENT_TYPES_H
