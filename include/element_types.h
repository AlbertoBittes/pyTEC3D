/**
 * @file element_types.h
 * @brief FE formulation descriptors and libMesh FE object containers (analogous to ANSYS ETYPE).
 *
 * Each ElementType instance corresponds to one ANSYS element-type entry.  It stores
 * the option array (OPTION1..N) that encodes the physical region and formulation
 * variant, and owns the libMesh FEBase / FEVectorBase objects (with their quadrature
 * rules) that evaluate shape functions during element-matrix construction.
 *
 * Variable ownership and DOF mapping live in VariablesMap; ElementType only stores
 * formulation metadata and FE evaluation objects.
 */
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
    libMesh::Order order;
    libMesh::FEFamily family;
    unsigned int dimension;
    libMesh::Order quadrature_order;
    unsigned int quadrature_dim; // If 0, use default quadrature
    bool is_infinite;
    bool is_vector; // NEW: true for FEVectorBase, false for FEBase
};
typedef std::unordered_map<unsigned int, std::unique_ptr<ElementType>> ElementTypesMap;
typedef std::map<GlobalVariableName, std::set<libMesh::subdomain_id_type>> VariableToElementIdsMapType;

/**
 * @brief FE formulation descriptor for one element type.
 *
 * Stores the OPTION array that encodes the physical-region enum and axisymmetry
 * flags, and owns a polymorphic array of FEBase / FEVectorBase objects used by
 * ElementBuilder to evaluate shape functions at quadrature points.
 *
 * Variable ownership and DOF mapping live in VariablesMap.
 */
class ElementType
{
public:
    enum class FEKind : unsigned char
    {
        Scalar,
        Vector
    };

    /** @brief Constructs the ElementType, builds quadrature rules, and populates the FE array. */
    ElementType(const std::vector<unsigned int> &OPTIONs,
                std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);

    /** @brief (Re)builds the FE object array from the provided collector callback. */
    void build_fe_array(std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);

    /** @brief Returns the OPTION array describing the formulation variant. */
    inline const std::vector<unsigned int> &get_OPTIONs() const { return OPTIONs; }

    /** @brief Returns the sizes of matrices needed for the nodal-solution post-processing step. */
    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_nodal_solution_matrix_idx_size_map() = 0;

    /** @brief Returns the sizes of matrices needed for the element-solution post-processing step. */
    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_element_solution_matrix_idx_size_map() = 0;

    /** @brief Returns a mutable reference to the scalar FE object at slot @p idx. */
    libMesh::FEBase &scalar_fe(std::size_t idx);
    /** @brief Returns a const reference to the scalar FE object at slot @p idx. */
    const libMesh::FEBase &scalar_fe(std::size_t idx) const;
    /** @brief Returns a mutable reference to the vector FE object at slot @p idx. */
    libMesh::FEVectorBase &vector_fe(std::size_t idx);
    /** @brief Returns a const reference to the vector FE object at slot @p idx. */
    const libMesh::FEVectorBase &vector_fe(std::size_t idx) const;

    std::vector<std::unique_ptr<libMesh::QBase>> quadrature_rules; ///< Owned quadrature rules (FE objects hold raw pointers into this).
    std::vector<std::unique_ptr<libMesh::FEAbstract>> fe_array;    ///< Polymorphic FE object storage.
    std::vector<FEKind> fe_kinds;                                  ///< Scalar/vector discriminator per FE slot.
    std::vector<libMesh::FEBase *> scalar_fe_slots;                ///< Non-owning scalar FE pointers (null for vector slots).
    std::vector<libMesh::FEVectorBase *> vector_fe_slots;          ///< Non-owning vector FE pointers (null for scalar slots).

protected:
    const std::vector<unsigned int> OPTIONs; ///< Formulation option codes (OPTION1..N).
};

/** @brief Base AV-formulation element type (2D and 3D variants derive from this). */
class AVElementType : public ElementType
{
public:
    AVElementType(const std::vector<unsigned int> &OPTIONs,
                  std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);
};

/** @brief Base class for all 3D AV formulations (nodal and edge). */
class AV3DElementType : public AVElementType
{
public:
    AV3DElementType(const std::vector<unsigned int> &OPTIONs,
                    std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);
};

/** @brief Base class for all 2D AV formulations (nodal). */
class AV2DElementType : public AVElementType
{
public:
    AV2DElementType(const std::vector<unsigned int> &OPTIONs,
                    std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector)
        : AVElementType(OPTIONs, fe_collector)
    {
    }
};

/** @brief 3D nodal (A-V) formulation element type using Lagrange scalar DOFs. */
class AV3DNODALType : public AV3DElementType
{
public:
    /** @brief Constructs the element type with the given OPTION array. */
    AV3DNODALType(const std::vector<unsigned int> &OPTIONs);

    /** @brief Chain constructor for derived types that extend the FE descriptor list. */
    AV3DNODALType(const std::vector<unsigned int> &OPTIONs,
                  std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);

    /** @brief Appends the standard AV3DNODAL FE descriptors (scalar Lagrange) to @p fe_structs. */
    static void add_av3dnodal_fe_structs(std::vector<FEDataShapeFunctions> &fe_structs,
                                         const std::vector<unsigned int> &OPTIONs);

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_nodal_solution_matrix_idx_size_map() override;
    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_element_solution_matrix_idx_size_map() override;
};

/** @brief 3D edge (Nedelec A-V) formulation element type. */
class AV3DEDGEType : public AV3DElementType
{
public:
    /** @brief Constructs the element type with the given OPTION array. */
    AV3DEDGEType(const std::vector<unsigned int> &OPTIONs);

    /** @brief Chain constructor for derived types that extend the FE descriptor list. */
    AV3DEDGEType(const std::vector<unsigned int> &OPTIONs,
                 std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);

    /** @brief Appends the standard AV3DEDGE FE descriptors (vector Nedelec) to @p fe_structs. */
    static void add_av3dedge_fe_structs(std::vector<FEDataShapeFunctions> &fe_structs,
                                        const std::vector<unsigned int> &OPTIONs);

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_nodal_solution_matrix_idx_size_map() override;
    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_element_solution_matrix_idx_size_map() override;
};

/** @brief 2D nodal (Az-V) formulation element type. */
class AV2DNODALType : public AV2DElementType
{
public:
    /** @brief Constructs the element type with the given OPTION array. */
    AV2DNODALType(const std::vector<unsigned int> &OPTIONs);

    /** @brief Chain constructor for derived types that extend the FE descriptor list. */
    AV2DNODALType(const std::vector<unsigned int> &OPTIONs,
                  std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector);

    /** @brief Appends the standard AV2DNODAL FE descriptors (scalar Lagrange) to @p fe_structs. */
    static void add_av2dnodal_fe_structs(std::vector<FEDataShapeFunctions> &fe_structs,
                                         const std::vector<unsigned int> &OPTIONs);

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_nodal_solution_matrix_idx_size_map() override;
    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_element_solution_matrix_idx_size_map() override;
};

#endif // ELEMENT_TYPES_H
