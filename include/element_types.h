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
// Define Gauss quadrature rules.

#include "fe_data_object.h"
#include "globals.h" // Include the header for global constants

using namespace libMesh;
// using namespace Globals;

class ElementType;

enum ElementTypeEnum : unsigned int
{
    SOLID97 = 97,
    SOLID236 = 236,
    PLANE13 = 13,
    PLANE233 = 233
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

enum DOFMatrixNames
{
    Ke,
    Ce,
};

enum NodalSolutionMatrixNames
{
    Me,
    AAe,
};
enum ElementSolutionMatrixNames
{
    JAe,
    JVe,
};
enum VectorNames
{
    Fe,
    Js,
};

struct VariableData
{
    Globals::VariableCategory category;
    bool is_active;
    Globals::APIVariableName api_label;
    FEFamily fe_family;
    Order order;
    unsigned int system_var_num = 0;
};

typedef std::unordered_map<unsigned int, std::vector<unsigned int>> ElementTypesMapCSV;
typedef std::unordered_map<unsigned int, std::unique_ptr<ElementType>> ElementTypesMap;
typedef std::map<Globals::GlobalVariableName, std::set<subdomain_id_type>> VariableToElementIdsMapType;
//                         label         is_active, physical components, is scalar, API label, API components, fe family, order
typedef std::map<Globals::GlobalVariableName, VariableData> VariableInformation;

class ElementType
{
public:
    // Constructor
    ElementType(const std::vector<unsigned int> &keyopts)
        : keyopts(keyopts)
    {
        physical_region = keyopts.empty() ? PhysicalRegionEnum::UNKNOWN : static_cast<PhysicalRegionEnum>(keyopts[1]);
    }

    // Virtual methods to be implemented by derived classes
    inline std::vector<Globals::GlobalVariableName> get_active_variables(Globals::VariableCategory variable_cat) const
    {
        std::vector<Globals::GlobalVariableName> active_vars;
        for (const auto &pair : variable_info)
        {
            const Globals::GlobalVariableName &var_name = pair.first;
            const VariableData &var_data = pair.second;

            if (var_data.category == variable_cat && var_data.is_active)
            {
                active_vars.push_back(var_name);
            }
        }
        return active_vars;
    }

    // Getters
    inline const std::vector<unsigned int> &get_keyopts() const { return keyopts; }
    // Getter for ElementTypeData
    inline const VariableInformation &get_variable_information() const
    {
        return variable_info;
    }
    inline VariableInformation &get_variable_information()
    {
        return variable_info;
    }
    inline unsigned int get_active_dof_vars_size() const
    {
        return active_dof_vars_size;
    }
    inline PhysicalRegionEnum get_physical_region() const
    {
        return physical_region;
    }
    // const std::map<std::string, std::tuple<unsigned int, unsigned int, unsigned int>> &get_matrix_idx_size_map() const { return matrix_idx_size_map; }
    // const std::map<std::string, std::tuple<unsigned int, unsigned int, unsigned int>> &get_derived_matrix_idx_size_map() const { return derived_matrix_idx_size_map; }
    // const std::map<std::string, std::tuple<unsigned int, unsigned int>> &get_vector_idx_size_map() const { return vector_idx_size_map; }

    // Virtual methods for populating maps
    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_dof_matrix_idx_size_map() = 0;

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_nodal_solution_matrix_idx_size_map() = 0;

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_element_solution_matrix_idx_size_map() = 0;

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_vector_idx_size_map() = 0;

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_nodal_solution_vector_idx_size_map() = 0;

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_element_solution_vector_idx_size_map() = 0;

    std::vector<std::unique_ptr<FEAbstract>> fe_array; // Finite element object

protected:
    const std::vector<unsigned int> keyopts; // Element type keyopts
    VariableInformation variable_info;
    unsigned int active_dof_vars_size;  // Number of active DOF variables
    PhysicalRegionEnum physical_region; // Flag to indicate the physical region of the element type
    // Maps for matrix and vector sizes
    // std::map<std::string, std::tuple<unsigned int, unsigned int, unsigned int>> matrix_idx_size_map;
    // std::map<std::string, std::tuple<unsigned int, unsigned int, unsigned int>> derived_matrix_idx_size_map;
    // std::map<std::string, std::tuple<unsigned int, unsigned int>> vector_idx_size_map;
};

class AVElementType : public ElementType
{
public:
    // Constructor
    AVElementType(const std::vector<unsigned int> &keyopts);
};

class AV3DElementType : public AVElementType
{
public:
    // Constructor
    AV3DElementType(const std::vector<unsigned int> &keyopts);

    std::unique_ptr<std::vector<std::vector<unsigned int>>> get_dof_matrix_idx_size_map() override;

    std::unique_ptr<std::vector<std::vector<unsigned int>>> get_nodal_solution_matrix_idx_size_map() override;

    std::unique_ptr<std::vector<std::vector<unsigned int>>> get_element_solution_matrix_idx_size_map() override;

    std::unique_ptr<std::vector<std::vector<unsigned int>>> get_vector_idx_size_map() override;

    std::unique_ptr<std::vector<std::vector<unsigned int>>> get_nodal_solution_vector_idx_size_map() override;

    std::unique_ptr<std::vector<std::vector<unsigned int>>> get_element_solution_vector_idx_size_map() override;
};

class AV2DElementType : public AVElementType
{
public:
    // Constructor
    AV2DElementType(const std::vector<unsigned int> &keyopts);

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_dof_matrix_idx_size_map() override;
    std::unique_ptr<std::vector<std::vector<unsigned int>>> get_nodal_solution_matrix_idx_size_map() override;
    std::unique_ptr<std::vector<std::vector<unsigned int>>> get_element_solution_matrix_idx_size_map() override;
    std::unique_ptr<std::vector<std::vector<unsigned int>>> get_vector_idx_size_map() override;
    std::unique_ptr<std::vector<std::vector<unsigned int>>> get_nodal_solution_vector_idx_size_map() override;
    std::unique_ptr<std::vector<std::vector<unsigned int>>> get_element_solution_vector_idx_size_map() override;
};

class SOLID97Type : public AV3DElementType
{
public:
    SOLID97Type(const std::vector<unsigned int> &keyopts);
};

class SOLID236Type : public AV3DElementType
{
public:
    SOLID236Type(const std::vector<unsigned int> &keyopts);
};

class PLANE13Type : public AV2DElementType
{
public:
    PLANE13Type(const std::vector<unsigned int> &keyopts);
};

class PLANE233Type : public AV2DElementType
{
public:
    PLANE233Type(const std::vector<unsigned int> &keyopts);
};

class ElementTypeFactory : public FEDataObjectMap<ElementTypesMap>
{
public:
    // Constructor: if csv_map_ptr is not null, call create_element_types
    ElementTypeFactory(ElementTypesMapCSV *csv_map_ptr = nullptr)
    {
        if (csv_map_ptr)
            create_element_types(csv_map_ptr);
    }

    // Factory method to create an ElementType object
    static std::unique_ptr<ElementType> create_element_type(ElementTypeEnum type,
                                                            const std::vector<unsigned int> &keyopts);

    VariableToElementIdsMapType associate_variables_with_element_IDs(Globals::VariableCategory variable_category) const;

    inline ElementType *get_element_type(subdomain_id_type element_id)
    {
        auto it = data_map->find(element_id);
        if (it != data_map->end())
        {
            return it->second.get();
        }
        else
        {
            libmesh_error_msg("Element type not found for subdomain ID: " + std::to_string(element_id));
        }
    }

    void create_element_types(ElementTypesMapCSV *csv_map);
};
#endif // ELEMENT_TYPES_H