#ifndef ELEMENT_TYPES_H
#define ELEMENT_TYPES_H

#include <unordered_map>
#include <set>
#include <string>
#include <memory>
#include <vector>
#include <libmesh/libmesh.h> // Include the libMesh header
// Define the Finite Element object.
#include "libmesh/fe.h"

// Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"
#include "libmesh/string_to_enum.h"
#include "fe_data_object.h"

using namespace libMesh;

class ElementType;

enum ElementTypeEnum : unsigned int
{
    SOLID97 = 97,
    SOLID236 = 236,
    PLANE13 = 13,
    PLANE233 = 233
};

enum VariableCategory
{
    DOFVariables,
    NodalSolutionVariables,
    ElementSolutionVariables,
    LoadsVariables,
    NodalBCsVariables,
};

enum DOFMatrixNames
{
    Ke,
    Ce,
};

enum SolutionMatrixNames
{
    Me,
    AAe,
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
    VariableCategory category;
    bool is_active;
    unsigned int physical_components;
    bool is_scalar;
    std::string api_label;
    unsigned int api_components;
    FEFamily fe_family;
    Order order;
};

typedef std::unordered_map<unsigned int, std::vector<unsigned int>> ElementTypesMapCSV;
typedef std::unordered_map<unsigned int, std::unique_ptr<ElementType>> ElementTypesMap;
typedef std::map<std::string, std::set<subdomain_id_type>> VariableToElementIdsMapType;
//                         label         is_active, physical components, is scalar, API label, API components, fe family, order
typedef std::map<std::string, VariableData> VariableInformation;

class ElementType
{
public:
    // Constructor
    ElementType(const std::vector<unsigned int> &keyopts)
        : keyopts(keyopts)
    {
    }

    // Virtual methods to be implemented by derived classes
    inline std::vector<std::string> get_active_variables(VariableCategory variable_cat) const
    {
        std::vector<std::string> active_vars;
        for (const auto &pair : variable_info)
        {
            const std::string &var_name = pair.first;
            const VariableData &var_data = pair.second;

            if (var_data.category == variable_cat && var_data.is_active)
            {
                active_vars.push_back(var_name);
            }
        }
        return active_vars;
    }

    inline void set_conductive()
    {
        auto it = variable_info.find("VOLT");
        is_conductive_ = (it != variable_info.end()) && it->second.is_active;
    }

    inline bool is_conductive() const { return is_conductive_; }

    // Getters
    inline const std::vector<unsigned int> &get_keyopts() const { return keyopts; }
    // Getter for ElementTypeData
    inline const VariableInformation &get_variable_information() const
    {
        return variable_info;
    }
    inline unsigned int get_active_dof_vars_size() const
    {
        return active_dof_vars_size;
    }
    // const std::map<std::string, std::tuple<unsigned int, unsigned int, unsigned int>> &get_matrix_idx_size_map() const { return matrix_idx_size_map; }
    // const std::map<std::string, std::tuple<unsigned int, unsigned int, unsigned int>> &get_derived_matrix_idx_size_map() const { return derived_matrix_idx_size_map; }
    // const std::map<std::string, std::tuple<unsigned int, unsigned int>> &get_vector_idx_size_map() const { return vector_idx_size_map; }

    // Virtual methods for populating maps
    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_dof_matrix_idx_size_map() = 0;

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_solution_matrix_idx_size_map() = 0;

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_vector_idx_size_map() = 0;

    std::vector<std::unique_ptr<FEAbstract>> fe_array; // Finite element object

protected:
    const std::vector<unsigned int> keyopts; // Element type keyopts
    VariableInformation variable_info;
    unsigned int active_dof_vars_size; // Number of active DOF variables
    bool is_conductive_;               // Flag to indicate if the element type is conductive
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

    std::unique_ptr<std::vector<std::vector<unsigned int>>> get_solution_matrix_idx_size_map() override;

    std::unique_ptr<std::vector<std::vector<unsigned int>>> get_vector_idx_size_map() override;
};

class AV2DElementType : public AVElementType
{
public:
    // Constructor
    AV2DElementType(const std::vector<unsigned int> &keyopts);

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_dof_matrix_idx_size_map()
    {
        libmesh_error_msg("get_dof_matrix_idx_size_map() is not implemented for AV2DElementType");
    };

    std::unique_ptr<std::vector<std::vector<unsigned int>>> get_solution_matrix_idx_size_map()
    {
        libmesh_error_msg("get_solution_matrix_idx_size_map() is not implemented for AV2DElementType");
    }

    std::unique_ptr<std::vector<std::vector<unsigned int>>> get_vector_idx_size_map()
    {
        libmesh_error_msg("get_vector_idx_size_map() is not implemented for AV2DElementType");
    }
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
    // Factory method to create an ElementType object
    static std::unique_ptr<ElementType> create_element_type(ElementTypeEnum type,
                                                            const std::vector<unsigned int> &keyopts);
    // static ElementTypeEnum string_to_element_type_enum(const std::string &type)
    // {
    //     // Convert the string type to ElementTypeEnum
    //     if (type == "SOLID97")
    //         return ElementTypeEnum::SOLID97;
    //     else if (type == "SOLID236")
    //         return ElementTypeEnum::SOLID236;
    //     else if (type == "PLANE13")
    //         return ElementTypeEnum::PLANE13;
    //     else if (type == "PLANE233")
    //         return ElementTypeEnum::PLANE233;
    //     else
    //         libmesh_error_msg("Unknown element type");
    // }

    VariableToElementIdsMapType associate_variables_with_element_IDs(VariableCategory variable_category) const;

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

    void create_element_types(const ElementTypesMapCSV &csv_map);
};
#endif // ELEMENT_TYPES_H