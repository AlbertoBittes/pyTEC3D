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
    LoadsVariables
};

enum DOFMatrixNames
{
    Ke,
    Ce,
    NUM_MATRICES
};

enum SolutionMatrixNames
{
    Me,
    AAe,
    JAe,
    JVe,
    NUM_SOLUTION_MATRICES
};
enum VectorNames
{
    Fe,
    Js,
    NUM_VECTORS
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
typedef std::unordered_map<std::string, std::set<subdomain_id_type>> VariableToElementIdsMapType;
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
    std::vector<std::string> get_active_variables(VariableCategory variable_cat) const
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

    // inline bool is_conductive() const
    // {
    //     // Check if "VOLT" exists in element_type_data["DOFVariables"]
    //     auto active_vars = this->getActiveVariables("DOFVariables");
    //     for (const auto &var : active_vars)
    //     {
    //         if (var == "VOLT")
    //         {
    //             return true; // "VOLT" found, element is conductive
    //         }
    //     }
    //     return false; // "DOFVariables" or "VOLT" not found
    // }

    // Getters
    const std::vector<unsigned int> &get_keyopts() const { return keyopts; }
    // Getter for ElementTypeData
    const VariableInformation &get_variable_information() const
    {
        return variable_info;
    }
    unsigned int get_active_dof_vars_size() const
    {
        return active_dof_vars_size;
    }
    // const std::map<std::string, std::tuple<unsigned int, unsigned int, unsigned int>> &get_matrix_idx_size_map() const { return matrix_idx_size_map; }
    // const std::map<std::string, std::tuple<unsigned int, unsigned int, unsigned int>> &get_derived_matrix_idx_size_map() const { return derived_matrix_idx_size_map; }
    // const std::map<std::string, std::tuple<unsigned int, unsigned int>> &get_vector_idx_size_map() const { return vector_idx_size_map; }

    // Virtual methods for populating maps
    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_dof_matrix_idx_size_map()
    {
        // Check this when you have periodic BCs
        unsigned int n_sf = fe_array[0]->n_shape_functions();
        std::unique_ptr<std::vector<std::vector<unsigned int>>> matrix_idx_size_map =
            std::make_unique<std::vector<std::vector<unsigned int>>>(DOFMatrixNames::NUM_MATRICES, std::vector<unsigned int>(3, 0));

        (*matrix_idx_size_map)[Ke] = {0, active_dof_vars_size * n_sf, active_dof_vars_size * n_sf};
        (*matrix_idx_size_map)[Ce] = {1, active_dof_vars_size * n_sf, active_dof_vars_size * n_sf};

        return matrix_idx_size_map;
    };
    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_solution_mat_idx_size_map() = 0;
    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_vector_idx_size_map() = 0;

    // Populate all maps
    void populate_all();

    std::vector<std::unique_ptr<FEAbstract>> fe_array; // Finite element object

protected:
    const std::vector<unsigned int> keyopts; // Element type keyopts
    VariableInformation variable_info;
    unsigned int active_dof_vars_size; // Number of active DOF variables
    // Maps for matrix and vector sizes
    // std::map<std::string, std::tuple<unsigned int, unsigned int, unsigned int>> matrix_idx_size_map;
    // std::map<std::string, std::tuple<unsigned int, unsigned int, unsigned int>> derived_matrix_idx_size_map;
    // std::map<std::string, std::tuple<unsigned int, unsigned int>> vector_idx_size_map;
};

class AVElementType : public ElementType
{
public:
    // Constructor
    AVElementType(const std::vector<unsigned int> &keyopts)
        : ElementType(keyopts) // Call the parent class constructor
    {
        // Additional initialization for AVElementType (if needed)
        variable_info["B_x"] = VariableData{VariableCategory::NodalSolutionVariables, true, 3, false, "B", 0, L2_LAGRANGE, FIRST};
        variable_info["B_y"] = VariableData{VariableCategory::NodalSolutionVariables, true, 3, false, "B", 1, L2_LAGRANGE, FIRST};
        variable_info["J_z"] = VariableData{VariableCategory::ElementSolutionVariables, true, 3, false, "J", 2, MONOMIAL, CONSTANT};
        variable_info["Js"] = VariableData{VariableCategory::LoadsVariables, true, 3, false, "Js", 3, FEFamily::INVALID_FE, Order::INVALID_ORDER};
    }
};

class AV3DElementType : public AVElementType
{
public:
    // Constructor
    AV3DElementType(const std::vector<unsigned int> &keyopts) : AVElementType(keyopts) // Call the parent class constructor
    {
        // Additional initialization for AVElementType (if needed)
        variable_info["B_z"] = VariableData{VariableCategory::NodalSolutionVariables, true, 3, false, "B", 2, L2_LAGRANGE, FIRST};
        variable_info["J_x"] = VariableData{VariableCategory::ElementSolutionVariables, true, 3, false, "J", 0, MONOMIAL, CONSTANT};
        variable_info["J_y"] = VariableData{VariableCategory::ElementSolutionVariables, true, 3, false, "J", 1, MONOMIAL, CONSTANT};
        variable_info["VOLT"] = VariableData{VariableCategory::DOFVariables, false, 1, true, "VOLT", 0, LAGRANGE, FIRST};
        variable_info["EMF"] = VariableData{VariableCategory::DOFVariables, false, 1, true, "EMF", 0, LAGRANGE, FIRST};

        FEType fe_type(FIRST, LAGRANGE);
        // Build a Finite Element object of the specified type.  Since the
        // FEBase::build() member dynamically creates memory we will
        // store the object as a std::unique_ptr<FEBase>.  This can be thought
        // of as a pointer that will clean up after itself.
        fe_array.push_back(FEBase::build(3, fe_type));
        // A 2nd order Gauss quadrature rule for numerical integration.
        QGauss *qrule = new QGauss(3, fe_type.default_quadrature_order());
        // Tell the finite element object to use our quadrature rule.
        fe_array[0]->attach_quadrature_rule(qrule);
    }
};

class AV2DElementType : public AVElementType
{
public:
    // Constructor
    AV2DElementType(const std::vector<unsigned int> &keyopts) : AVElementType(keyopts) // Call the parent class constructor
    {
        variable_info["VOLT"] = VariableData{VariableCategory::DOFVariables, false, 1, true, "VOLT", 0, LAGRANGE, FIRST};
        // Additional initialization for AVElementType (if needed)
        FEType fe_type(FIRST, LAGRANGE);
        // Build a Finite Element object of the specified type.  Since the
        // FEBase::build() member dynamically creates memory we will
        // store the object as a std::unique_ptr<FEBase>.  This can be thought
        // of as a pointer that will clean up after itself.
        fe_array.push_back(FEBase::build(2, fe_type));
        // A 2nd order Gauss quadrature rule for numerical integration.
        QGauss *qrule = new QGauss(2, fe_type.default_quadrature_order());
        // Tell the finite element object to use our quadrature rule.
        fe_array[0]->attach_quadrature_rule(qrule);
    }
};

class SOLID97Type : public AV3DElementType
{
public:
    SOLID97Type(const std::vector<unsigned int> &keyopts);

    std::unique_ptr<std::vector<std::vector<unsigned int>>> get_dof_matrix_idx_size_map() override;
    std::unique_ptr<std::vector<std::vector<unsigned int>>> get_solution_mat_idx_size_map() override;
    std::unique_ptr<std::vector<std::vector<unsigned int>>> get_vector_idx_size_map() override;
};

class SOLID236Type : public AV3DElementType
{
public:
    SOLID236Type(const std::vector<unsigned int> &keyopts);

    std::unique_ptr<std::vector<std::vector<unsigned int>>> get_dof_matrix_idx_size_map() override;
    std::unique_ptr<std::vector<std::vector<unsigned int>>> get_solution_mat_idx_size_map() override;
    std::unique_ptr<std::vector<std::vector<unsigned int>>> get_vector_idx_size_map() override;
};

class PLANE13Type : public AV2DElementType
{
public:
    PLANE13Type(const std::vector<unsigned int> &keyopts);

    std::unique_ptr<std::vector<std::vector<unsigned int>>> get_dof_matrix_idx_size_map() override;
    std::unique_ptr<std::vector<std::vector<unsigned int>>> get_solution_mat_idx_size_map() override;
    std::unique_ptr<std::vector<std::vector<unsigned int>>> get_vector_idx_size_map() override;
};

class PLANE233Type : public AV2DElementType
{
public:
    PLANE233Type(const std::vector<unsigned int> &keyopts);

    std::unique_ptr<std::vector<std::vector<unsigned int>>> get_dof_matrix_idx_size_map() override;
    std::unique_ptr<std::vector<std::vector<unsigned int>>> get_solution_mat_idx_size_map() override;
    std::unique_ptr<std::vector<std::vector<unsigned int>>> get_vector_idx_size_map() override;
};

class ElementTypesCollector
{
public:
    ElementTypesCollector() {}

    ElementTypesCollector(ElementTypesMapCSV &elementTypesMapCSV)
    {
        collectElementTypes(elementTypesMapCSV);
    }

    ElementTypesMap &getElementTypes()
    {
        return elementTypes;
    }

    ElementTypesMapCSV *readElementTypes(const std::string &filename);

    VariableToElementIdsMapType associateVariablesWithElementIds(std::string variables_set) const;

    void collectElementTypes(ElementTypesMapCSV &elementTypesMapCSV);

    ElementTypesMap elementTypes;

    std::set<std::string> getLoadsLabels() const;
};

class ElementTypeFactory
{
public:
    // Factory method to create an ElementType object
    static std::unique_ptr<ElementType> create_element_type(ElementTypeEnum type, const std::vector<unsigned int> &keyopts)
    {
        // Convert the string type to ElementTypeEnum
        switch (type)
        {
        case ElementTypeEnum::SOLID97:
            return std::make_unique<SOLID97Type>(keyopts);
        case ElementTypeEnum::SOLID236:
            return std::make_unique<SOLID236Type>(keyopts);
        case ElementTypeEnum::PLANE13:
            return std::make_unique<PLANE13Type>(keyopts);
        case ElementTypeEnum::PLANE233:
            return std::make_unique<PLANE233Type>(keyopts);
        default:
            libmesh_error_msg("Unknown element type");
        }
    }
    static ElementTypeEnum string_to_element_type_enum(const std::string &type)
    {
        // Convert the string type to ElementTypeEnum
        if (type == "SOLID97")
            return ElementTypeEnum::SOLID97;
        else if (type == "SOLID236")
            return ElementTypeEnum::SOLID236;
        else if (type == "PLANE13")
            return ElementTypeEnum::PLANE13;
        else if (type == "PLANE233")
            return ElementTypeEnum::PLANE233;
        else
            libmesh_error_msg("Unknown element type");
    }
};
#endif // ELEMENT_TYPES_H