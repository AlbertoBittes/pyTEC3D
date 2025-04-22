#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <string>
#include <memory>
#include <vector>
#include <limits> // for std::numeric_limits
#include "../include/element_types.h"

// Base ElementType class

ElementType::ElementType(const ElementTypeProperties &properties, const std::string &order = "default_order", const std::string &fe_family = "default_family", unsigned int quadrature_order = 0, unsigned int dimension = 0)
    : properties(properties), order(order), fe_family(fe_family), quadrature_order(quadrature_order), dimension(dimension), az_dofs(0) {}

std::string ElementType::getVariableFEFamily(std::string variable)
{
    auto vars = this->getVariables("DOFVariables");
    if (std::find(vars.begin(), vars.end(), variable) != vars.end())
    {
        if (variable == "V")
            return "LAGRANGE";
        else
            return fe_family;
    }
    else
    {
        auto nodal_deriv_vars = this->getVariables("NodalDerivedQuantitiesVariables");
        if (std::find(nodal_deriv_vars.begin(), nodal_deriv_vars.end(), variable) != nodal_deriv_vars.end())
        {
            return "L2_LAGRANGE";
        }
        else
        {
            auto elem_deriv_vars = this->getVariables("ElementDerivedQuantitiesVariables");
            if (std::find(elem_deriv_vars.begin(), elem_deriv_vars.end(), variable) != elem_deriv_vars.end())
            {
                return "MONOMIAL";
            }
            else
            {
                libmesh_error_msg("Invalid variable set: " + variable);
            }
        }
    }
}

std::string ElementType::getVariableOrder(std::string variable)
{
    auto vars = this->getVariables("DOFVariables");
    if (std::find(vars.begin(), vars.end(), variable) != vars.end())
    {
        return order;
    }
    else
    {
        auto nodal_deriv_vars = this->getVariables("NodalDerivedQuantitiesVariables");
        if (std::find(nodal_deriv_vars.begin(), nodal_deriv_vars.end(), variable) != nodal_deriv_vars.end())
        {
            return "FIRST";
        }
        else
        {
            auto elem_deriv_vars = this->getVariables("ElementDerivedQuantitiesVariables");
            if (std::find(elem_deriv_vars.begin(), elem_deriv_vars.end(), variable) != elem_deriv_vars.end())
            {
                return "CONSTANT";
            }
            else
            {
                libmesh_error_msg("Invalid variable set: " + variable);
            }
        }
    }
}

inline const ScalarShapeFunctions &ElementType::get_scalar_shape_functions() const
{
    // initialize the ShapeFunctions struct
    // and fill it with the shape functions from the finite element object
    // and the finite element vector object
    return ScalarShapeFunctions{fe->get_JxW(), fe->get_phi(), fe->get_dphi()};
}

inline const VectorShapeFunctions &ElementType::get_vector_shape_functions() const
{
    // initialize the ShapeFunctions struct
    // and fill it with the shape functions from the finite element object
    // and the finite element vector object
    return VectorShapeFunctions{fe_vec->get_JxW(), fe_vec->get_phi(), fe_vec->get_dphi(), fe_vec->get_curl_phi()};
}

// think if this can be don with n_dofs
void ElementType::populate_all()
{
    this->populate_matrix_idx_size_map();
    this->populate_derived_matrix_idx_size_map();
    this->populate_vector_idx_size_map();
}

// SOLID97

SOLID97Type::SOLID97Type(const ElementTypeProperties &properties)
    : ElementType(properties, "FIRST", "LAGRANGE", 0, 3)
{
    az_dofs = 3;

    FEType fe_type(Utility::string_to_enum<Order>(order), Utility::string_to_enum<FEFamily>(fe_family));
    // Build a Finite Element object of the specified type.  Since the
    // FEBase::build() member dynamically creates memory we will
    // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.
    fe = FEBase::build(dimension, fe_type);
    // A 2nd order Gauss quadrature rule for numerical integration.
    QGauss qrule(dimension, fe_type.default_quadrature_order());
    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule(&qrule);

    FEType fe_type_vec(Utility::string_to_enum<Order>(order), Utility::string_to_enum<FEFamily>(fe_family));
    // Build a Finite Element object of the specified type.  Since the
    // FEBase::build() member dynamically creates memory we will
    // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.
    fe_vec = FEVectorBase::build(dimension, fe_type_vec);
    // A 2nd order Gauss quadrature rule for numerical integration.
    QGauss qrule(dimension, fe_type_vec.default_quadrature_order());
    // Tell the finite element object to use our quadrature rule.
    fe_vec->attach_quadrature_rule(&qrule);
}

std::vector<std::string> SOLID97Type::getElementLoads() const
{
    std::vector<std::string> loads = {"Js"};
    return loads;
}

void SOLID97Type::populate_matrix_idx_size_map()
{
    // Check this when you have periodic BCs
    matrix_idx_size_map["Ke"] = {0, this->getActiveVariables("DOFVariables").size() * fe->n_shape_functions(), this->getActiveVariables("DOFVariables").size() * fe->n_shape_functions()};
    matrix_idx_size_map["Ce"] = {1, this->getActiveVariables("DOFVariables").size() * fe->n_shape_functions(), this->getActiveVariables("DOFVariables").size() * fe->n_shape_functions()};
}
void SOLID97Type::populate_derived_matrix_idx_size_map()
{
    // Check this when you have periodic BCs
    derived_matrix_idx_size_map["Me"] = {0, fe_vec->n_shape_functions(), fe_vec->n_shape_functions()};
    derived_matrix_idx_size_map["AAe"] = {1, fe_vec->n_shape_functions(), fe_vec->n_shape_functions()};
    derived_matrix_idx_size_map["JAe"] = {2, fe->n_quadrature_points(), fe_vec->n_shape_functions()};
    derived_matrix_idx_size_map["JVe"] = {3, 3 * fe->n_quadrature_points(), fe_vec->n_shape_functions()};
}
void SOLID97Type::populate_vector_idx_size_map()
{
    vector_idx_size_map["Fe"] = {0, this->getActiveVariables("DOFVariables").size() * fe->n_shape_functions()};
}
std::vector<std::string> SOLID97Type::getActiveVariables(std::string variables_set) const
{
    if (variables_set == "DOFVariables")
    {
        std::vector<std::string> variables = {"A_x", "A_y", "A_z"};
        if (properties.keyopts[0] == 1)
        {
            variables.push_back("V");
        }
        return variables;
    }
    else if (variables_set == "NodalDerivedQuantitiesVariables")
    {
        std::vector<std::string> variables = {"B_x", "B_y", "B_z"};
        return variables;
    }
    else if (variables_set == "ElementDerivedQuantitiesVariables")
    {
        std::vector<std::string> variables = {"J_x", "J_y", "J_z"};
        return variables;
    }
    else
    {
        libmesh_error_msg("Invalid variables set: " + variables_set);
    }
}

std::vector<std::string> SOLID97Type::getVariables(std::string variables_set) const
{
    if (variables_set == "DOFVariables")
    {
        std::vector<std::string> variables = {"A_x", "A_y", "A_z", "V"};
        return variables;
    }
    else if (variables_set == "NodalDerivedQuantitiesVariables")
    {
        std::vector<std::string> variables = {"B_x", "B_y", "B_z"};
        return variables;
    }
    else if (variables_set == "ElementDerivedQuantitiesVariables")
    {
        std::vector<std::string> variables = {"J_x", "J_y", "J_z"};
        return variables;
    }
    else
    {
        libmesh_error_msg("Invalid variables set: " + variables_set);
    }
}

// SOLID 236

SOLID236Type::SOLID236Type(const ElementTypeProperties &properties)
    : ElementType(properties, "FIRST", "NEDELEC_ONE", 0, 3)
{
    az_dofs = 1;

    FEType fe_type(Utility::string_to_enum<Order>(order), LAGRANGE);
    // Build a Finite Element object of the specified type.  Since the
    // FEBase::build() member dynamically creates memory we will
    // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.
    fe = FEBase::build(dimension, fe_type);
    // A 2nd order Gauss quadrature rule for numerical integration.
    QGauss qrule(dimension, fe_type.default_quadrature_order());
    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule(&qrule);

    FEType fe_type_vec(Utility::string_to_enum<Order>(order), Utility::string_to_enum<FEFamily>(fe_family));
    // Build a Finite Element object of the specified type.  Since the
    // FEBase::build() member dynamically creates memory we will
    // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.
    fe_vec = FEVectorBase::build(dimension, fe_type_vec);
    // A 2nd order Gauss quadrature rule for numerical integration.
    QGauss qrule(dimension, fe_type_vec.default_quadrature_order());
    // Tell the finite element object to use our quadrature rule.
    fe_vec->attach_quadrature_rule(&qrule);
}

std::vector<std::string> SOLID236Type::getElementLoads() const
{
    std::vector<std::string> loads = {"Js"};
    return loads;
}

std::vector<std::string> SOLID236Type::getActiveVariables(std::string variables_set) const
{
    if (variables_set == "DOFVariables")
    {
        std::vector<std::string> variables = {"Az"};
        if (properties.keyopts[0] == 1)
        {
            variables.push_back("V");
        }
        return variables;
    }
    else if (variables_set == "NodalDerivedQuantitiesVariables")
    {
        std::vector<std::string> variables = {"B_x", "B_y", "B_z"};
        return variables;
    }
    else if (variables_set == "ElementDerivedQuantitiesVariables")
    {
        std::vector<std::string> variables = {"J_x", "J_y", "J_z"};
        return variables;
    }
    else
    {
        libmesh_error_msg("Invalid variables set: " + variables_set);
    }
}

std::vector<std::string> SOLID236Type::getVariables(std::string variables_set) const
{
    if (variables_set == "DOFVariables")
    {
        std::vector<std::string> variables = {"Az", "V"};
        return variables;
    }
    else if (variables_set == "NodalDerivedQuantitiesVariables")
    {
        std::vector<std::string> variables = {"B_x", "B_y", "B_z"};
        return variables;
    }
    else if (variables_set == "ElementDerivedQuantitiesVariables")
    {
        std::vector<std::string> variables = {"J_x", "J_y", "J_z"};
        return variables;
    }
    else
    {
        libmesh_error_msg("Invalid variables set: " + variables_set);
    }
}
// Dummy implementations for SOLID236Type

void SOLID236Type::populate_matrix_idx_size_map()
{
    libmesh_error_msg("SOLID236Type::populate_matrix_idx_size_map is not implemented.");
}

void SOLID236Type::populate_derived_matrix_idx_size_map()
{
    libmesh_error_msg("SOLID236Type::populate_derived_matrix_idx_size_map is not implemented.");
}

void SOLID236Type::populate_vector_idx_size_map()
{
    libmesh_error_msg("SOLID236Type::populate_vector_idx_size_map is not implemented.");
}
// PLANE13

PLANE13Type::PLANE13Type(const ElementTypeProperties &properties)
    : ElementType(properties, "FIRST", "LAGRANGE", 0, 2)
{
    az_dofs = 1;

    FEType fe_type(Utility::string_to_enum<Order>(order), LAGRANGE);
    // Build a Finite Element object of the specified type.  Since the
    // FEBase::build() member dynamically creates memory we will
    // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.
    fe = FEBase::build(dimension, fe_type);
    // A 2nd order Gauss quadrature rule for numerical integration.
    QGauss qrule(dimension, fe_type.default_quadrature_order());
    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule(&qrule);
} // Default dimension for PLANE13 is 2

std::vector<std::string> PLANE13Type::getElementLoads() const
{
    std::vector<std::string> loads = {"Js"};
    return loads;
}

void PLANE13Type::populate_matrix_idx_size_map()
{
    matrix_idx_size_map["Ke"] = {0, fe->n_shape_functions(), fe->n_shape_functions()};
    matrix_idx_size_map["Ce"] = {1, fe->n_shape_functions(), fe->n_shape_functions()};
}
void PLANE13Type::populate_derived_matrix_idx_size_map()
{
    derived_matrix_idx_size_map["Me"] = {0, 2 * fe->n_shape_functions(), 2 * fe->n_shape_functions()};
    derived_matrix_idx_size_map["AAe"] = {1, 2 * fe->n_shape_functions(), fe->n_shape_functions()};
    derived_matrix_idx_size_map["JAe"] = {2, fe->n_quadrature_points(), fe->n_shape_functions()};
}
void PLANE13Type::populate_vector_idx_size_map()
{
    vector_idx_size_map["Fe"] = {0, fe->n_shape_functions()};
}

std::vector<std::string> PLANE13Type::getActiveVariables(std::string variables_set) const
{
    if (variables_set == "DOFVariables")
    {
        std::vector<std::string> variables = {"A_z"};
        // if (properties.keyopts[0] == 1)
        // {
        //     variables.push_back("V");
        // }
        return variables;
    }
    else if (variables_set == "NodalDerivedQuantitiesVariables")
    {
        std::vector<std::string> variables = {"B_x", "B_y"};
        return variables;
    }
    else if (variables_set == "ElementDerivedQuantitiesVariables")
    {
        std::vector<std::string> variables = {"J_z"};
        return variables;
    }
    else
    {
        libmesh_error_msg("Invalid variables set: " + variables_set);
    }
}

std::vector<std::string> PLANE13Type::getVariables(std::string variables_set) const
{
    if (variables_set == "DOFVariables")
    {
        std::vector<std::string> variables = {"A_z", "V"};
        return variables;
    }
    else if (variables_set == "NodalDerivedQuantitiesVariables")
    {
        std::vector<std::string> variables = {"B_x", "B_y"};
        return variables;
    }
    else if (variables_set == "ElementDerivedQuantitiesVariables")
    {
        std::vector<std::string> variables = {"J_z"};
        return variables;
    }
    else
    {
        libmesh_error_msg("Invalid variables set: " + variables_set);
    }
}

// PLANE233

PLANE233Type::PLANE233Type(const ElementTypeProperties &properties)
    : ElementType(properties, "FIRST", "NEDELEC_ONE", 0, 2)
{
    az_dofs = 1;

    FEType fe_type(Utility::string_to_enum<Order>(order), LAGRANGE);
    // Build a Finite Element object of the specified type.  Since the
    // FEBase::build() member dynamically creates memory we will
    // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.
    fe = FEBase::build(dimension, fe_type);
    // A 2nd order Gauss quadrature rule for numerical integration.
    QGauss qrule(dimension, fe_type.default_quadrature_order());
    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule(&qrule);
} // Default dimension for PLANE13 is 2

std::vector<std::string> PLANE233Type::getElementLoads() const
{
    std::vector<std::string> loads = {"Js"};
    return loads;
}

std::vector<std::string> PLANE233Type::getActiveVariables(std::string variables_set) const
{
    if (variables_set == "DOFVariables")
    {
        std::vector<std::string> variables = {"Az"};
        if (properties.keyopts[0] == 1)
        {
            variables.push_back("V");
        }
        return variables;
    }
    else if (variables_set == "NodalDerivedQuantitiesVariables")
    {
        std::vector<std::string> variables = {"B_x", "B_y"};
        return variables;
    }
    else if (variables_set == "ElementDerivedQuantitiesVariables")
    {
        std::vector<std::string> variables = {"J_z"};
        return variables;
    }
    else
    {
        libmesh_error_msg("Invalid variables set: " + variables_set);
    }
}

std::vector<std::string> PLANE233Type::getVariables(std::string variables_set) const
{
    if (variables_set == "DOFVariables")
    {
        std::vector<std::string> variables = {"Az", "V"};
        return variables;
    }
    else if (variables_set == "NodalDerivedQuantitiesVariables")
    {
        std::vector<std::string> variables = {"B_x", "B_y"};
        return variables;
    }
    else if (variables_set == "ElementDerivedQuantitiesVariables")
    {
        std::vector<std::string> variables = {"J_z"};
        return variables;
    }
    else
    {
        libmesh_error_msg("Invalid variables set: " + variables_set);
    }
}

void PLANE233Type::populate_matrix_idx_size_map()
{
    libmesh_error_msg("PLANE233Type::populate_matrix_idx_size_map is not implemented.");
}

void PLANE233Type::populate_derived_matrix_idx_size_map()
{
    libmesh_error_msg("PLANE233Type::populate_derived_matrix_idx_size_map is not implemented.");
}

void PLANE233Type::populate_vector_idx_size_map()
{
    libmesh_error_msg("PLANE233Type::populate_vector_idx_size_map is not implemented.");
}

// elementTpesCollector

ElementTypesMapCSV &ElementTypesCollector::readElementTypes(const std::string &filename)
{
    ElementTypesMapCSV elementTypesMap;
    std::ifstream file(filename);

    if (!file.is_open())
    {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::string line;
    std::getline(file, line); // Skip the header line

    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string ansys_element_type_str, element_type;
        std::vector<int> keyopts(9, 0); // Initialize with 9 zeros
        unsigned int ansys_element_type;

        if (!(std::getline(iss, ansys_element_type_str, ',') && std::getline(iss, element_type, ',')))
        {
            throw std::runtime_error("Error parsing line: " + line);
        }

        ansys_element_type = std::stoul(ansys_element_type_str);

        for (int i = 0; i < 9; ++i)
        {
            std::string keyopt_str;
            if (std::getline(iss, keyopt_str, ',') && !keyopt_str.empty())
            {
                keyopts[i] = std::stoi(keyopt_str);
            }
        }

        ElementTypeProperties properties = {element_type, keyopts};
        elementTypesMap[ansys_element_type] = properties;
    }

    file.close();
    return elementTypesMap;
}

void ElementTypesCollector::collectElementTypes(ElementTypesMapCSV &elementTypesMapCSV)
{
    for (const auto &pair : elementTypesMapCSV)
    {
        const auto &properties = pair.second;

        // Create the appropriate ElementType object
        if (properties.element_type == "SOLID97")
        {
            elementTypes[pair.first] = std::make_shared<SOLID97Type>(SOLID97Type(properties));
        }
        else if (properties.element_type == "SOLID236")
        {
            elementTypes[pair.first] = std::make_shared<SOLID236Type>(SOLID236Type(properties));
        }
        else if (properties.element_type == "PLANE13")
        {
            elementTypes[pair.first] = std::make_shared<PLANE13Type>(PLANE13Type(properties));
        }
        else if (properties.element_type == "PLANE233")
        {
            elementTypes[pair.first] = std::make_shared<PLANE233Type>(PLANE233Type(properties));
        }
        // Add more conditions for other ElementType subclasses as needed
    }
}

VariableToElementIdsMapType &ElementTypesCollector::associateVariablesWithElementIds(std::string variables_set) const
{
    VariableToElementIdsMapType variableToElementIds;

    for (const auto &pair : elementTypes)
    {
        subdomain_id_type elementId = pair.first;
        const auto &elementType = pair.second;

        for (const auto &variable : elementType->getActiveVariables(variables_set))
        {
            variableToElementIds[variable].insert(elementId);
        }
    }

    return variableToElementIds;
}

std::set<std::string> ElementTypesCollector::getLoadsLabels() const
{
    std::set<std::string> loads_labels;
    for (const auto &elem_type : elementTypes)
    {
        auto loads = elem_type.second->getElementLoads();
        for (const auto &load : loads)
        {
            loads_labels.insert(load);
        }
    }
    return loads_labels;
}