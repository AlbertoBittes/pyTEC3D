#ifndef ELEMENT_TYPES_H
#define ELEMENT_TYPES_H

#include <unordered_map>
#include <set>
#include <string>
#include <memory>
#include <vector>
#include <libmesh/libmesh.h> // Include the libMesh header
using namespace libMesh;

struct ElementTypeProperties
{
    std::string element_type;
    std::vector<int> keyopts;
};

typedef std::unordered_map<unsigned int, ElementTypeProperties> ElementTypesMapCSV;

ElementTypesMapCSV readElementTypes(const std::string &filename);

class ElementType
{
public:
    ElementType(const ElementTypeProperties &properties, const std::string &order = "default_order", const std::string &fe_family = "default_family", unsigned int quadrature_order = 0, unsigned int dimension = 0)
        : properties(properties), order(order), fe_family(fe_family), quadrature_order(quadrature_order), dimension(dimension), az_dofs(0) {}

    virtual std::vector<std::string> getActiveVariables(std::string variables_set) const = 0;

    virtual std::vector<std::string> getVariables(std::string variables_set) const = 0;

    std::string getVariableFEFamily(std::string variable)
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

    std::string getVariableOrder(std::string variable)
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
    ElementTypeProperties properties;
    std::string order;
    std::string fe_family;
    unsigned int quadrature_order;
    unsigned int dimension;
    unsigned int az_dofs;
    std::unique_ptr<FEBase> fe;
    std::unique_ptr<FEVectorBase> fe_vec;
};

class SOLID97Type : public ElementType
{
public:
    SOLID97Type(const ElementTypeProperties &properties)
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

    std::vector<std::string> getActiveVariables(std::string variables_set) const override
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

    std::vector<std::string> getVariables(std::string variables_set) const override
    {
        if (variables_set == "DOFVariables")
        {
            std::vector<std::string> variables = {"A_x", "A_y", "A_z", "v"};
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
};

class SOLID236Type : public ElementType
{
public:
    SOLID236Type(const ElementTypeProperties &properties)
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

    std::vector<std::string> getActiveVariables(std::string variables_set) const override
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

    std::vector<std::string> getVariables(std::string variables_set) const override
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
};

class PLANE13Type : public ElementType
{
public:
    PLANE13Type(const ElementTypeProperties &properties)
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

    std::vector<std::string> getActiveVariables(std::string variables_set) const override
    {
        if (variables_set == "DOFVariables")
        {
            std::vector<std::string> variables = {"A_z"};
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

    std::vector<std::string> getVariables(std::string variables_set) const override
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
};

class PLANE233Type : public ElementType
{
public:
    PLANE233Type(const ElementTypeProperties &properties)
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

    std::vector<std::string> getActiveVariables(std::string variables_set) const override
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

    std::vector<std::string> getVariables(std::string variables_set) const override
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
};

class ElementTypesCollector
{
public:
    typedef std::unordered_map<unsigned int, ElementType> ElementTypesMap;
    typedef std::unordered_map<std::string, std::set<subdomain_id_type>> VariableToElementIdsMapType;

    ElementTypesCollector(ElementTypesMapCSV &elementTypesMapCSV)
    {
        collectElementTypes(elementTypesMapCSV);
    }

    const ElementTypesMap &getElementTypes() const
    {
        return elementTypes;
    }

    VariableToElementIdsMapType associateVariablesWithElementIds(std::string variables_set) const
    {
        VariableToElementIdsMapType variableToElementIds;

        for (const auto &pair : elementTypes)
        {
            subdomain_id_type elementId = pair.first;
            const auto &elementType = pair.second;

            for (const auto &variable : elementType.getActiveVariables(variables_set))
            {
                variableToElementIds[variable].insert(elementId);
            }
        }

        return variableToElementIds;
    }

private:
    void collectElementTypes(ElementTypesMapCSV elementTypesMapCSV)
    {
        for (const auto &pair : elementTypesMapCSV)
        {
            const auto &properties = pair.second;

            // Create the appropriate ElementType object
            if (properties.element_type == "SOLID97")
            {
                elementTypes[pair.first] = SOLID97Type(properties);
            }
            else if (properties.element_type == "SOLID236")
            {
                elementTypes[pair.first] = SOLID236Type(properties);
            }
            else if (properties.element_type == "PLANE13")
            {
                elementTypes[pair.first] = PLANE13Type(properties);
            }
            else if (properties.element_type == "PLANE233")
            {
                elementTypes[pair.first] = PLANE233Type(properties);
            }
            // Add more conditions for other ElementType subclasses as needed
        }
    }

    ElementTypesMap elementTypes;
};
#endif // ELEMENT_TYPES_H