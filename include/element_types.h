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

struct ElementTypeProperties
{
    std::string element_type;
    std::vector<int> keyopts;
};

struct ScalarShapeFunctions
{
    const std::vector<Real> &JxW;
    const std::vector<std::vector<Real>> &phi;
    const std::vector<std::vector<RealGradient>> &dphi;
};

struct VectorShapeFunctions
{
    const std::vector<Real> &JxW;
    const std::vector<std::vector<RealGradient>> &phi_vec;
    const std::vector<std::vector<RealTensor>> &dphi_vec;
    const std::vector<std::vector<RealGradient>> &curl_phi_vec;
};

typedef std::unordered_map<unsigned int, ElementTypeProperties> ElementTypesMapCSV;
typedef std::unordered_map<unsigned int, std::shared_ptr<ElementType>> ElementTypesMap;
typedef std::unordered_map<std::string, std::set<subdomain_id_type>> VariableToElementIdsMapType;

ElementTypesMapCSV readElementTypes(const std::string &filename);

class ElementType
{
public:
    ElementType(const ElementTypeProperties &properties, const std::string &order = "default_order", const std::string &fe_family = "default_family", unsigned int quadrature_order = 0, unsigned int dimension = 0)
        : properties(properties), order(order), fe_family(fe_family), quadrature_order(quadrature_order), dimension(dimension), az_dofs(0) {}

    virtual std::vector<std::string> getActiveVariables(std::string variables_set) const = 0;

    virtual std::vector<std::string> getVariables(std::string variables_set) const = 0;

    virtual bool is_conductive() const = 0;

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

    ScalarShapeFunctions get_scalar_shape_functions() const
    {
        // initialize the ShapeFunctions struct
        // and fill it with the shape functions from the finite element object
        // and the finite element vector object
        return ScalarShapeFunctions{fe->get_JxW(), fe->get_phi(), fe->get_dphi()};
    }

    VectorShapeFunctions get_vector_shape_functions() const
    {
        // initialize the ShapeFunctions struct
        // and fill it with the shape functions from the finite element object
        // and the finite element vector object
        return VectorShapeFunctions{fe_vec->get_JxW(), fe_vec->get_phi(), fe_vec->get_dphi(), fe_vec->get_curl_phi()};
    }

    virtual void populate_matrix_idx_size_map() = 0;
    virtual void populate_derived_matrix_idx_size_map() = 0;
    virtual void populate_vector_idx_size_map() = 0;

    ElementTypeProperties properties;
    std::string order;
    std::string fe_family;
    unsigned int quadrature_order;
    unsigned int dimension;
    unsigned int az_dofs;
    std::unique_ptr<FEBase> fe;
    std::unique_ptr<FEVectorBase> fe_vec;
    std::map<std::string, std::tuple<unsigned int, unsigned int, unsigned int>> matrix_idx_size_map;
    std::map<std::string, std::tuple<unsigned int, unsigned int, unsigned int>> derived_matrix_idx_size_map;
    std::map<std::string, std::tuple<unsigned int, unsigned int>> vector_idx_size_map;
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

    void populate_matrix_idx_size_map() override
    {
        matrix_idx_size_map["Ke"] = {0, this->getActiveVariables("DOFVariables").size() * fe->n_shape_functions(), this->getActiveVariables("DOFVariables").size() * fe->n_shape_functions()};
        matrix_idx_size_map["Ce"] = {1, this->getActiveVariables("DOFVariables").size() * fe->n_shape_functions(), this->getActiveVariables("DOFVariables").size() * fe->n_shape_functions()};
    }
    void populate_derived_matrix_idx_size_map() override
    {
        derived_matrix_idx_size_map["Me"] = {0, fe_vec->n_shape_functions(), fe_vec->n_shape_functions()};
        derived_matrix_idx_size_map["AAe"] = {1, fe_vec->n_shape_functions(), fe_vec->n_shape_functions()};
        derived_matrix_idx_size_map["JAe"] = {2, fe->n_quadrature_points(), fe_vec->n_shape_functions()};
        derived_matrix_idx_size_map["JVe"] = {3, 3 * fe->n_quadrature_points(), fe_vec->n_shape_functions()};
    }
    void populate_vector_idx_size_map() override
    {
        vector_idx_size_map["Fe"] = {0, this->getActiveVariables("DOFVariables").size() * fe->n_shape_functions()};
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

    bool is_conductive() const override
    {
        return properties.keyopts[0] == 1;
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
    bool is_conductive() const override
    {
        return properties.keyopts[0] == 1;
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

    void populate_matrix_idx_size_map() override
    {
        matrix_idx_size_map["Ke"] = {0, fe->n_shape_functions(), fe->n_shape_functions()};
        matrix_idx_size_map["Ce"] = {1, fe->n_shape_functions(), fe->n_shape_functions()};
    }
    void populate_derived_matrix_idx_size_map() override
    {
        derived_matrix_idx_size_map["Me"] = {0, 2 * fe->n_shape_functions(), 2 * fe->n_shape_functions()};
        derived_matrix_idx_size_map["AAe"] = {1, 2 * fe->n_shape_functions(), fe->n_shape_functions()};
        derived_matrix_idx_size_map["JAe"] = {2, fe->n_quadrature_points(), fe->n_shape_functions()};
    }
    void populate_vector_idx_size_map() override
    {
        vector_idx_size_map["Fe"] = {0, fe->n_shape_functions()};
    }
    std::vector<std::string> getActiveVariables(std::string variables_set) const override
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
    bool is_conductive() const override
    {
        return properties.keyopts[0] == 1;
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
    bool is_conductive() const override
    {
        return properties.keyopts[0] == 1;
    }
};

class ElementTypesCollector
{
public:
    ElementTypesCollector(ElementTypesMapCSV &elementTypesMapCSV)
    {
        collectElementTypes(elementTypesMapCSV);
    }

    const ElementTypesMap &getElementTypes() const
    {
        return elementTypes;
    }

    VariableToElementIdsMapType associateVariablesWithElementIds(std::string variables_set) const;

    void collectElementTypes(ElementTypesMapCSV elementTypesMapCSV);

    ElementTypesMap elementTypes;
};
#endif // ELEMENT_TYPES_H