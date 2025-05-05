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
// think if this can be don with n_dofs
void ElementType::populate_all()
{
    this->populate_matrix_idx_size_map();
    this->populate_derived_matrix_idx_size_map();
    this->populate_vector_idx_size_map();
}

// SOLID97

SOLID97Type::SOLID97Type(const ElementTypeKeyopts &keyopts)
    : AV3DElementType(keyopts)
{
    element_type_data["DOFVariables"]["A_x"] =
        std::make_tuple(true, 3, false, "A", 0, "LAGRANGE", "FIRST");
    element_type_data["DOFVariables"]["A_y"] =
        std::make_tuple(true, 3, false, "A", 1, "LAGRANGE", "FIRST");
    element_type_data["DOFVariables"]["A_z"] =
        std::make_tuple(true, 3, false, "A", 2, "LAGRANGE", "FIRST");
    element_type_data["DOFVariables"]["CURR"] =
        std::make_tuple(false, 1, false, "CURR", 0, "LAGRANGE", "FIRST");
    if (keyopts.keyopts[0] == 1 || keyopts.keyopts[0] == 4 || keyopts.keyopts[0] == 5 || keyopts.keyopts[0] == 6)
        std::get<0>(element_type_data["DOFVariables"]["VOLT"]) = true; // Modify the first member
    if (keyopts.keyopts[0] == 2 || keyopts.keyopts[0] == 3)
        std::get<0>(element_type_data["DOFVariables"]["CURR"]) = true; // Modify the first member
    if (keyopts.keyopts[0] == 3)
        std::get<0>(element_type_data["DOFVariables"]["EMF"]) = true; // Modify the first member

    FEType fe_type_vec(FIRST, LAGRANGE_VEC);
    // Build a Finite Element object of the specified type.  Since the
    // FEBase::build() member dynamically creates memory we will
    // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.
    fe_vec = FEVectorBase::build(3, fe_type_vec);
    // A 2nd order Gauss quadrature rule for numerical integration.
    QGauss qrule_vec(3, fe_type_vec.default_quadrature_order());
    // Tell the finite element object to use our quadrature rule.
    fe_vec->attach_quadrature_rule(&qrule_vec);
}

void SOLID97Type::populate_matrix_idx_size_map()
{
    auto active_vars = this->getActiveVariables("DOFVariables");
    // Check this when you have periodic BCs
    matrix_idx_size_map["Ke"] = {0, active_vars.size() * fe->n_shape_functions(), active_vars.size() * fe->n_shape_functions()};
    matrix_idx_size_map["Ce"] = {1, active_vars.size() * fe->n_shape_functions(), active_vars.size() * fe->n_shape_functions()};
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
// SOLID 236

SOLID236Type::SOLID236Type(const ElementTypeKeyopts &keyopts)
    : AV3DElementType(keyopts)
{
    element_type_data["DOFVariables"]["A"] =
        std::make_tuple(true, 3, false, "A", 2, "NEDELEC_ONE", "FIRST");
    if (keyopts.keyopts[0] == 1 || keyopts.keyopts[0] == 2)
        std::get<0>(element_type_data["DOFVariables"]["VOLT"]) = true; // Modify the first member
    if (keyopts.keyopts[0] == 2)
        std::get<0>(element_type_data["DOFVariables"]["EMF"]) = true; // Modify the first member

    FEType fe_type_vec(FIRST, NEDELEC_ONE);
    // Build a Finite Element object of the specified type.  Since the
    // FEBase::build() member dynamically creates memory we will
    // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.
    fe_vec = FEVectorBase::build(3, fe_type_vec);
    // A 2nd order Gauss quadrature rule for numerical integration.
    QGauss qrule_vec(3, fe_type_vec.default_quadrature_order());
    // Tell the finite element object to use our quadrature rule.
    fe_vec->attach_quadrature_rule(&qrule_vec);
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

PLANE13Type::PLANE13Type(const ElementTypeKeyopts &keyopts)
    : AV2DElementType(keyopts)
{
    element_type_data["DOFVariables"]["A_z"] =
        std::make_tuple(true, 1, true, "A", 2, "LAGRANGE", "FIRST");

    if (keyopts.keyopts[0] == 6)
        std::get<0>(element_type_data["DOFVariables"]["VOLT"]) = true; // Modify the first member
} // Default dimension for PLANE13 is 2

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

// PLANE233

PLANE233Type::PLANE233Type(const ElementTypeKeyopts &keyopts)
    : AV2DElementType(keyopts)
{
} // Default dimension for PLANE13 is 2

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

void ElementTypesCollector::collectElementTypes(ElementTypesMapCSV &elementTypesMapCSV)
{
    for (const auto &pair : elementTypesMapCSV)
    {
        const auto &keyopts = pair.second;

        // Create the appropriate ElementType object
        if (keyopts.element_type == "SOLID97")
        {
            elementTypes[pair.first] = std::make_shared<SOLID97Type>(SOLID97Type(keyopts));
        }
        else if (keyopts.element_type == "SOLID236")
        {
            elementTypes[pair.first] = std::make_shared<SOLID236Type>(SOLID236Type(keyopts));
        }
        else if (keyopts.element_type == "PLANE13")
        {
            elementTypes[pair.first] = std::make_shared<PLANE13Type>(PLANE13Type(keyopts));
        }
        else if (keyopts.element_type == "PLANE233")
        {
            elementTypes[pair.first] = std::make_shared<PLANE233Type>(PLANE233Type(keyopts));
        }
        // Add more conditions for other ElementType subclasses as needed
    }
}

VariableToElementIdsMapType *ElementTypesCollector::associateVariablesWithElementIds(std::string variables_set) const
{
    VariableToElementIdsMapType *variableToElementIds = new VariableToElementIdsMapType();
    std::string sett = "DOFVariables";
    for (const auto &pair : elementTypes)
    {
        subdomain_id_type elementId = pair.first;
        const auto &elementType = pair.second;

        auto active_vars = elementType->getActiveVariables(sett);
        for (const auto &var : active_vars)
        {
            (*variableToElementIds)[var].insert(elementId);
        }
    }
    return variableToElementIds;
}