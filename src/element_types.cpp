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

SOLID97Type::SOLID97Type(const std::vector<unsigned int> &keyopts)
    : AV3DElementType(keyopts)
{

    variable_info["A_x"] = VariableData{VariableCategory::DOFVariables, true, 3, false, "A", 0, LAGRANGE, FIRST};
    variable_info["A_y"] = VariableData{VariableCategory::DOFVariables, true, 3, false, "A", 1, LAGRANGE, FIRST};
    variable_info["A_z"] = VariableData{VariableCategory::DOFVariables, true, 3, false, "A", 2, LAGRANGE, FIRST};
    variable_info["CURR"] = VariableData{VariableCategory::DOFVariables, false, 1, false, "CURR", 0, LAGRANGE, FIRST};

    if (keyopts[1] == 1 || keyopts[1] == 4 || keyopts[1] == 5 || keyopts[1] == 6)
        variable_info["VOLT"].is_active = true; // Modify the first member
    if (keyopts[1] == 2 || keyopts[1] == 3)
        variable_info["CURR"].is_active = true; // Modify the first member
    if (keyopts[1] == 3)
        variable_info["EMF"].is_active = true; // Modify the first member

    FEType fe_type_vec(FIRST, LAGRANGE_VEC);
    // Build a Finite Element object of the specified type.  Since the
    // FEBase::build() member dynamically creates memory we will
    // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.
    std::unique_ptr<FEBase> fe_vec = FEBase::build(3, fe_type_vec);
    // A 2nd order Gauss quadrature rule for numerical integration.
    QGauss *qrule_vec = new QGauss(3, fe_type_vec.default_quadrature_order());
    // Tell the finite element object to use our quadrature rule.
    fe_vec->attach_quadrature_rule(qrule_vec);

    fe_array.push_back(std::move(fe_vec));
}

std::unique_ptr<std::vector<std::vector<unsigned int>>> SOLID97Type::get_dof_matrix_idx_size_map()
{
    // Check this when you have periodic BCs
    unsigned int n_sf = fe_array[0]->n_shape_functions();
    std::unique_ptr<std::vector<std::vector<unsigned int>>> matrix_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(DOFMatrixNames::NUM_MATRICES, std::vector<unsigned int>(3, 0));

    (*matrix_idx_size_map)[0] = {0, active_dof_vars_size * n_sf, active_dof_vars_size * n_sf};
    (*matrix_idx_size_map)[1] = {1, active_dof_vars_size * n_sf, active_dof_vars_size * n_sf};

    return matrix_idx_size_map;
}
std::unique_ptr<std::vector<std::vector<unsigned int>>> SOLID97Type::get_solution_mat_idx_size_map()
{
    // Check this when you have periodic BCs
    unsigned int n_sf = fe_array[0]->n_shape_functions();
    unsigned int n_sf_vec = fe_array[1]->n_shape_functions();
    unsigned int n_qp = fe_array[0]->n_quadrature_points();
    std::unique_ptr<std::vector<std::vector<unsigned int>>> solution_matrix_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(DOFMatrixNames::NUM_MATRICES, std::vector<unsigned int>(3, 0));

    (*solution_matrix_idx_size_map)[0] = {0, n_sf_vec, n_sf_vec};
    (*solution_matrix_idx_size_map)[1] = {1, n_sf_vec, n_sf_vec};
    (*solution_matrix_idx_size_map)[2] = {2, n_qp, n_sf_vec};
    (*solution_matrix_idx_size_map)[3] = {3, 3 * n_qp, n_sf_vec};

    return solution_matrix_idx_size_map;
}
std::unique_ptr<std::vector<std::vector<unsigned int>>> SOLID97Type::get_vector_idx_size_map()
{
    unsigned int n_sf = fe_array[0]->n_shape_functions();
    std::unique_ptr<std::vector<std::vector<unsigned int>>> vector_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(DOFMatrixNames::NUM_MATRICES, std::vector<unsigned int>(2, 0));

    (*vector_idx_size_map)[0] = {0, active_dof_vars_size * n_sf};

    return vector_idx_size_map;
}
// SOLID 236

SOLID236Type::SOLID236Type(const std::vector<unsigned int> &keyopts)
    : AV3DElementType(keyopts)
{
    // variable_info["DOFVariables"]["A"] =
    //     std::make_tuple(true, 3, false, "A", 2, "NEDELEC_ONE", "FIRST");
    // if (keyopts[1] == 1 || keyopts[1] == 2)
    //     std::get<0>(variable_info["DOFVariables"]["VOLT"]) = true; // Modify the first member
    // if (keyopts[1] == 2)
    //     std::get<0>(variable_info["DOFVariables"]["EMF"]) = true; // Modify the first member

    libmesh_error_msg("SOLID236Type::SOLID236Type is not implemented.");
    // FEType fe_type_vec(FIRST, NEDELEC_ONE);
    // // Build a Finite Element object of the specified type.  Since the
    // // FEBase::build() member dynamically creates memory we will
    // // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // // of as a pointer that will clean up after itself.
    // fe_vec = FEVectorBase::build(3, fe_type_vec);
    // // A 2nd order Gauss quadrature rule for numerical integration.
    // QGauss *qrule_vec = new QGauss(3, fe_type_vec.default_quadrature_order());
    // // Tell the finite element object to use our quadrature rule.
    // fe_vec->attach_quadrature_rule(qrule_vec);

    // auto sf_vec = this->get_vector_shape_functions();
}

// Dummy implementations for SOLID236Type

std::unique_ptr<std::vector<std::vector<unsigned int>>> SOLID236Type::populate_matrix_idx_size_map()
{
    libmesh_error_msg("SOLID236Type::populate_matrix_idx_size_map is not implemented.");
}

std::unique_ptr<std::vector<std::vector<unsigned int>>> SOLID236Type::populate_derived_matrix_idx_size_map()
{
    libmesh_error_msg("SOLID236Type::populate_derived_matrix_idx_size_map is not implemented.");
}

std::unique_ptr<std::vector<std::vector<unsigned int>>> SOLID236Type::populate_vector_idx_size_map()
{
    libmesh_error_msg("SOLID236Type::populate_vector_idx_size_map is not implemented.");
}
// PLANE13

PLANE13Type::PLANE13Type(const std::vector<unsigned int> &keyopts)
    : AV2DElementType(keyopts)
{
    variable_info["A_z"] = {VariableCategory::DOFVariables, true, 1, true, "A", 2, LAGRANGE, FIRST};

    if (keyopts[1] == 6)
        variable_info["VOLT"].is_active = true; // Modify the first member
} // Default dimension for PLANE13 is 2

std::unique_ptr<std::vector<std::vector<unsigned int>>> PLANE13Type::get_dof_matrix_idx_size_map()
{
    // Check this when you have periodic BCs
    unsigned int n_sf = fe_array[0]->n_shape_functions();
    std::unique_ptr<std::vector<std::vector<unsigned int>>> matrix_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(DOFMatrixNames::NUM_MATRICES, std::vector<unsigned int>(3, 0));

    (*matrix_idx_size_map)[Ke] = {0, active_dof_vars_size * n_sf, active_dof_vars_size * n_sf};
    (*matrix_idx_size_map)[Ce] = {1, active_dof_vars_size * n_sf, active_dof_vars_size * n_sf};

    return matrix_idx_size_map;
}
std::unique_ptr<std::vector<std::vector<unsigned int>>> PLANE13Type::get_solution_mat_idx_size_map()
{
    // Check this when you have periodic BCs
    unsigned int n_sf = fe_array[0]->n_shape_functions();
    unsigned int n_sf_vec = fe_array[1]->n_shape_functions();
    unsigned int n_qp = fe_array[0]->n_quadrature_points();
    std::unique_ptr<std::vector<std::vector<unsigned int>>> solution_matrix_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(SolutionMatrixNames::NUM_SOLUTION_MATRICES, std::vector<unsigned int>(3, 0));

    (*solution_matrix_idx_size_map)[Me] = {0, n_sf_vec, n_sf_vec};
    (*solution_matrix_idx_size_map)[AAe] = {1, n_sf_vec, n_sf_vec};
    (*solution_matrix_idx_size_map)[JAe] = {2, n_qp, n_sf_vec};
    (*solution_matrix_idx_size_map)[JVe] = {3, 3 * n_qp, n_sf_vec};

    return solution_matrix_idx_size_map;
}
std::unique_ptr<std::vector<std::vector<unsigned int>>> PLANE13Type::get_vector_idx_size_map()
{
    unsigned int n_sf = fe_array[0]->n_shape_functions();
    std::unique_ptr<std::vector<std::vector<unsigned int>>> vector_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(VectorNames::NUM_VECTORS, std::vector<unsigned int>(2, 0));

    (*vector_idx_size_map)[0] = {0, active_dof_vars_size * n_sf};

    return vector_idx_size_map;
}

// PLANE233

PLANE233Type::PLANE233Type(const std::vector<unsigned int> &keyopts)
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
        if (keyopts[0] == 97)
        {
            elementTypes[pair.first] = std::make_shared<SOLID97Type>(SOLID97Type(keyopts));
        }
        else if (keyopts[0] == 236)
        {
            elementTypes[pair.first] = std::make_shared<SOLID236Type>(SOLID236Type(keyopts));
        }
        else if (keyopts[0] == 13)
        {
            elementTypes[pair.first] = std::make_shared<PLANE13Type>(PLANE13Type(keyopts));
        }
        else if (keyopts[0] == 233)
        {
            elementTypes[pair.first] = std::make_shared<PLANE233Type>(PLANE233Type(keyopts));
        }
        // Add more conditions for other ElementType subclasses as needed
    }
}

VariableToElementIdsMapType ElementTypesCollector::associateVariablesWithElementIds(std::string variables_set) const
{
    VariableToElementIdsMapType variableToElementIds;
    for (const auto &pair : elementTypes)
    {
        subdomain_id_type elementId = pair.first;
        const auto &elementType = pair.second;

        auto active_vars = elementType->get_active_variables(variables_set);
        for (const auto &var : active_vars)
        {
            (variableToElementIds)[var].insert(elementId);
        }
    }
    return variableToElementIds;
}