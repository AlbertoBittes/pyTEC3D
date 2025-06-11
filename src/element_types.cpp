#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <string>
#include <memory>
#include <vector>
#include <limits> // for std::numeric_limits
#include "element_types.h"

using namespace Globals;
AVElementType::AVElementType(const std::vector<unsigned int> &keyopts)
    : ElementType(keyopts) // Call the parent class constructor
{
    // Additional initialization for AVElementType (if needed)
    variable_info[GlobalVariableName::B_x] = VariableData{VariableCategory::NodalSolutionVariables, true, APIVariableName::B, L2_LAGRANGE, FIRST};
    variable_info[GlobalVariableName::B_y] = VariableData{VariableCategory::NodalSolutionVariables, true, APIVariableName::B, L2_LAGRANGE, FIRST};
    variable_info[GlobalVariableName::J_z] = VariableData{VariableCategory::ElementSolutionVariables, true, APIVariableName::J, MONOMIAL, CONSTANT};
    variable_info[GlobalVariableName::Js_z] = VariableData{VariableCategory::LoadsVariables, true, APIVariableName::Js, MONOMIAL, CONSTANT};
}

AV2DElementType::AV2DElementType(const std::vector<unsigned int> &keyopts) : AVElementType(keyopts) // Call the parent class constructor
{
    variable_info[GlobalVariableName::VOLT] = VariableData{VariableCategory::DOFVariables, false, APIVariableName::VOLT, LAGRANGE, FIRST};
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

AV3DElementType::AV3DElementType(const std::vector<unsigned int> &keyopts) : AVElementType(keyopts) // Call the parent class constructor
{
    // Additional initialization for AVElementType (if needed)
    variable_info[GlobalVariableName::B_z] = VariableData{VariableCategory::NodalSolutionVariables, true, APIVariableName::B, L2_LAGRANGE, FIRST};
    variable_info[GlobalVariableName::J_x] = VariableData{VariableCategory::ElementSolutionVariables, true, APIVariableName::J, MONOMIAL, CONSTANT};
    variable_info[GlobalVariableName::J_y] = VariableData{VariableCategory::ElementSolutionVariables, true, APIVariableName::J, MONOMIAL, CONSTANT};
    variable_info[GlobalVariableName::VOLT] = VariableData{VariableCategory::DOFVariables, false, APIVariableName::VOLT, LAGRANGE, FIRST};
    variable_info[GlobalVariableName::EMF] = VariableData{VariableCategory::DOFVariables, false, APIVariableName::EMF, LAGRANGE, FIRST};
    variable_info[GlobalVariableName::Js_x] = VariableData{VariableCategory::LoadsVariables, true, APIVariableName::Js, MONOMIAL, CONSTANT};
    variable_info[GlobalVariableName::Js_y] = VariableData{VariableCategory::LoadsVariables, true, APIVariableName::Js, MONOMIAL, CONSTANT};

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

std::unique_ptr<std::vector<std::vector<unsigned int>>> AV3DElementType::get_dof_matrix_idx_size_map()
{
    // Check this when you have periodic BCs
    unsigned int n_sf = fe_array[0]->n_shape_functions();
    std::unique_ptr<std::vector<std::vector<unsigned int>>> matrix_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(2, std::vector<unsigned int>(2, 0));

    (*matrix_idx_size_map)[Ke] = {active_dof_vars_size * n_sf, active_dof_vars_size * n_sf};
    (*matrix_idx_size_map)[Ce] = {active_dof_vars_size * n_sf, active_dof_vars_size * n_sf};

    return matrix_idx_size_map;
};

std::unique_ptr<std::vector<std::vector<unsigned int>>> AV3DElementType::get_nodal_solution_matrix_idx_size_map()
{
    // Check this when you have periodic BCs
    unsigned int n_sf = fe_array[0]->n_shape_functions();
    unsigned int n_sf_vec = fe_array[1]->n_shape_functions();
    unsigned int n_qp = fe_array[0]->n_quadrature_points();
    std::unique_ptr<std::vector<std::vector<unsigned int>>> solution_matrix_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(2, std::vector<unsigned int>(2, 0));

    (*solution_matrix_idx_size_map)[0] = {n_sf_vec, n_sf_vec};
    (*solution_matrix_idx_size_map)[1] = {n_sf_vec, n_sf_vec};

    return solution_matrix_idx_size_map;
}

std::unique_ptr<std::vector<std::vector<unsigned int>>> AV3DElementType::get_element_solution_matrix_idx_size_map()
{
    // Check this when you have periodic BCs
    unsigned int n_sf = fe_array[0]->n_shape_functions();
    unsigned int n_sf_vec = fe_array[1]->n_shape_functions();
    unsigned int n_qp = fe_array[0]->n_quadrature_points();
    std::unique_ptr<std::vector<std::vector<unsigned int>>> solution_matrix_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(2, std::vector<unsigned int>(2, 0));

    (*solution_matrix_idx_size_map)[0] = {n_qp, n_sf_vec};
    (*solution_matrix_idx_size_map)[1] = {3 * n_qp, n_sf_vec};

    return solution_matrix_idx_size_map;
}

std::unique_ptr<std::vector<std::vector<unsigned int>>> AV3DElementType::get_vector_idx_size_map()
{
    unsigned int n_sf = fe_array[0]->n_shape_functions();
    std::unique_ptr<std::vector<std::vector<unsigned int>>> vector_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(1, std::vector<unsigned int>(1, 0));

    (*vector_idx_size_map)[0] = {active_dof_vars_size * n_sf};

    return vector_idx_size_map;
}
// SOLID97

SOLID97Type::SOLID97Type(const std::vector<unsigned int> &keyopts)
    : AV3DElementType(keyopts)
{

    variable_info[GlobalVariableName::A_x] = VariableData{VariableCategory::DOFVariables, true, APIVariableName::A, LAGRANGE, FIRST};
    variable_info[GlobalVariableName::A_y] = VariableData{VariableCategory::DOFVariables, true, APIVariableName::A, LAGRANGE, FIRST};
    variable_info[GlobalVariableName::A_z] = VariableData{VariableCategory::DOFVariables, true, APIVariableName::A, LAGRANGE, FIRST};
    variable_info[GlobalVariableName::CURR] = VariableData{VariableCategory::DOFVariables, false, APIVariableName::CURR, LAGRANGE, FIRST};

    if (keyopts[1] == 1 || keyopts[1] == 4 || keyopts[1] == 5 || keyopts[1] == 6)
        variable_info[GlobalVariableName::VOLT].is_active = true; // Modify the first member
    if (keyopts[1] == 2 || keyopts[1] == 3)
        variable_info[GlobalVariableName::CURR].is_active = true; // Modify the first member
    if (keyopts[1] == 3)
        variable_info[GlobalVariableName::EMF].is_active = true; // Modify the first member

    FEType fe_type_vec(FIRST, LAGRANGE_VEC);
    // Build a Finite Element object of the specified type.  Since the
    // FEBase::build() member dynamically creates memory we will
    // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.
    std::unique_ptr<FEVectorBase> fe_vec = FEVectorBase::build(3, fe_type_vec);
    // A 2nd order Gauss quadrature rule for numerical integration.
    QGauss *qrule_vec = new QGauss(3, fe_type_vec.default_quadrature_order());
    // Tell the finite element object to use our quadrature rule.
    fe_vec->attach_quadrature_rule(qrule_vec);

    fe_array.push_back(std::move(fe_vec));

    active_dof_vars_size = 0;
    for (const auto &pair : variable_info)
    {
        const auto &var_data = pair.second;
        if (var_data.is_active && var_data.category == VariableCategory::DOFVariables)
            ++active_dof_vars_size;
    }
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

PLANE233Type::PLANE233Type(const std::vector<unsigned int> &keyopts)
    : AV2DElementType(keyopts)
{
    // TODO: Add initialization if needed
}
// Dummy implementations for SOLID236Type

// std::unique_ptr<std::vector<std::vector<unsigned int>>> SOLID236Type::populate_matrix_idx_size_map()
// {
//     libmesh_error_msg("SOLID236Type::populate_matrix_idx_size_map is not implemented.");
// }

// std::unique_ptr<std::vector<std::vector<unsigned int>>> SOLID236Type::populate_derived_matrix_idx_size_map()
// {
//     libmesh_error_msg("SOLID236Type::populate_derived_matrix_idx_size_map is not implemented.");
// }

// std::unique_ptr<std::vector<std::vector<unsigned int>>> SOLID236Type::populate_vector_idx_size_map()
// {
//     libmesh_error_msg("SOLID236Type::populate_vector_idx_size_map is not implemented.");
// }
// PLANE13

PLANE13Type::PLANE13Type(const std::vector<unsigned int> &keyopts)
    : AV2DElementType(keyopts)
{
    variable_info[GlobalVariableName::A_z] = {VariableCategory::DOFVariables, true, APIVariableName::A, LAGRANGE, FIRST};

    if (keyopts[1] == 6)
        variable_info[GlobalVariableName::VOLT].is_active = true; // Modify the first member

    active_dof_vars_size = 0;
    for (const auto &pair : variable_info)
    {
        const auto &var_data = pair.second;
        if (var_data.is_active && var_data.category == VariableCategory::DOFVariables)
            ++active_dof_vars_size;
    }
} // Default dimension for PLANE13 is 2

// std::unique_ptr<std::vector<std::vector<unsigned int>>> PLANE13Type::get_dof_matrix_idx_size_map()
// {
//     // Check this when you have periodic BCs
//     unsigned int n_sf = fe_array[0]->n_shape_functions();
//     std::unique_ptr<std::vector<std::vector<unsigned int>>> matrix_idx_size_map =
//         std::make_unique<std::vector<std::vector<unsigned int>>>(DOFMatrixNames::NUM_MATRICES, std::vector<unsigned int>(3, 0));

//     (*matrix_idx_size_map)[Ke] = {0, active_dof_vars_size * n_sf, active_dof_vars_size * n_sf};
//     (*matrix_idx_size_map)[Ce] = {1, active_dof_vars_size * n_sf, active_dof_vars_size * n_sf};

//     return matrix_idx_size_map;
// }
// std::unique_ptr<std::vector<std::vector<unsigned int>>> PLANE13Type::get_solution_matrix_idx_size_map()
// {
//     // Check this when you have periodic BCs
//     unsigned int n_sf = fe_array[0]->n_shape_functions();
//     unsigned int n_sf_vec = fe_array[1]->n_shape_functions();
//     unsigned int n_qp = fe_array[0]->n_quadrature_points();
//     std::unique_ptr<std::vector<std::vector<unsigned int>>> solution_matrix_idx_size_map =
//         std::make_unique<std::vector<std::vector<unsigned int>>>(SolutionMatrixNames::NUM_SOLUTION_MATRICES, std::vector<unsigned int>(3, 0));

//     (*solution_matrix_idx_size_map)[Me] = {0, n_sf_vec, n_sf_vec};
//     (*solution_matrix_idx_size_map)[AAe] = {1, n_sf_vec, n_sf_vec};
//     (*solution_matrix_idx_size_map)[JAe] = {2, n_qp, n_sf_vec};
//     (*solution_matrix_idx_size_map)[JVe] = {3, 3 * n_qp, n_sf_vec};

//     return solution_matrix_idx_size_map;
// }
// std::unique_ptr<std::vector<std::vector<unsigned int>>> PLANE13Type::get_vector_idx_size_map()
// {
//     unsigned int n_sf = fe_array[0]->n_shape_functions();
//     std::unique_ptr<std::vector<std::vector<unsigned int>>> vector_idx_size_map =
//         std::make_unique<std::vector<std::vector<unsigned int>>>(VectorNames::NUM_VECTORS, std::vector<unsigned int>(2, 0));

//     (*vector_idx_size_map)[0] = {0, active_dof_vars_size * n_sf};

//     return vector_idx_size_map;
// }

// PLANE233

// PLANE233Type::PLANE233Type(const std::vector<unsigned int> &keyopts)
//     : AV2DElementType(keyopts)
// {
// } // Default dimension for PLANE13 is 2

// void PLANE233Type::populate_matrix_idx_size_map()
// {
//     libmesh_error_msg("PLANE233Type::populate_matrix_idx_size_map is not implemented.");
// }

// void PLANE233Type::populate_derived_matrix_idx_size_map()
// {
//     libmesh_error_msg("PLANE233Type::populate_derived_matrix_idx_size_map is not implemented.");
// }

// void PLANE233Type::populate_vector_idx_size_map()
// {
//     libmesh_error_msg("PLANE233Type::populate_vector_idx_size_map is not implemented.");
// }

VariableToElementIdsMapType ElementTypeFactory::associate_variables_with_element_IDs(VariableCategory variable_category) const
{
    VariableToElementIdsMapType variableToElementIds;
    for (const auto &pair : *data_map)
    {
        subdomain_id_type elementId = pair.first;
        const auto &elementType = pair.second;

        auto active_vars = elementType->get_active_variables(variable_category);
        for (const auto &var : active_vars)
        {
            (variableToElementIds)[var].insert(elementId);
        }
    }
    return variableToElementIds;
}

std::unique_ptr<ElementType> ElementTypeFactory::create_element_type(ElementTypeEnum type, const std::vector<unsigned int> &keyopts)
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

void ElementTypeFactory::create_element_types(ElementTypesMapCSV *csv_map)
{
    data_map->clear();
    for (const auto &pair : *csv_map)
    {
        unsigned int id = pair.first;
        const std::vector<unsigned int> &keyopts = pair.second;
        // You may need to determine the type from keyopts or another source.
        // Here, we assume the first entry in keyopts is the type enum value.
        if (keyopts.empty())
            continue;
        ElementTypeEnum type_enum = static_cast<ElementTypeEnum>(keyopts[0]);
        (*data_map)[id] = create_element_type(type_enum, keyopts);
        // Optionally, call set_conductive after construction:
        (*data_map)[id]->set_conductive();
    }
}