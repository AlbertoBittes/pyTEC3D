#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <string>
#include <memory>
#include <vector>
#include <limits> // for std::numeric_limits
#include "element_types.h"
#include <libmesh/inf_fe.h>

using namespace Globals;

void ElementType::set_active_dof_vars_size()
{
    active_dof_vars_size = 0;
    for (const auto &pair : variable_info)
    {
        const auto &var_data = pair.second;
        if (var_data.is_active && var_data.category == VariableCategory::DOFVariables)
            ++active_dof_vars_size;
    }
}

AVElementType::AVElementType(const std::vector<unsigned int> &keyopts)
    : ElementType(keyopts) // Call the parent class constructor
{
    set_variable_info();
    set_fe_array();
}

void AVElementType::set_variable_info()
{
    // Additional initialization for AVElementType (if needed)
    variable_info[GlobalVariableName::B_x] = VariableData{VariableCategory::NodalSolutionVariables, true, APIVariableName::B, L2_LAGRANGE, order};
    variable_info[GlobalVariableName::B_y] = VariableData{VariableCategory::NodalSolutionVariables, true, APIVariableName::B, L2_LAGRANGE, order};
    variable_info[GlobalVariableName::J_z] = VariableData{VariableCategory::ElementSolutionVariables, true, APIVariableName::J, MONOMIAL, CONSTANT};
    variable_info[GlobalVariableName::Js_z] = VariableData{VariableCategory::LoadsVariables, true, APIVariableName::Js, MONOMIAL, CONSTANT};
}

void AVElementType::set_fe_array() {};

AV2DElementType::AV2DElementType(const std::vector<unsigned int> &keyopts) : AVElementType(keyopts) // Call the parent class constructor
{
    // Additional initialization for AVElementType (if needed)
    set_variable_info();
    set_fe_array();
}

void AV2DElementType::set_variable_info()
{
    // Additional initialization for AVElementType (if needed)
    variable_info[GlobalVariableName::VOLT] = VariableData{VariableCategory::DOFVariables, false, APIVariableName::VOLT, LAGRANGE, order};
}

void AV2DElementType::set_fe_array()
{
    FEType fe_type(order, LAGRANGE);
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

std::unique_ptr<std::vector<std::vector<unsigned int>>> AV2DElementType::get_dof_matrix_idx_size_map()
{
    // Check this when you have periodic BCs
    unsigned int n_sf = fe_array[0]->n_shape_functions();
    std::unique_ptr<std::vector<std::vector<unsigned int>>> matrix_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(2, std::vector<unsigned int>(2, 0));

    (*matrix_idx_size_map)[0] = {n_sf, n_sf};
    (*matrix_idx_size_map)[1] = {n_sf, n_sf};

    return matrix_idx_size_map;
}

std::unique_ptr<std::vector<std::vector<unsigned int>>> AV2DElementType::get_nodal_solution_matrix_idx_size_map()
{
    // Check this when you have periodic BCs
    unsigned int n_sf = fe_array[0]->n_shape_functions();
    unsigned int n_qp = fe_array[0]->n_quadrature_points();
    std::unique_ptr<std::vector<std::vector<unsigned int>>> solution_matrix_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(2, std::vector<unsigned int>(2, 0));

    (*solution_matrix_idx_size_map)[0] = {2 * n_sf, 2 * n_sf};
    (*solution_matrix_idx_size_map)[1] = {2 * n_sf, n_sf};

    return solution_matrix_idx_size_map;
}

std::unique_ptr<std::vector<std::vector<unsigned int>>> AV2DElementType::get_element_solution_matrix_idx_size_map()
{
    // Check this when you have periodic BCs
    unsigned int n_sf = fe_array[0]->n_shape_functions();
    unsigned int n_qp = fe_array[0]->n_quadrature_points();
    std::unique_ptr<std::vector<std::vector<unsigned int>>> solution_matrix_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(1, std::vector<unsigned int>(2, 0));

    (*solution_matrix_idx_size_map)[0] = {n_qp, n_sf};

    return solution_matrix_idx_size_map;
}

std::unique_ptr<std::vector<std::vector<unsigned int>>> AV2DElementType::get_vector_idx_size_map()
{
    unsigned int n_sf = fe_array[0]->n_shape_functions();
    std::unique_ptr<std::vector<std::vector<unsigned int>>> vector_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(1, std::vector<unsigned int>(1, 0));

    (*vector_idx_size_map)[0] = {active_dof_vars_size * n_sf};

    return vector_idx_size_map;
}

std::unique_ptr<std::vector<std::vector<unsigned int>>> AV2DElementType::get_nodal_solution_vector_idx_size_map()
{
    unsigned int n_sf = fe_array[0]->n_shape_functions();
    std::unique_ptr<std::vector<std::vector<unsigned int>>> vector_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(1, std::vector<unsigned int>(1, 0));

    (*vector_idx_size_map)[0] = {2 * n_sf};

    return vector_idx_size_map;
}

std::unique_ptr<std::vector<std::vector<unsigned int>>> AV2DElementType::get_element_solution_vector_idx_size_map()
{
    std::unique_ptr<std::vector<std::vector<unsigned int>>> vector_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(1, std::vector<unsigned int>(1, 0));

    (*vector_idx_size_map)[0] = {1};

    return vector_idx_size_map;
}

AV3DElementType::AV3DElementType(const std::vector<unsigned int> &keyopts) : AVElementType(keyopts) // Call the parent class constructor
{
    set_variable_info();
    set_fe_array();
}

void AV3DElementType::set_variable_info()
{
    // Additional initialization for AVElementType (if needed)
    variable_info[GlobalVariableName::B_z] = VariableData{VariableCategory::NodalSolutionVariables, true, APIVariableName::B, L2_LAGRANGE, order};
    variable_info[GlobalVariableName::J_x] = VariableData{VariableCategory::ElementSolutionVariables, true, APIVariableName::J, MONOMIAL, CONSTANT};
    variable_info[GlobalVariableName::J_y] = VariableData{VariableCategory::ElementSolutionVariables, true, APIVariableName::J, MONOMIAL, CONSTANT};
    variable_info[GlobalVariableName::VOLT] = VariableData{VariableCategory::DOFVariables, false, APIVariableName::VOLT, LAGRANGE, order};
    variable_info[GlobalVariableName::EMF] = VariableData{VariableCategory::DOFVariables, false, APIVariableName::EMF, LAGRANGE, order};
    variable_info[GlobalVariableName::Js_x] = VariableData{VariableCategory::LoadsVariables, true, APIVariableName::Js, MONOMIAL, CONSTANT};
    variable_info[GlobalVariableName::Js_y] = VariableData{VariableCategory::LoadsVariables, true, APIVariableName::Js, MONOMIAL, CONSTANT};
}

void AV3DElementType::set_fe_array()
{
    FEType fe_type(order, LAGRANGE);
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

    (*matrix_idx_size_map)[0] = {active_dof_vars_size * n_sf, active_dof_vars_size * n_sf};
    (*matrix_idx_size_map)[1] = {active_dof_vars_size * n_sf, active_dof_vars_size * n_sf};

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

    (*solution_matrix_idx_size_map)[0] = {3 * n_qp, n_sf_vec};
    (*solution_matrix_idx_size_map)[1] = {3 * n_qp, n_sf};

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

std::unique_ptr<std::vector<std::vector<unsigned int>>> AV3DElementType::get_nodal_solution_vector_idx_size_map()
{
    unsigned int n_sf = fe_array[0]->n_shape_functions();
    std::unique_ptr<std::vector<std::vector<unsigned int>>> vector_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(1, std::vector<unsigned int>(1, 0));

    (*vector_idx_size_map)[0] = {3 * n_sf};

    return vector_idx_size_map;
}

std::unique_ptr<std::vector<std::vector<unsigned int>>> AV3DElementType::get_element_solution_vector_idx_size_map()
{
    std::unique_ptr<std::vector<std::vector<unsigned int>>> vector_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(1, std::vector<unsigned int>(1, 0));

    (*vector_idx_size_map)[0] = {3};

    return vector_idx_size_map;
}

// SOLID97

SOLID97Type::SOLID97Type(const std::vector<unsigned int> &keyopts)
    : AV3DElementType(keyopts)
{
    set_variable_info();
    set_fe_array();
    set_active_dof_vars_size();
}

void SOLID97Type::set_variable_info()
{
    variable_info[GlobalVariableName::A_x] = VariableData{VariableCategory::DOFVariables, true, APIVariableName::A, LAGRANGE, order};
    variable_info[GlobalVariableName::A_y] = VariableData{VariableCategory::DOFVariables, true, APIVariableName::A, LAGRANGE, order};
    variable_info[GlobalVariableName::A_z] = VariableData{VariableCategory::DOFVariables, true, APIVariableName::A, LAGRANGE, order};
    variable_info[GlobalVariableName::CURR] = VariableData{VariableCategory::DOFVariables, false, APIVariableName::CURR, LAGRANGE, order};

    if (physical_region == SOLID_CONDUCTOR || physical_region == SOLID_CONDUCTOR_CIRCUIT_COUPLED)
        variable_info[GlobalVariableName::VOLT].is_active = true; // Modify the first member
    if (physical_region == STRANDED_COIL_CIRCUIT_COUPLED || physical_region == STRANDED_COIL_VOLTAGE_FED)
        variable_info[GlobalVariableName::CURR].is_active = true; // Modify the first member
    if (physical_region == STRANDED_COIL_CIRCUIT_COUPLED || physical_region == SOLID_CONDUCTOR_CIRCUIT_COUPLED)
        variable_info[GlobalVariableName::EMF].is_active = true; // Modify the first member
}

void SOLID97Type::set_fe_array()
{
    FEType fe_type_vec(order, LAGRANGE_VEC);
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

    FEType fe_type_face(order, LAGRANGE);

    std::unique_ptr<FEBase> fe_face = FEBase::build(3, fe_type_face);

    // Boundary integration requires another quadrature rule,
    // with dimensionality one less than the dimensionality
    // of the element.
    // In 1D, the Clough and Gauss quadrature rules are identical.
    QGauss *qface = new QGauss(2, fe_type_face.default_quadrature_order());

    // Tell the finite element object to use our
    // quadrature rule.
    fe_face->attach_quadrature_rule(qface);

    fe_array.push_back(std::move(fe_face));
}

// SOLID 236

SOLID236Type::SOLID236Type(const std::vector<unsigned int> &keyopts)
    : AV3DElementType(keyopts)
{
    // variable_info["DOFVariables"]["A"] =
    //     std::make_tuple(true, 3, false, "A", 2, "NEDELEC_ONE", "FIRST");
    // if (physical_region == 1 || physical_region == 2)
    //     std::get<0>(variable_info["DOFVariables"]["VOLT"]) = true; // Modify the first member
    // if (physical_region == 2)
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
    libmesh_error_msg("PLANE233Type::PLANE233Type is not implemented.");
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
    set_variable_info();
    set_fe_array();
    set_active_dof_vars_size();
} // Default dimension for PLANE13 is 2

void PLANE13Type::set_variable_info()
{
    variable_info[GlobalVariableName::A_z] = {VariableCategory::DOFVariables, true, APIVariableName::A, LAGRANGE, order};

    if (physical_region == STRANDED_COIL)
        variable_info[GlobalVariableName::VOLT].is_active = true; // Modify the first member
    if (physical_region == STRANDED_COIL || physical_region == SOLID_CONDUCTOR_CIRCUIT_COUPLED || physical_region == STRANDED_COIL_CIRCUIT_COUPLED)
        variable_info[GlobalVariableName::EMF].is_active = true; // Modify the first member
    if (physical_region == STRANDED_COIL_CIRCUIT_COUPLED || physical_region == STRANDED_COIL_VOLTAGE_FED || physical_region == SOLID_CONDUCTOR_CIRCUIT_COUPLED)
        variable_info[GlobalVariableName::CURR].is_active = true;
}

void PLANE13Type::set_fe_array() {};

INFIN110Type::INFIN110Type(const std::vector<unsigned int> &keyopts)
    : PLANE13Type(keyopts)
{
    set_variable_info();
    set_fe_array();
    set_active_dof_vars_size();
}

void INFIN110Type::set_variable_info() {};

void INFIN110Type::set_fe_array()
{
    // Build a Finite Element object of the specified type.  Since the
    // FEBase::build() member dynamically creates memory we will
    // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.
    // InfMapType inf_map_type = static_cast<InfMapType>(keyopts[5]);
    FEType fe_type(order, LAGRANGE);
    // Build a Finite Element object of the specified type.
    std::unique_ptr<FEBase> fe = FEBase::build_InfFE(2, fe_type);
    // A 3rd order Gauss quadrature rule for numerical integration.
    QGauss *qrule_vec = new QGauss(2, fe_type.default_quadrature_order());
    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule(qrule_vec);

    // Replace the first element of fe_array with the new pointer
    if (fe_array.empty())
        fe_array.push_back(std::move(fe));
    else
        fe_array[0] = std::move(fe);
};

INFIN111Type::INFIN111Type(const std::vector<unsigned int> &keyopts)
    : SOLID97Type(keyopts)
{
    set_variable_info();
    set_fe_array();
    set_active_dof_vars_size();
}

void INFIN111Type::set_variable_info() {};

void INFIN111Type::set_fe_array()
{
    FEType fe_type(order, LAGRANGE);
    // Build a Finite Element object of the specified type.
    std::unique_ptr<FEBase> fe = FEBase::build_InfFE(3, fe_type);
    // A 3rd order Gauss quadrature rule for numerical integration.
    QGauss *qrule_vec = new QGauss(3, fe_type.default_quadrature_order());
    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule(qrule_vec);

    // Replace the first element of fe_array with the new pointer
    if (fe_array.empty())
        fe_array.push_back(std::move(fe));
    else
        fe_array[0] = std::move(fe);
};
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
