#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <string>
#include <memory>
#include <vector>
#include <limits> // for std::numeric_limits
#include "element_types.h"
#include "element_builder.h"
#include <libmesh/inf_fe.h>

using namespace Globals;

void ElementType::set_variables_order()
{
    if (order != FIRST)
        for (auto &pair : variable_info)
        {
            if (pair.second.fe_family == LAGRANGE || pair.second.fe_family == L2_LAGRANGE)
            {
                pair.second.order = order;
            }
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
    variable_info[GlobalVariableName::B_x] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::B_x);
    variable_info[GlobalVariableName::B_y] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::B_y);
    variable_info[GlobalVariableName::H_x] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::H_x);
    variable_info[GlobalVariableName::H_y] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::H_y);
    variable_info[GlobalVariableName::Jt_z] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::Jt_z);
    variable_info[GlobalVariableName::Be_x] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::Be_x);
    variable_info[GlobalVariableName::Be_y] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::Be_y);
    variable_info[GlobalVariableName::He_x] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::He_x);
    variable_info[GlobalVariableName::He_y] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::He_y);
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
    variable_info[GlobalVariableName::CSG_z] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::CSG_z);
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

std::unique_ptr<std::vector<std::vector<unsigned int>>> AV2DElementType::get_solution_matrix_idx_size_map()
{
    // Check this when you have periodic BCs
    unsigned int n_sf = fe_array[0]->n_shape_functions();
    unsigned int n_qp = fe_array[0]->n_quadrature_points();
    std::unique_ptr<std::vector<std::vector<unsigned int>>> solution_matrix_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(5, std::vector<unsigned int>(2, 0));

    // Me
    (*solution_matrix_idx_size_map)[ElementBuilder::PLANE13SolutionMatrices::Me_PLANE13] = {2 * n_sf, 2 * n_sf};
    // AAe
    (*solution_matrix_idx_size_map)[ElementBuilder::PLANE13SolutionMatrices::AAe_PLANE13] = {2 * n_sf, n_sf};
    // JAe
    (*solution_matrix_idx_size_map)[ElementBuilder::PLANE13SolutionMatrices::JAe_PLANE13] = {n_qp, n_sf};
    // JHe
    // (*solution_matrix_idx_size_map)[3] = {n_qp, 2 * n_sf};
    // BAe
    (*solution_matrix_idx_size_map)[ElementBuilder::PLANE13SolutionMatrices::BAe_PLANE13] = {2 * n_qp, n_sf};

    return solution_matrix_idx_size_map;
}

AV3DElementType::AV3DElementType(const std::vector<unsigned int> &keyopts) : AVElementType(keyopts) // Call the parent class constructor
{
    set_variable_info();
    set_fe_array();
}

void AV3DElementType::set_variable_info()
{
    // Additional initialization for AVElementType (if needed)
    variable_info[GlobalVariableName::B_z] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::B_z);
    variable_info[GlobalVariableName::H_z] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::H_z);
    variable_info[GlobalVariableName::Jt_x] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::Jt_x);
    variable_info[GlobalVariableName::Jt_y] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::Jt_y);
    variable_info[GlobalVariableName::Be_z] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::Be_z);
    variable_info[GlobalVariableName::He_z] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::He_z);
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

// SOLID97

SOLID97Type::SOLID97Type(const std::vector<unsigned int> &keyopts)
    : AV3DElementType(keyopts)
{
    set_variable_info();
    set_fe_array();
    set_variables_order();
}

void SOLID97Type::set_variable_info()
{
    variable_info[GlobalVariableName::A_x] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::A_x);
    variable_info[GlobalVariableName::A_y] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::A_y);
    variable_info[GlobalVariableName::A_z] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::A_z);

    variable_info[GlobalVariableName::CSG_x] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::CSG_x);
    variable_info[GlobalVariableName::CSG_y] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::CSG_y);
    variable_info[GlobalVariableName::CSG_z] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::CSG_z);

    if (physical_region == SOLID_CONDUCTOR || physical_region == SOLID_CONDUCTOR_CIRCUIT_COUPLED)
    {
        variable_info[GlobalVariableName::VOLT] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::VOLT);
        variable_info[GlobalVariableName::AMPS] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::AMPS);
    }
    if (physical_region == STRANDED_COIL_CIRCUIT_COUPLED || physical_region == STRANDED_COIL_VOLTAGE_FED)
    {
        variable_info[GlobalVariableName::CURR] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::CURR);
        variable_info[GlobalVariableName::VLTG] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::VLTG);
    }
    if (physical_region == STRANDED_COIL_CIRCUIT_COUPLED || physical_region == SOLID_CONDUCTOR_CIRCUIT_COUPLED)
    {
        variable_info[GlobalVariableName::EMF] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::EMF);
        variable_info[GlobalVariableName::CURT] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::CURT);
    }
    if (physical_region == STRANDED_CONDUCTOR)
    {
        variable_info[GlobalVariableName::Js_x] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::Js_x);
        variable_info[GlobalVariableName::Js_y] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::Js_y);
        variable_info[GlobalVariableName::Js_z] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::Js_z);
    }
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

std::unique_ptr<std::vector<std::vector<unsigned int>>> SOLID97Type::get_solution_matrix_idx_size_map()
{
    // Check this when you have periodic BCs
    unsigned int n_sf = fe_array[0]->n_shape_functions();
    unsigned int n_sf_vec = fe_array[1]->n_shape_functions();
    unsigned int n_qp = fe_array[0]->n_quadrature_points();
    std::unique_ptr<std::vector<std::vector<unsigned int>>> solution_matrix_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(5, std::vector<unsigned int>(2, 0));

    // Me
    (*solution_matrix_idx_size_map)[ElementBuilder::SOLID97SolutionMatrices::Me_SOLID97] = {n_sf_vec, n_sf_vec};
    // AAe
    (*solution_matrix_idx_size_map)[ElementBuilder::SOLID97SolutionMatrices::AAe_SOLID97] = {n_sf_vec, n_sf_vec};
    // JAe
    (*solution_matrix_idx_size_map)[ElementBuilder::SOLID97SolutionMatrices::JAe_SOLID97] = {3 * n_qp, n_sf_vec};
    // JVe
    (*solution_matrix_idx_size_map)[ElementBuilder::SOLID97SolutionMatrices::JVe_SOLID97] = {3 * n_qp, n_sf};
    // BAe
    (*solution_matrix_idx_size_map)[ElementBuilder::SOLID97SolutionMatrices::BAe_SOLID97] = {3 * n_qp, n_sf_vec};

    return solution_matrix_idx_size_map;
}

// SOLID 236

SOLID236Type::SOLID236Type(const std::vector<unsigned int> &keyopts)
    : AV3DElementType(keyopts)
{
    set_variable_info();
    set_fe_array();
    set_variables_order();
}

void SOLID236Type::set_variable_info()
{
    variable_info[GlobalVariableName::Az] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::Az);
    variable_info[GlobalVariableName::Az].fe_family = NEDELEC_ONE;
    // variable_info[GlobalVariableName::CSG_x] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::CSG_x);
    // variable_info[GlobalVariableName::CSG_y] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::CSG_y);
    // variable_info[GlobalVariableName::CSG_z] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::CSG_z);

    if (physical_region == SOLID_CONDUCTOR || physical_region == SOLID_CONDUCTOR_CIRCUIT_COUPLED)
    {
        variable_info[GlobalVariableName::VOLT] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::VOLT);
        variable_info[GlobalVariableName::AMPS] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::AMPS);
    }
    if (physical_region == STRANDED_COIL)
    {
        variable_info[GlobalVariableName::EMF] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::EMF);
        variable_info[GlobalVariableName::CURT] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::CURT);
    }
    if (physical_region == STRANDED_CONDUCTOR)
    {
        variable_info[GlobalVariableName::Js_x] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::Js_x);
        variable_info[GlobalVariableName::Js_y] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::Js_y);
        variable_info[GlobalVariableName::Js_z] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::Js_z);
    }

    variable_info[GlobalVariableName::TREE] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::TREE);
}

void SOLID236Type::set_fe_array()
{
    FEType fe_type_vec(FIRST, NEDELEC_ONE);
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

    std::unique_ptr<FEVectorBase> fe_vec_copy = FEVectorBase::build(3, fe_type_vec);

    QGauss *qrule_vec_copy = new QGauss(3, fe_type_vec.default_quadrature_order());

    fe_vec_copy->attach_quadrature_rule(qrule_vec_copy);

    fe_array.push_back(std::move(fe_vec_copy));

    FEType fe_type_vec_lagrange(order, LAGRANGE_VEC);

    std::unique_ptr<FEVectorBase> fe_vec_a_lagrange = FEVectorBase::build(3, fe_type_vec_lagrange);

    QGauss *qrule_vec_a_lagrange = new QGauss(3, fe_type_vec_lagrange.default_quadrature_order());

    fe_vec_a_lagrange->attach_quadrature_rule(qrule_vec_a_lagrange);

    fe_array.push_back(std::move(fe_vec_a_lagrange));
}

std::unique_ptr<std::vector<std::vector<unsigned int>>> SOLID236Type::get_solution_matrix_idx_size_map()
{
    // TEMPORARY !!!!!!!!!!!!!!!!!!!!!!!!     /////
    // Check this when you have periodic BCs
    unsigned int n_sf = fe_array[0]->n_shape_functions();
    unsigned int n_sf_vec = 3 * n_sf;
    unsigned int n_sf_edge = fe_array[1]->n_shape_functions();
    unsigned int n_qp = fe_array[0]->n_quadrature_points();
    std::unique_ptr<std::vector<std::vector<unsigned int>>> solution_matrix_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(5, std::vector<unsigned int>(2, 0));

    // Me
    (*solution_matrix_idx_size_map)[ElementBuilder::SOLID97SolutionMatrices::Me_SOLID97] = {n_sf_vec, n_sf_vec};
    // AAe
    (*solution_matrix_idx_size_map)[ElementBuilder::SOLID97SolutionMatrices::AAe_SOLID97] = {n_sf_vec, n_sf_edge};
    // JAe
    (*solution_matrix_idx_size_map)[ElementBuilder::SOLID97SolutionMatrices::JAe_SOLID97] = {3 * n_qp, n_sf_edge};
    // JVe
    (*solution_matrix_idx_size_map)[ElementBuilder::SOLID97SolutionMatrices::JVe_SOLID97] = {3 * n_qp, n_sf};
    // BAe
    (*solution_matrix_idx_size_map)[ElementBuilder::SOLID97SolutionMatrices::BAe_SOLID97] = {3 * n_qp, n_sf_edge};

    return solution_matrix_idx_size_map;
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
    set_variables_order();
} // Default dimension for PLANE13 is 2

void PLANE13Type::set_variable_info()
{
    variable_info[GlobalVariableName::A_z] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::A_z);
    variable_info[GlobalVariableName::CSG_z] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::CSG_z);

    if (physical_region == STRANDED_COIL)
    {
        variable_info[GlobalVariableName::VOLT] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::VOLT);
        variable_info[GlobalVariableName::AMPS] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::AMPS);
    }
    if (physical_region == STRANDED_COIL || physical_region == SOLID_CONDUCTOR_CIRCUIT_COUPLED || physical_region == STRANDED_COIL_CIRCUIT_COUPLED)
    {
        variable_info[GlobalVariableName::EMF] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::EMF);
        variable_info[GlobalVariableName::CURT] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::CURT);
    }
    if (physical_region == STRANDED_COIL_CIRCUIT_COUPLED || physical_region == STRANDED_COIL_VOLTAGE_FED || physical_region == SOLID_CONDUCTOR_CIRCUIT_COUPLED)
    {
        variable_info[GlobalVariableName::CURR] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::CURR);
        variable_info[GlobalVariableName::VLTG] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::VLTG);
    }
    if (physical_region == STRANDED_CONDUCTOR)
    {
        variable_info[GlobalVariableName::Js_z] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, GlobalVariableName::Js_z);
    }
}

void PLANE13Type::set_fe_array() {};

INFIN110Type::INFIN110Type(const std::vector<unsigned int> &keyopts)
    : PLANE13Type(keyopts)
{
    set_variable_info();
    set_fe_array();
    set_variables_order();
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
    set_variables_order();
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
