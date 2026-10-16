/**
 * @file variables_map.h
 * @brief The libMesh variables of the systems: labels, grid array names and components.
 *
 * Every label a field can have as a libMesh variable is generated from fields.json and element_types.json
 * (registry::VariableId): a scalar field F is the variable F; a vector field F is the vector variable F
 * (LAGRANGE_VEC, L2_LAGRANGE_VEC or NEDELEC_ONE) with the components of the element dimension, a single component
 * F_z (out-of-plane, 2D), or, for element fields, the constant scalar variables F_x, F_y, F_z.
 * The grid arrays are named by the api label: F for the solution, _F for the Dirichlet values.
 *
 * An H(curl) (NEDELEC_ONE) variable F is told by its family. Its DOFs live on the edges, one per edge node
 * (the midside nodes of the grid): its Dirichlet values are component 2 of _F on those nodes
 * (edge_constraint_component).
 */
#ifndef VARIABLES_MAP_H
#define VARIABLES_MAP_H

#include <functional>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <libmesh/enum_fe_family.h>
#include <libmesh/enum_order.h>
#include <libmesh/id_types.h>

#include "pytec3d_system.h"
#include "registries.h"

namespace libMesh
{
    class System;
}

namespace pytec3d
{
    class ElementType;
}

class ElementTypesFactory;

using GlobalVariableName = pytec3d::registry::VariableId;

using VariableToElementIdsMapType = std::map<GlobalVariableName, std::set<libMesh::subdomain_id_type>>;

struct VariableInfo
{
    GlobalVariableName global_name;
    std::string global_label;     ///< Of the libMesh variable.
    std::string api_label;        ///< Of the grid arrays.
    unsigned int component_index; ///< Physical component of a component variable; 0 otherwise.
    bool vector;                  ///< A vector variable: one DOF per component on each node.
};

class VariablesMap
{
public:
    /** @brief One physical component of a vector field in a system. */
    struct VectorComponent
    {
        unsigned int var_number = libMesh::invalid_uint;
        unsigned int component_index = libMesh::invalid_uint; ///< Physical component, 0=x.
        unsigned int dof_component = 0;                        ///< Component of the DOF on a node: Node::dof_number(sys, var, dof_component).
        std::string api_label;
        std::string global_label;
    };

    static const VariableInfo &get_info_by_key(GlobalVariableName key);
    static std::vector<GlobalVariableName> get_keys_by_predicate(
        const std::function<bool(const VariableInfo &)> &pred);
    static GlobalVariableName get_key_by_predicate_one(
        const std::function<bool(const VariableInfo &)> &pred);

    /** @brief True when variable @p var of @p system is an H(curl) (NEDELEC_ONE) variable. */
    static bool is_edge_variable(const libMesh::System &system, unsigned int var);

    /** @brief True when @p system has an H(curl) variable. */
    static bool has_edge_variable(const libMesh::System &system);

    /** @brief Component of _F that holds the Dirichlet value of the edge DOFs of an H(curl) variable F. */
    static constexpr unsigned int edge_constraint_component = 2;

    /** @brief Point array with the Dirichlet values of the scalar or edge variable @p var of @p system. */
    static std::string constraint_array_name(const libMesh::System &system, unsigned int var);

    /** @brief Labels of the variables of @p system that @p label names: the variable itself, or the variables of a field. */
    static std::vector<std::string> get_global_labels_for_label(const libMesh::System &system, const std::string &label);

    /** @brief One libMesh variable of a system, as resolved from an element type. */
    struct SystemVariable
    {
        GlobalVariableName key;
        libMesh::FEFamily family;
        libMesh::Order order;
    };

    /** @brief libMesh variables the fields of the element type active in @p domain contribute to @p system_type. */
    static std::vector<SystemVariable> get_system_variables_for_element(const pytec3d::ElementType &element,
                                                                        pytec3d::registry::Domain domain,
                                                                        pyTEC3DSystemBase::SystemType system_type);

    /** @brief For every variable of @p system_type, the subdomains whose element type contains it. */
    static VariableToElementIdsMapType associate_variables_with_element_IDs(const ElementTypesFactory &factory,
                                                                            pyTEC3DSystemBase::SystemType system_type);

    /** @brief The nodal vector fields of @p system with their components, by api label (vector and per-component
     * variables; H(curl) variables are not nodal vectors). */
    static std::map<std::string, std::vector<VectorComponent>> get_vector_variables_for_system(const libMesh::System &system);

    /** @brief Variable label -> point array of its nodal loads, for the variables of @p system that have one. */
    static std::map<std::string, std::string> get_dof_variable_nodal_load_map_for_system(const libMesh::System &system);
};

#endif // VARIABLES_MAP_H
