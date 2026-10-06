/**
 * @file element_field.h
 * @brief A physical field of an element type: DOF unknown, nodal field or element field.
 */
#pragma once

#include <algorithm>
#include <string>
#include <vector>

#include <libmesh/enum_fe_family.h>
#include <libmesh/enum_order.h>

namespace libMesh
{
    class FEAbstract;
}

namespace pytec3d
{
    enum class FieldKind : unsigned char
    {
        Dof,     ///< Unknown solved for by the system.
        Nodal,   ///< Post-processed field, L2-projected onto nodal shape functions.
        Element, ///< Post-processed field, averaged over the quadrature points.
    };

    /** @brief One libMesh variable of a field: a scalar component, or a whole H(curl) field. */
    struct LibMeshVariable
    {
        std::string label;
        libMesh::FEFamily family = libMesh::LAGRANGE;
        libMesh::Order order = libMesh::FIRST;
    };

    struct Field
    {
        std::string name;
        FieldKind kind = FieldKind::Dof;
        unsigned int index = 0; ///< Position among the fields of its kind.

        bool physical_vector = false;                ///< Declared as a vector field (a single active component makes it a scalar).
        bool vector_output = false;                  ///< Evaluated with a vector FE: always every component of the dimension.
        std::vector<unsigned int> active_components; ///< 0=x, 1=y, 2=z; every component for a vector output.

        libMesh::FEFamily family = libMesh::LAGRANGE;
        libMesh::Order order = libMesh::FIRST;
        std::vector<unsigned int> domains; ///< Physical regions where the field exists; empty means all.
        std::string nodal_loads;           ///< VTK point array that loads this DOF field, if any.

        libMesh::FEAbstract *fe = nullptr; ///< Shared by all fields of the same family and order; null for element fields.
        std::vector<LibMeshVariable> libmesh_variables;

        bool active_in(unsigned int domain) const
        {
            return domains.empty() || std::find(domains.begin(), domains.end(), domain) != domains.end();
        }
    };
} // namespace pytec3d
