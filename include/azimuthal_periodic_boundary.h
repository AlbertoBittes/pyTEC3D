// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

// <h1> Systems Example 9 - Linear elasticity problem with periodic constraints </h1>
//
// In this example we illustrate periodic constraints with a linear elasticity example.
// We consider a sector of a circular domain, and hence we must impose an "azimuthal"
// periodic condition, which entails a dof-transformation between the two periodic
// boundaries to account for the rotation of the coordinate system.
//
// The baseline code is from systems_of_equations_ex6, and we edit to use a different mesh
// and to impose the periodic boundary condition.

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>

// libMesh includes
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/perf_log.h"
#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"
#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/solver_configuration.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/periodic_boundaries.h"
#include "libmesh/periodic_boundary.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"

#ifndef AZIMUTHAL_PERIODIC_BOUNDARY_H
#define AZIMUTHAL_PERIODIC_BOUNDARY_H

// Bring in everything from the libMesh namespace
using namespace libMesh;

#ifdef LIBMESH_ENABLE_PERIODIC
// Here we define the azimuthal periodic boundary condition class.
// This class assumes that (u,v,w) are variables 0, 1, and 2 in
// the System, but this could be made more general if needed.
class AzimuthalPeriodicBoundary : public PeriodicBoundaryBase
{
public:
    AzimuthalPeriodicBoundary(
        Point center,
        Point axis,
        Real angle)
        : PeriodicBoundaryBase(),
          _center(center),
          _axis(axis),
          _theta(angle)
    {
        set_variable(0);
        // set_variable(1);
        // set_variable(2);
        set_up_rotation_matrix();
    }

    AzimuthalPeriodicBoundary(const AzimuthalPeriodicBoundary &o, TransformationType t = FORWARD)
        : PeriodicBoundaryBase(o),
          _center(o._center),
          _axis(o._axis),
          _theta(o._theta)
    {
        if (t == INVERSE)
        {
            std::swap(myboundary, pairedboundary);
            _theta *= -1.0;
        }

        set_variable(0);
        // set_variable(1);
        // set_variable(2);
        set_up_rotation_matrix();
    }

    virtual ~AzimuthalPeriodicBoundary() = default;

    void set_up_rotation_matrix();

    virtual Point get_corresponding_pos(const Point &pt) const override
    {
        // Note that since _theta defines the angle from "paired boundary" to
        // "my boundary", and we want the inverse of that here, we multiply
        // by the transpose of the transformation matrix below.
        return _R.left_multiply(pt - _center) + _center;
    }

    virtual std::unique_ptr<PeriodicBoundaryBase> clone(TransformationType t = FORWARD) const override
    {
        return std::make_unique<AzimuthalPeriodicBoundary>(*this, t);
    }

private:
    // Define the properties needed for the azimuthal periodic boundary.
    // Note that _theta specifies the angle from "paired boundary" to
    // "my boundary".
    Point _center;
    Point _axis;
    Real _theta;
    RealTensorValue _R;
};

#endif // LIBMESH_ENABLE_PERIODIC

#endif // AZIMUTHAL_PERIODIC_BOUNDARY_H