// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// Local includes
#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/enum_to_string.h"

namespace libMesh
{

  template <>
  RealGradient FE<3, NEDELEC_ONE>::shape(const Elem *elem,
                                         const Order order,
                                         const unsigned int i,
                                         const Point &p,
                                         const bool add_p_level)
  {
#if LIBMESH_DIM == 3
    libmesh_assert(elem);

    const Order totalorder = order + add_p_level * elem->p_level();
    libmesh_assert_less(i, n_dofs(elem->type(), totalorder));

    const char sign = elem->positive_edge_orientation(i) ? 1 : -1;

    const Real xi = p(0);
    const Real eta = p(1);
    const Real zeta = p(2);

    switch (totalorder)
    {
      // linear Nedelec (first kind) shape functions
    case FIRST:
    {
      switch (elem->type())
      {
      case HEX20:
      case HEX27:
      {
        // Even with a loose inverse_map tolerance we ought to
        // be nearly on the element interior in master
        // coordinates
        libmesh_assert_less_equal(std::fabs(xi), 1.0 + 10 * TOLERANCE);
        libmesh_assert_less_equal(std::fabs(eta), 1.0 + 10 * TOLERANCE);
        libmesh_assert_less_equal(std::fabs(zeta), 1.0 + 10 * TOLERANCE);

        switch (i)
        {
        case 0:
          return sign * RealGradient(-0.125 * (1.0 - eta - zeta + eta * zeta), 0.0, 0.0);
        case 1:
          return sign * RealGradient(0.0, -0.125 * (1.0 + xi - zeta - xi * zeta), 0.0);
        case 2:
          return sign * RealGradient(0.125 * (1.0 + eta - zeta - eta * zeta), 0.0, 0.0);
        case 3:
          return sign * RealGradient(0.0, -0.125 * (1.0 - xi - zeta + xi * zeta), 0.0);
        case 4:
          return sign * RealGradient(0.0, 0.0, -0.125 * (1.0 - xi - eta + xi * eta));
        case 5:
          return sign * RealGradient(0.0, 0.0, -0.125 * (1.0 + xi - eta - xi * eta));
        case 6:
          return sign * RealGradient(0.0, 0.0, -0.125 * (1.0 + xi + eta + xi * eta));
        case 7:
          return sign * RealGradient(0.0, 0.0, -0.125 * (1.0 - xi + eta - xi * eta));
        case 8:
          return sign * RealGradient(-0.125 * (1.0 - eta + zeta - eta * zeta), 0.0, 0.0);
        case 9:
          return sign * RealGradient(0.0, -0.125 * (1.0 + xi + zeta + xi * zeta), 0.0);
        case 10:
          return sign * RealGradient(0.125 * (1.0 + eta + zeta + eta * zeta), 0.0, 0.0);
        case 11:
          return sign * RealGradient(0.0, -0.125 * (1.0 - xi + zeta - xi * zeta), 0.0);

        default:
          libmesh_error_msg("Invalid i = " << i);
        }
      }

      case TET10:
      case TET14:
      {
        switch (i)
        {
        case 0:
          return sign * RealGradient(-1.0 + eta + zeta, -xi, -xi);
        case 1:
          return sign * RealGradient(eta, -xi, 0.0);
        case 2:
          return sign * RealGradient(-eta, -1.0 + xi + zeta, -eta);
        case 3:
          return sign * RealGradient(-zeta, -zeta, -1.0 + xi + eta);
        case 4:
          return sign * RealGradient(zeta, 0.0, -xi);
        case 5:
          return sign * RealGradient(0.0, zeta, -eta);

        default:
          libmesh_error_msg("Invalid i = " << i);
        }
      }
      case PRISM15:
      case PRISM18:
      {
        // Reference prism: triangle in (xi,eta) and zeta in [-1,1]
        libmesh_assert_greater_equal(xi, -10. * TOLERANCE);
        libmesh_assert_greater_equal(eta, -10. * TOLERANCE);
        libmesh_assert_less_equal(xi + eta, 1.0 + 10. * TOLERANCE);
        libmesh_assert_less_equal(std::fabs(zeta), 1.0 + 10. * TOLERANCE);

        const Real f0 = 0.5 * (1.0 - zeta); // bottom blend
        const Real f1 = 0.5 * (1.0 + zeta); // top blend

        switch (i)
        {
        // Bottom face edges: (0,1), (1,2), (0,2)
        case 0:
          return sign * RealGradient(f0 * (1.0 - eta), f0 * xi, 0.0);
        case 1:
          return sign * RealGradient(f0 * (-eta), f0 * xi, 0.0);
        case 2:
          return sign * RealGradient(f0 * (eta), f0 * (1.0 - xi), 0.0);

        // Vertical edges: (0,3), (1,4), (2,5)
        // 0.5 factor ensures unit tangential integral on zeta in [-1,1]
        case 3:
          return sign * RealGradient(0.0, 0.0, 1.0 - xi - eta);
        case 4:
          return sign * RealGradient(0.0, 0.0, xi);
        case 5:
          return sign * RealGradient(0.0, 0.0, eta);

        // Top face edges: (3,4), (4,5), (3,5)
        case 6:
          return sign * RealGradient(f1 * (1.0 - eta), f1 * xi, 0.0);
        case 7:
          return sign * RealGradient(f1 * (-eta), f1 * xi, 0.0);
        case 8:
          return sign * RealGradient(f1 * (eta), f1 * (1.0 - xi), 0.0);

        default:
          libmesh_error_msg("Invalid i = " << i);
        }
      }
      case PYRAMID13:
      case PYRAMID14:
      {
        libmesh_assert_less(i, 8);

        const Real xi = p(0);
        const Real eta = p(1);
        const Real zeta = p(2);

        // Reference pyramid: base zeta=0, apex zeta=1
        libmesh_assert_less_equal(std::fabs(xi), 1.0 + 10. * TOLERANCE);
        libmesh_assert_less_equal(std::fabs(eta), 1.0 + 10. * TOLERANCE);
        libmesh_assert_greater_equal(zeta, -10. * TOLERANCE);
        libmesh_assert_less_equal(zeta, 1.0 + 10. * TOLERANCE);

        Real a = 1.0 - zeta;
        if (std::fabs(a) < 1e-14)
          a = (a >= 0.0) ? 1e-14 : -1e-14;

        const Real inva = 1.0 / a;
        const Real inva2 = inva * inva;

        const Real l0 = 0.25 * (a - xi - eta + xi * eta * inva);
        const Real l1 = 0.25 * (a + xi - eta - xi * eta * inva);
        const Real l2 = 0.25 * (a + xi + eta + xi * eta * inva);
        const Real l3 = 0.25 * (a - xi + eta - xi * eta * inva);
        const Real l4 = zeta;

        const RealGradient g0(0.25 * (-1.0 + eta * inva), 0.25 * (-1.0 + xi * inva), 0.25 * (-1.0 + xi * eta * inva2));
        const RealGradient g1(0.25 * (1.0 - eta * inva), 0.25 * (-1.0 - xi * inva), 0.25 * (-1.0 - xi * eta * inva2));
        const RealGradient g2(0.25 * (1.0 + eta * inva), 0.25 * (1.0 + xi * inva), 0.25 * (-1.0 + xi * eta * inva2));
        const RealGradient g3(0.25 * (-1.0 - eta * inva), 0.25 * (1.0 - xi * inva), 0.25 * (-1.0 - xi * eta * inva2));
        const RealGradient g4(0.0, 0.0, 1.0);

        RealGradient v;
        switch (i)
        {
        // Base edges
        case 0:
          v = l0 * g1 - l1 * g0;
          break; // (0,1)
        case 1:
          v = l1 * g2 - l2 * g1;
          break; // (1,2)
        case 2:
          v = l2 * g3 - l3 * g2;
          break; // (2,3)
        case 3:
          v = l3 * g0 - l0 * g3;
          break; // (3,0)
        // Side edges to apex
        case 4:
          v = l0 * g4 - l4 * g0;
          break; // (0,4)
        case 5:
          v = l1 * g4 - l4 * g1;
          break; // (1,4)
        case 6:
          v = l2 * g4 - l4 * g2;
          break; // (2,4)
        case 7:
          v = l3 * g4 - l4 * g3;
          break; // (3,4)
        default:
          libmesh_error_msg("Invalid i = " << i);
        }

        if (!elem->positive_edge_orientation(i))
          v = -v;
        return v;
      }

      default:
        libmesh_error_msg("ERROR: Unsupported 3D element type!: " << Utility::enum_to_string(elem->type()));
      }
    }

      // unsupported order
    default:
      libmesh_error_msg("ERROR: Unsupported 3D FE order!: " << totalorder);
    }

#else // LIBMESH_DIM != 3
    libmesh_ignore(elem, order, i, p, add_p_level);
    libmesh_not_implemented();
#endif
  }

  template <>
  RealGradient FE<3, NEDELEC_ONE>::shape(const ElemType,
                                         const Order,
                                         const unsigned int,
                                         const Point &)
  {
    libmesh_error_msg("Nedelec elements require the element type \nbecause edge orientation is needed.");
    return RealGradient();
  }

  template <>
  RealGradient FE<3, NEDELEC_ONE>::shape(const FEType fet,
                                         const Elem *elem,
                                         const unsigned int i,
                                         const Point &p,
                                         const bool add_p_level)
  {
    return FE<3, NEDELEC_ONE>::shape(elem, fet.order, i, p, add_p_level);
  }

  template <>
  RealGradient FE<3, NEDELEC_ONE>::shape_deriv(const Elem *elem,
                                               const Order order,
                                               const unsigned int i,
                                               const unsigned int j,
                                               const Point &p,
                                               const bool add_p_level)
  {
#if LIBMESH_DIM == 3
    libmesh_assert(elem);
    libmesh_assert_less(j, 3);

    const Order totalorder = order + add_p_level * elem->p_level();
    libmesh_assert_less(i, n_dofs(elem->type(), totalorder));

    const char sign = elem->positive_edge_orientation(i) ? 1 : -1;

    const Real xi = p(0);
    const Real eta = p(1);
    const Real zeta = p(2);

    switch (totalorder)
    {
      // linear Nedelec (first kind) shape function first derivatives
    case FIRST:
    {
      switch (elem->type())
      {
      case HEX20:
      case HEX27:
      {
        // Even with a loose inverse_map tolerance we ought to
        // be nearly on the element interior in master
        // coordinates
        libmesh_assert_less_equal(std::fabs(xi), 1.0 + 10 * TOLERANCE);
        libmesh_assert_less_equal(std::fabs(eta), 1.0 + 10 * TOLERANCE);
        libmesh_assert_less_equal(std::fabs(zeta), 1.0 + 10 * TOLERANCE);

        switch (j)
        {
          // d()/dxi
        case 0:
        {
          switch (i)
          {
          case 0:
          case 2:
          case 8:
          case 10:
            return RealGradient();
          case 1:
            return sign * RealGradient(0.0, -0.125 * (1.0 - zeta));
          case 3:
            return sign * RealGradient(0.0, -0.125 * (-1.0 + zeta));
          case 4:
            return sign * RealGradient(0.0, 0.0, -0.125 * (-1.0 + eta));
          case 5:
            return sign * RealGradient(0.0, 0.0, -0.125 * (1.0 - eta));
          case 6:
            return sign * RealGradient(0.0, 0.0, -0.125 * (1.0 + eta));
          case 7:
            return sign * RealGradient(0.0, 0.0, -0.125 * (-1.0 - eta));
          case 9:
            return sign * RealGradient(0.0, -0.125 * (1.0 + zeta), 0.0);
          case 11:
            return sign * RealGradient(0.0, -0.125 * (-1.0 - zeta), 0.0);

          default:
            libmesh_error_msg("Invalid i = " << i);
          } // switch(i)

        } // j = 0

          // d()/deta
        case 1:
        {
          switch (i)
          {
          case 1:
          case 3:
          case 9:
          case 11:
            return RealGradient();
          case 0:
            return sign * RealGradient(-0.125 * (-1.0 + zeta), 0.0, 0.0);
          case 2:
            return sign * RealGradient(0.125 * (1.0 - zeta), 0.0, 0.0);
          case 4:
            return sign * RealGradient(0.0, 0.0, -0.125 * (-1.0 + xi));
          case 5:
            return sign * RealGradient(0.0, 0.0, -0.125 * (-1.0 - xi));
          case 6:
            return sign * RealGradient(0.0, 0.0, -0.125 * (1.0 + xi));
          case 7:
            return sign * RealGradient(0.0, 0.0, -0.125 * (1.0 - xi));
          case 8:
            return sign * RealGradient(-0.125 * (-1.0 - zeta), 0.0, 0.0);
          case 10:
            return sign * RealGradient(0.125 * (1.0 + zeta), 0.0, 0.0);

          default:
            libmesh_error_msg("Invalid i = " << i);
          } // switch(i)

        } // j = 1

          // d()/dzeta
        case 2:
        {
          switch (i)
          {
          case 4:
          case 5:
          case 6:
          case 7:
            return RealGradient();
          case 0:
            return sign * RealGradient(-0.125 * (-1.0 + eta), 0.0, 0.0);
          case 1:
            return sign * RealGradient(0.0, -0.125 * (-1.0 - xi), 0.0);
          case 2:
            return sign * RealGradient(0.125 * (-1.0 - eta), 0.0, 0.0);
          case 3:
            return sign * RealGradient(0.0, -0.125 * (-1.0 + xi), 0.0);
          case 8:
            return sign * RealGradient(-0.125 * (1.0 - eta), 0.0, 0.0);
          case 9:
            return sign * RealGradient(0.0, -0.125 * (1.0 + xi), 0.0);
          case 10:
            return sign * RealGradient(0.125 * (1.0 + eta), 0.0, 0.0);
          case 11:
            return sign * RealGradient(0.0, -0.125 * (1.0 - xi), 0.0);

          default:
            libmesh_error_msg("Invalid i = " << i);
          } // switch(i)

        } // j = 2

        default:
          libmesh_error_msg("Invalid j = " << j);
        }
      }

      case TET10:
      case TET14:
      {
        switch (j)
        {
          // d()/dxi
        case 0:
        {
          switch (i)
          {
          case 0:
            return sign * RealGradient(0.0, -1.0, -1.0);
          case 1:
            return sign * RealGradient(0.0, -1.0, 0.0);
          case 2:
            return sign * RealGradient(0.0, 1.0, 0.0);
          case 3:
            return sign * RealGradient(0.0, 0.0, 1.0);
          case 4:
            return sign * RealGradient(0.0, 0.0, -1.0);
          case 5:
            return RealGradient();

          default:
            libmesh_error_msg("Invalid i = " << i);
          } // switch(i)

        } // j = 0

          // d()/deta
        case 1:
        {
          switch (i)
          {
          case 0:
          case 1:
            return sign * RealGradient(1.0, 0.0, 0.0);
          case 2:
            return sign * RealGradient(-1.0, 0.0, -1.0);
          case 3:
            return sign * RealGradient(0.0, 0.0, 1.0);
          case 4:
            return RealGradient();
          case 5:
            return sign * RealGradient(0.0, 0.0, -1.0);

          default:
            libmesh_error_msg("Invalid i = " << i);
          } // switch(i)

        } // j = 1

          // d()/dzeta
        case 2:
        {
          switch (i)
          {
          case 0:
            return sign * RealGradient(1.0, 0.0, 0.0);
          case 1:
            return RealGradient();
          case 2:
          case 5:
            return sign * RealGradient(0.0, 1.0, 0.0);
          case 3:
            return sign * RealGradient(-1.0, -1.0, 0.0);
          case 4:
            return sign * RealGradient(1.0, 0.0, 0.0);

          default:
            libmesh_error_msg("Invalid i = " << i);
          } // switch(i)

        } // j = 2

        default:
          libmesh_error_msg("Invalid j = " << j);
        }
      }
      case PRISM15:
      case PRISM18:
      {
        libmesh_assert_less(i, 9);

        const Real xi = p(0);
        const Real eta = p(1);
        const Real zeta = p(2);

        libmesh_assert_greater_equal(xi, -10. * TOLERANCE);
        libmesh_assert_greater_equal(eta, -10. * TOLERANCE);
        libmesh_assert_less_equal(xi + eta, 1.0 + 10. * TOLERANCE);
        libmesh_assert_less_equal(std::fabs(zeta), 1.0 + 10. * TOLERANCE);

        const Real f0 = 0.5 * (1.0 - zeta);
        const Real f1 = 0.5 * (1.0 + zeta);
        const Real df0 = -0.5;
        const Real df1 = 0.5;

        RealGradient dv;
        switch (j)
        {
          // d()/dxi
        case 0:
        {
          switch (i)
          {
          case 0:
            return sign * RealGradient(0.0, f0, 0.0);
          case 1:
            return sign * RealGradient(0.0, -1.0, 0.0);
          case 2:
            return sign * RealGradient(0.0, 1.0, 0.0);
          case 3:
            return sign * RealGradient(0.0, 0.0, 1.0);
          case 4:
            return sign * RealGradient(0.0, 0.0, -1.0);
          case 5:
            return RealGradient();

          default:
            libmesh_error_msg("Invalid i = " << i);
          } // switch(i)

        } // j = 0

          // d()/deta
        case 1:
        {
          switch (i)
          {
          case 0:
            return sign * RealGradient(-f0, 0.0, 0.0);
          case 1:
            return sign * RealGradient(1.0, 0.0, 0.0);
          case 2:
            return sign * RealGradient(-1.0, 0.0, -1.0);
          case 3:
            return sign * RealGradient(0.0, 0.0, 1.0);
          case 4:
            return RealGradient();
          case 5:
            return sign * RealGradient(0.0, 0.0, -1.0);

          default:
            libmesh_error_msg("Invalid i = " << i);
          } // switch(i)

        } // j = 1

          // d()/dzeta
        case 2:
        {
          switch (i)
          {
          case 0:
            return sign * RealGradient(df0 * (1.0 - eta), df0 * xi, 0.0);
          case 1:
            return RealGradient();
          case 2:
          case 5:
            return sign * RealGradient(0.0, 1.0, 0.0);
          case 3:
            return sign * RealGradient(-1.0, -1.0, 0.0);
          case 4:
            return sign * RealGradient(1.0, 0.0, 0.0);

          default:
            libmesh_error_msg("Invalid i = " << i);
          } // switch(i)

        } // j = 2

        default:
          libmesh_error_msg("Invalid j = " << j);
        }
        switch (i)
        {
        case 0:
          switch (j)
          {
          case 1:
            dv =
                break;
          default:
            dv =
                break;
          }
          break;

        case 1:
          switch (j)
          {
          case 0:
            dv = RealGradient(0.0, f0, 0.0);
            break;
          case 1:
            dv = RealGradient(-f0, 0.0, 0.0);
            break;
          default:
            dv = RealGradient(df0 * (-eta), df0 * xi, 0.0);
            break;
          }
          break;

        case 2:
          switch (j)
          {
          case 0:
            dv = RealGradient(0.0, -f0, 0.0);
            break;
          case 1:
            dv = RealGradient(f0, 0.0, 0.0);
            break;
          default:
            dv = RealGradient(df0 * (eta), df0 * (1.0 - xi), 0.0);
            break;
          }
          break;

        case 3:
          switch (j)
          {
          case 0:
            dv = RealGradient(0.0, 0.0, -0.5);
            break;
          case 1:
            dv = RealGradient(0.0, 0.0, -0.5);
            break;
          default:
            dv = RealGradient();
            break;
          }
          break;

        case 4:
          switch (j)
          {
          case 0:
            dv = RealGradient(0.0, 0.0, 0.5);
            break;
          default:
            dv = RealGradient();
            break;
          }
          break;

        case 5:
          switch (j)
          {
          case 1:
            dv = RealGradient(0.0, 0.0, 0.5);
            break;
          default:
            dv = RealGradient();
            break;
          }
          break;

        case 6:
          switch (j)
          {
          case 0:
            dv = RealGradient(0.0, f1, 0.0);
            break;
          case 1:
            dv = RealGradient(-f1, 0.0, 0.0);
            break;
          default:
            dv = RealGradient(df1 * (1.0 - eta), df1 * xi, 0.0);
            break;
          }
          break;

        case 7:
          switch (j)
          {
          case 0:
            dv = RealGradient(0.0, f1, 0.0);
            break;
          case 1:
            dv = RealGradient(-f1, 0.0, 0.0);
            break;
          default:
            dv = RealGradient(df1 * (-eta), df1 * xi, 0.0);
            break;
          }
          break;

        case 8:
          switch (j)
          {
          case 0:
            dv = RealGradient(0.0, -f1, 0.0);
            break;
          case 1:
            dv = RealGradient(f1, 0.0, 0.0);
            break;
          default:
            dv = RealGradient(df1 * (eta), df1 * (1.0 - xi), 0.0);
            break;
          }
          break;

        default:
          libmesh_error_msg("Invalid i = " << i);
        }

        if (!elem->positive_edge_orientation(i))
          dv = -dv;
        return dv;
      }

      case PYRAMID13:
      case PYRAMID14:
      {
        libmesh_assert_less(i, 8);

        const Real xi = p(0);
        const Real eta = p(1);
        const Real zeta = p(2);

        libmesh_assert_less_equal(std::fabs(xi), 1.0 + 10. * TOLERANCE);
        libmesh_assert_less_equal(std::fabs(eta), 1.0 + 10. * TOLERANCE);
        libmesh_assert_greater_equal(zeta, -10. * TOLERANCE);
        libmesh_assert_less_equal(zeta, 1.0 + 10. * TOLERANCE);

        Real a = 1.0 - zeta;
        if (std::fabs(a) < 1e-14)
          a = (a >= 0.0) ? 1e-14 : -1e-14;

        const Real inva = 1.0 / a;
        const Real inva2 = inva * inva;
        const Real inva3 = inva2 * inva;

        const Real l0 = 0.25 * (a - xi - eta + xi * eta * inva);
        const Real l1 = 0.25 * (a + xi - eta - xi * eta * inva);
        const Real l2 = 0.25 * (a + xi + eta + xi * eta * inva);
        const Real l3 = 0.25 * (a - xi + eta - xi * eta * inva);
        const Real l4 = zeta;

        const RealGradient g0(0.25 * (-1.0 + eta * inva), 0.25 * (-1.0 + xi * inva), 0.25 * (-1.0 + xi * eta * inva2));
        const RealGradient g1(0.25 * (1.0 - eta * inva), 0.25 * (-1.0 - xi * inva), 0.25 * (-1.0 - xi * eta * inva2));
        const RealGradient g2(0.25 * (1.0 + eta * inva), 0.25 * (1.0 + xi * inva), 0.25 * (-1.0 + xi * eta * inva2));
        const RealGradient g3(0.25 * (-1.0 - eta * inva), 0.25 * (1.0 - xi * inva), 0.25 * (-1.0 - xi * eta * inva2));
        const RealGradient g4(0.0, 0.0, 1.0);

        RealGradient dg0, dg1, dg2, dg3, dg4;
        dg4 = RealGradient();
        if (j == 0)
        {
          dg0 = RealGradient(0.0, 0.25 * inva, 0.25 * eta * inva2);
          dg1 = RealGradient(0.0, -0.25 * inva, -0.25 * eta * inva2);
          dg2 = RealGradient(0.0, 0.25 * inva, 0.25 * eta * inva2);
          dg3 = RealGradient(0.0, -0.25 * inva, -0.25 * eta * inva2);
        }
        else if (j == 1)
        {
          dg0 = RealGradient(0.25 * inva, 0.0, 0.25 * xi * inva2);
          dg1 = RealGradient(-0.25 * inva, 0.0, -0.25 * xi * inva2);
          dg2 = RealGradient(0.25 * inva, 0.0, 0.25 * xi * inva2);
          dg3 = RealGradient(-0.25 * inva, 0.0, -0.25 * xi * inva2);
        }
        else
        {
          dg0 = RealGradient(0.25 * eta * inva2, 0.25 * xi * inva2, 0.5 * xi * eta * inva3);
          dg1 = RealGradient(-0.25 * eta * inva2, -0.25 * xi * inva2, -0.5 * xi * eta * inva3);
          dg2 = RealGradient(0.25 * eta * inva2, 0.25 * xi * inva2, 0.5 * xi * eta * inva3);
          dg3 = RealGradient(-0.25 * eta * inva2, -0.25 * xi * inva2, -0.5 * xi * eta * inva3);
        }

        auto deriv_whitney = [&](const Real la, const Real lb,
                                 const RealGradient &ga, const RealGradient &gb,
                                 const RealGradient &dga, const RealGradient &dgb) -> RealGradient
        {
          const Real dla = ga(j);
          const Real dlb = gb(j);
          return dla * gb + la * dgb - dlb * ga - lb * dga;
        };

        RealGradient dv;
        switch (i)
        {
        case 0:
          dv = deriv_whitney(l0, l1, g0, g1, dg0, dg1);
          break; // (0,1)
        case 1:
          dv = deriv_whitney(l1, l2, g1, g2, dg1, dg2);
          break; // (1,2)
        case 2:
          dv = deriv_whitney(l2, l3, g2, g3, dg2, dg3);
          break; // (2,3)
        case 3:
          dv = -(deriv_whitney(l0, l3, g0, g3, dg0, dg3));
          break; // N = l3*g0 - l0*g3
        case 4:
          dv = deriv_whitney(l0, l4, g0, g4, dg0, dg4);
          break; // (0,4)
        case 5:
          dv = deriv_whitney(l1, l4, g1, g4, dg1, dg4);
          break; // (1,4)
        case 6:
          dv = deriv_whitney(l2, l4, g2, g4, dg2, dg4);
          break; // (2,4)
        case 7:
          dv = deriv_whitney(l3, l4, g3, g4, dg3, dg4);
          break; // (3,4)
        default:
          libmesh_error_msg("Invalid i = " << i);
        }

        libmesh_ignore(l4);

        if (!elem->positive_edge_orientation(i))
          dv = -dv;
        return dv;
      }
      default:
        libmesh_error_msg("ERROR: Unsupported 3D element type!: " << Utility::enum_to_string(elem->type()));
      }
    }
      // unsupported order
    default:
      libmesh_error_msg("ERROR: Unsupported 3D FE order!: " << totalorder);
    }

#else // LIBMESH_DIM != 3
    libmesh_ignore(elem, order, i, j, p, add_p_level);
    libmesh_not_implemented();
#endif
  }

  template <>
  RealGradient FE<3, NEDELEC_ONE>::shape_deriv(const ElemType,
                                               const Order,
                                               const unsigned int,
                                               const unsigned int,
                                               const Point &)
  {
    libmesh_error_msg("Nedelec elements require the element type \nbecause edge orientation is needed.");
    return RealGradient();
  }

  template <>
  RealGradient FE<3, NEDELEC_ONE>::shape_deriv(const FEType fet,
                                               const Elem *elem,
                                               const unsigned int i,
                                               const unsigned int j,
                                               const Point &p,
                                               const bool add_p_level)
  {
    return FE<3, NEDELEC_ONE>::shape_deriv(elem, fet.order, i, j, p, add_p_level);
  }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  template <>
  RealGradient FE<3, NEDELEC_ONE>::shape_second_deriv(const Elem *elem,
                                                      const Order order,
                                                      const unsigned int i,
                                                      const unsigned int j,
                                                      const Point &,
                                                      const bool add_p_level)
  {
#if LIBMESH_DIM == 3

    libmesh_assert(elem);

    // j = 0 ==> d^2 phi / dxi^2
    // j = 1 ==> d^2 phi / dxi deta
    // j = 2 ==> d^2 phi / deta^2
    // j = 3 ==> d^2 phi / dxi dzeta
    // j = 4 ==> d^2 phi / deta dzeta
    // j = 5 ==> d^2 phi / dzeta^2
    libmesh_assert_less(j, 6);

    const Order totalorder = order + add_p_level * elem->p_level();
    libmesh_assert_less(i, n_dofs(elem->type(), totalorder));

    const char sign = elem->positive_edge_orientation(i) ? 1 : -1;

    switch (totalorder)
    {
      // linear Nedelec (first kind) shape function second derivatives
    case FIRST:
    {
      switch (elem->type())
      {
      case HEX20:
      case HEX27:
      {
        switch (j)
        {
          // d^2()/dxi^2, d^2()/deta^2, d^2()/dzeta^2
        case 0:
        case 2:
        case 5:
          return RealGradient();

          // d^2()/dxideta
        case 1:
        {
          switch (i)
          {
          case 0:
          case 1:
          case 2:
          case 3:
          case 8:
          case 9:
          case 10:
          case 11:
            return RealGradient();
          case 4:
          case 6:
            return sign * RealGradient(0.0, 0.0, -0.125);
          case 5:
          case 7:
            return sign * RealGradient(0.0, 0.0, 0.125);

          default:
            libmesh_error_msg("Invalid i = " << i);
          } // switch(i)

        } // j = 1

          // d^2()/dxidzeta
        case 3:
        {
          switch (i)
          {
          case 0:
          case 2:
          case 4:
          case 5:
          case 6:
          case 7:
          case 8:
          case 10:
            return RealGradient();
          case 1:
          case 3:
          case 11:
            return sign * RealGradient(0.0, 0.125);
          case 9:
            return sign * RealGradient(0.0, -0.125, 0.0);

          default:
            libmesh_error_msg("Invalid i = " << i);
          } // switch(i)

        } // j = 3

          // d^2()/detadzeta
        case 4:
        {
          switch (i)
          {
          case 0:
            return sign * RealGradient(-0.125, 0.0, 0.0);
          case 1:
          case 3:
          case 4:
          case 5:
          case 6:
          case 7:
          case 9:
          case 11:
            return RealGradient();
          case 2:
          case 8:
          case 10:
            return sign * RealGradient(0.125, 0.0, 0.0);

          default:
            libmesh_error_msg("Invalid i = " << i);
          } // switch(i)

        } // j = 4

        default:
          libmesh_error_msg("Invalid j = " << j);
        }
      }

        // All second derivatives for linear tets are zero.
      case TET10:
      case TET14:
        return RealGradient();

      default:
        libmesh_error_msg("ERROR: Unsupported 3D element type!: " << Utility::enum_to_string(elem->type()));

      } // switch(type)
    }

      // unsupported order
    default:
      libmesh_error_msg("ERROR: Unsupported 3D FE order!: " << totalorder);
    }

#else // LIBMESH_DIM != 3
    libmesh_assert(true || p(0));
    libmesh_ignore(elem, order, i, j, add_p_level);
    libmesh_not_implemented();
#endif
  }

  template <>
  RealGradient FE<3, NEDELEC_ONE>::shape_second_deriv(const ElemType,
                                                      const Order,
                                                      const unsigned int,
                                                      const unsigned int,
                                                      const Point &)
  {
    libmesh_error_msg("Nedelec elements require the element type \nbecause edge orientation is needed.");
    return RealGradient();
  }

  template <>
  RealGradient FE<3, NEDELEC_ONE>::shape_second_deriv(const FEType fet,
                                                      const Elem *elem,
                                                      const unsigned int i,
                                                      const unsigned int j,
                                                      const Point &p,
                                                      const bool add_p_level)
  {
    return FE<3, NEDELEC_ONE>::shape_second_deriv(elem, fet.order, i, j, p, add_p_level);
  }

#endif

} // namespace libMesh
