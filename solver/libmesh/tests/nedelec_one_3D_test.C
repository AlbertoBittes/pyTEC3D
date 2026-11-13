// Tests of the NEDELEC_ONE shape functions on 3D cells (libMesh patched by
// solver/libmesh/patches/patch_nedelec_one_3D.py):
//
//   derivatives   shape_deriv() against central differences of shape(), for every permutation
//                 of the element vertices (every orientation of its edges and faces)
//   conformity    continuity of the tangential components of the physical shape functions across
//                 every interior face, on meshes of one element type and on a hybrid mesh, with
//                 the elements permuted; a shape function of a degree of freedom the neighbor does
//                 not share must have a vanishing tangential trace on the face
//   convergence   curl curl E + E = f on a skewed, permuted unit cube, with E a smooth field and
//                 natural boundary conditions: the L2 errors of E and curl E must decrease as h^p
//
// Usage: nedelec_one_3D_test [derivatives] [conformity] [convergence] [-ksp_type preonly -pc_type lu]

#include "libmesh/libmesh.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_map.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/tensor_value.h"

#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace libMesh;

namespace
{
const Parallel::Communicator * comm_world = nullptr;
unsigned int n_failures = 0;
unsigned int n_checks = 0;

void check(const bool ok, const std::string & what)
{
  ++n_checks;
  if (!ok)
    {
      ++n_failures;
      if (n_failures < 50)
        std::cout << "  FAIL: " << what << std::endl;
    }
}

std::string name(const ElemType t, const Order o)
{
  return Utility::enum_to_string(t) + " " + Utility::enum_to_string(o);
}

/// The element types with NEDELEC_ONE shape functions of order o
std::vector<ElemType> types(const Order o)
{
  if (o == FIRST)
    return {TET10, TET14, HEX20, HEX27, PRISM15, PRISM18, PRISM20, PRISM21,
            PYRAMID13, PYRAMID14, PYRAMID18};
  return {TET14, HEX27, PRISM20, PRISM21, PYRAMID18};
}

/// An affine skew, so that no symmetry of the mesh masks a wrong orientation
void skew(MeshBase & mesh)
{
  for (auto & node : mesh.node_ptr_range())
    {
      const Point p = *node;
      (*node)(0) = 1.0*p(0) + 0.2*p(1) + 0.1*p(2);
      (*node)(1) = 0.1*p(0) + 1.1*p(1) - 0.15*p(2);
      (*node)(2) = 0.05*p(0) + 0.1*p(1) + 0.9*p(2);
    }
}

/// Renumbers the vertices of every element, by its id and a round, without moving them
void permute(MeshBase & mesh, const unsigned int round)
{
  for (auto * elem : mesh.element_ptr_range())
    elem->permute((round + 7*elem->id()) % elem->n_permutations());
  mesh.find_neighbors();
}

// ------------------------------------------------------------------------------------------------
void test_derivatives()
{
  std::cout << "derivatives: shape_deriv against central differences of shape" << std::endl;
  const Real h = 1.e-6;

  for (const Order order : {FIRST, SECOND})
    for (const ElemType t : types(order))
      {
        Mesh mesh(*comm_world);
        MeshTools::Generation::build_cube(mesh, 1, 1, 1, 0., 1., 0., 1., 0., 1., t);
        skew(mesh);
        Real worst = 0;
        unsigned int n = 0;

        for (auto * elem : mesh.element_ptr_range())
          {
            QGauss qrule(3, FOURTH);
            qrule.init(*elem);
            for (const auto perm : make_range(elem->n_permutations()))
              {
                elem->permute(perm);
                const unsigned int n_dofs = FEInterface::n_dofs(FEType(order, NEDELEC_ONE), elem);
                for (const Point & p : qrule.get_points())
                  for (const auto i : make_range(n_dofs))
                    for (const auto j : make_range(3u))
                      {
                        Point dp;
                        dp(j) = h;
                        const RealGradient fd =
                          (FE<3,NEDELEC_ONE>::shape(elem, order, i, p + dp) -
                           FE<3,NEDELEC_ONE>::shape(elem, order, i, p - dp)) / (2*h);
                        const RealGradient d = FE<3,NEDELEC_ONE>::shape_deriv(elem, order, i, j, p);
                        const Real err = (fd - d).norm() / std::max(1., d.norm());
                        worst = std::max(worst, err);
                        ++n;
                        check(err < 1.e-6, name(t, order) + " i=" + std::to_string(i) +
                              " j=" + std::to_string(j));
                      }
              }
            // one element per type is enough for the reference shapes, but cover all its
            // orientations
            break;
          }
        std::cout << "  " << std::setw(18) << std::left << name(t, order)
                  << " " << n << " values, worst relative error " << worst << std::endl;
      }
}

/// The inverse transpose of the (constant) Jacobian of an affine element, from its vertices
RealTensor affine_inverse_transpose_jacobian(const Elem & elem)
{
  for (const auto a : make_range(1u, elem.n_vertices()))
    for (const auto b : make_range(a + 1, elem.n_vertices()))
      for (const auto c : make_range(b + 1, elem.n_vertices()))
        {
          const Point dxi[3] = {elem.master_point(a) - elem.master_point(0),
                                elem.master_point(b) - elem.master_point(0),
                                elem.master_point(c) - elem.master_point(0)};
          const Point dx[3] = {elem.point(a) - elem.point(0), elem.point(b) - elem.point(0),
                               elem.point(c) - elem.point(0)};
          RealTensor Xi, Xp;
          for (const auto r : make_range(3u))
            for (const auto k : make_range(3u))
              {
                Xi(r, k) = dxi[k](r);
                Xp(r, k) = dx[k](r);
              }
          if (std::abs(Xi.det()) > 1.e-8)
            return (Xp * Xi.inverse()).inverse().transpose();
        }
  libmesh_error_msg("degenerate element");
}

// ------------------------------------------------------------------------------------------------
/// Tangential conformity across the interior faces of mesh
void check_conformity(MeshBase & mesh, const Order order, const std::string & label)
{
  const FEType fe_type(order, NEDELEC_ONE);
  Real worst_shared = 0, worst_unshared = 0;
  unsigned int n_faces = 0;

  for (const auto round : make_range(3u))
    {
      permute(mesh, round);

      EquationSystems es(mesh);
      System & sys = es.add_system<System>("conformity");
      sys.add_variable("E", fe_type);
      es.init();
      const DofMap & dof_map = sys.get_dof_map();


      for (const auto * elem : mesh.active_element_ptr_range())
        {
          std::vector<dof_id_type> dofs;
          dof_map.dof_indices(elem, dofs);

          for (const auto s : elem->side_index_range())
            {
              const Elem * neighbor = elem->neighbor_ptr(s);
              if (!neighbor)
                continue;
              ++n_faces;

              std::vector<dof_id_type> neighbor_dofs;
              dof_map.dof_indices(neighbor, neighbor_dofs);

              std::unique_ptr<const Elem> side = elem->build_side_ptr(s);
              QGauss side_rule(2, FIFTH);
              side_rule.init(*side);

              // The reference points of the side quadrature points in both elements, from the
              // weights of the shared vertices: exact, since the faces are affine in reference
              // coordinates (inverse_map() would add the error of the geometry map)
              const ElemType vertex_type = side->n_vertices() == 3 ? TRI3 : QUAD4;
              std::vector<Point> points, neighbor_points;
              for (const auto & q : side_rule.get_points())
                {
                  Point a, b;
                  for (const auto v : make_range(side->n_vertices()))
                    {
                      const Real w = FE<2,LAGRANGE>::shape(vertex_type, FIRST, v, q);
                      const unsigned int n = elem->local_side_node(s, v);
                      a += w * elem->master_point(n);
                      b += w * neighbor->master_point(neighbor->local_node(elem->node_id(n)));
                    }
                  points.push_back(a);
                  neighbor_points.push_back(b);
                  check((FEMap::map(3, elem, a) - FEMap::map(3, neighbor, b)).norm() < 1.e-12,
                        label + ": side point construction");
                }
              // The physical shape functions, covariant: J^-T phi, with the exact Jacobian of the
              // affine elements (that of the geometry map of PYRAMID18 is not exact)
              const RealTensor JiT = affine_inverse_transpose_jacobian(*elem),
                               JiT_neighbor = affine_inverse_transpose_jacobian(*neighbor);
              std::vector<std::vector<RealGradient>> phi(dofs.size()), phi_neighbor(neighbor_dofs.size());
              for (const auto i : index_range(dofs))
                for (const auto & a : points)
                  phi[i].push_back(JiT * FE<3,NEDELEC_ONE>::shape(elem, order, i, a));
              for (const auto j : index_range(neighbor_dofs))
                for (const auto & b : neighbor_points)
                  phi_neighbor[j].push_back(JiT_neighbor * FE<3,NEDELEC_ONE>::shape(neighbor, order, j, b));

              // The unit normal of the (planar, after the affine skew) side
              const Point n0 = side->point(0), n1 = side->point(1), n2 = side->point(2);
              Point normal = (n1 - n0).cross(n2 - n0);
              normal /= normal.norm();

              for (const auto i : index_range(dofs))
                {
                  const auto it = std::find(neighbor_dofs.begin(), neighbor_dofs.end(), dofs[i]);
                  for (const auto q : index_range(points))
                    {
                      const RealGradient t = normal.cross(phi[i][q]);
                      if (it == neighbor_dofs.end())
                        {
                          worst_unshared = std::max(worst_unshared, t.norm());
                          check(t.norm() < 1.e-11, label + ": unshared dof with tangential trace");
                        }
                      else
                        {
                          const auto j = std::distance(neighbor_dofs.begin(), it);
                          const RealGradient tn = normal.cross(phi_neighbor[j][q]);
                          const Real err = (t - tn).norm() / std::max(1., t.norm());
                          worst_shared = std::max(worst_shared, err);
                          check(err < 1.e-11, label + ": tangential jump");
                        }
                    }
                }
            }
        }
    }
  std::cout << "  " << std::setw(30) << std::left << label << " " << n_faces
            << " face visits, worst tangential jump " << worst_shared
            << ", worst unshared trace " << worst_unshared << std::endl;
}

/// A hexahedron with a pyramid on its top and bottom faces, tetrahedra on the pyramid's triangles,
/// a prism on its side and a tetrahedron on the prism's top triangle: every pair of cell types
/// sharing a face, with second order (complete) elements
void build_hybrid(MeshBase & mesh)
{
  const std::vector<Point> points =
    {{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1},   // 0-7 hex
     {0.5,0.45,1.7}, {0.45,0.55,-0.75},                                         // 8, 9 apexes
     {0.5,-0.6,1.4}, {1.6,0.5,1.45}, {0.5,1.6,1.35}, {-0.6,0.5,1.4},            // 10-13 tets on top
     {1.8,0.4,0}, {1.8,0.4,1},                                                // 14, 15 prism
     {1.45,0.3,1.7}};                                                           // 16 tet on prism
  for (const auto i : index_range(points))
    mesh.add_point(points[i], i);

  const std::vector<std::pair<ElemType, std::vector<dof_id_type>>> elems =
    {{HEX8,     {0, 1, 2, 3, 4, 5, 6, 7}},
     {PYRAMID5, {4, 5, 6, 7, 8}},
     {PYRAMID5, {0, 3, 2, 1, 9}},
     {TET4,     {4, 5, 8, 10}},
     {TET4,     {5, 6, 8, 11}},
     {TET4,     {6, 7, 8, 12}},
     {TET4,     {7, 4, 8, 13}},
     {PRISM6,   {1, 2, 14, 5, 6, 15}},
     {TET4,     {5, 6, 15, 16}}};
  for (const auto & [type, nodes] : elems)
    {
      Elem * elem = mesh.add_elem(Elem::build(type));
      for (const auto n : index_range(nodes))
        elem->set_node(n, mesh.node_ptr(nodes[n]));
      // a positive volume
      if (type == TET4 || type == PRISM6)
        {
          const Point a = elem->point(1) - elem->point(0), b = elem->point(2) - elem->point(0),
                      c = elem->point(type == TET4 ? 3 : 3) - elem->point(0);
          if (a.cross(b) * c < 0)
            {
              if (type == TET4)
                {
                  Node * tmp = elem->node_ptr(1);
                  elem->set_node(1, elem->node_ptr(2));
                  elem->set_node(2, tmp);
                }
              else
                for (const unsigned int k : {0u, 3u})
                  {
                    Node * tmp = elem->node_ptr(k + 1);
                    elem->set_node(k + 1, elem->node_ptr(k + 2));
                    elem->set_node(k + 2, tmp);
                  }
            }
        }
    }
  mesh.prepare_for_use();
  mesh.all_complete_order();
}

void test_conformity()
{
  std::cout << "conformity: tangential continuity across interior faces" << std::endl;
  for (const Order order : {FIRST, SECOND})
    {
      for (const ElemType t : types(order))
        {
          Mesh mesh(*comm_world);
          MeshTools::Generation::build_cube(mesh, 2, 2, 2, 0., 1., 0., 1., 0., 1., t);
          skew(mesh);
          mesh.allow_renumbering(false);
          check_conformity(mesh, order, name(t, order));
        }
      Mesh mesh(*comm_world);
      build_hybrid(mesh);
      skew(mesh);
      mesh.allow_renumbering(false);
      check_conformity(mesh, order, std::string("hybrid ") + Utility::enum_to_string(order));
    }
}

// ------------------------------------------------------------------------------------------------
// The exact field: a solenoidal part with curl curl E_s = 2 pi^2 E_s and a gradient
const Real pi = libMesh::pi;
const Point g_dir(0.5, 0.3, -0.2);

RealGradient exact(const Point & x)
{
  const Real sx = std::sin(pi*x(0)), sy = std::sin(pi*x(1)), sz = std::sin(pi*x(2));
  const Real g = std::exp(g_dir * x);
  return RealGradient(sy*sz, sz*sx, sx*sy) + g*g_dir;
}

RealGradient exact_curl(const Point & x)
{
  const Real sx = std::sin(pi*x(0)), sy = std::sin(pi*x(1)), sz = std::sin(pi*x(2));
  const Real cx = std::cos(pi*x(0)), cy = std::cos(pi*x(1)), cz = std::cos(pi*x(2));
  return pi*RealGradient(sx*(cy - cz), sy*(cz - cx), sz*(cx - cy));
}

RealGradient source(const Point & x)
{
  const Real sx = std::sin(pi*x(0)), sy = std::sin(pi*x(1)), sz = std::sin(pi*x(2));
  const Real g = std::exp(g_dir * x);
  return (2*pi*pi + 1)*RealGradient(sy*sz, sz*sx, sx*sy) + g*g_dir;
}

void assemble_curl_curl(EquationSystems & es, const std::string &)
{
  const MeshBase & mesh = es.get_mesh();
  LinearImplicitSystem & sys = es.get_system<LinearImplicitSystem>("curl_curl");
  const DofMap & dof_map = sys.get_dof_map();
  const FEType fe_type = dof_map.variable_type(0);

  std::unique_ptr<FEVectorBase> fe = FEVectorBase::build(3, fe_type);
  QGauss qrule(3, static_cast<Order>(2*fe_type.order + 3));
  fe->attach_quadrature_rule(&qrule);
  const auto & phi = fe->get_phi();
  const auto & curl_phi = fe->get_curl_phi();
  const auto & JxW = fe->get_JxW();
  const auto & xyz = fe->get_xyz();

  std::unique_ptr<FEVectorBase> fe_face = FEVectorBase::build(3, fe_type);
  QGauss qface(2, static_cast<Order>(2*fe_type.order + 3));
  fe_face->attach_quadrature_rule(&qface);
  const auto & phi_face = fe_face->get_phi();
  const auto & JxW_face = fe_face->get_JxW();
  const auto & xyz_face = fe_face->get_xyz();
  const auto & normals = fe_face->get_normals();

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;
  std::vector<dof_id_type> dofs;

  for (const auto * elem : mesh.active_local_element_ptr_range())
    {
      dof_map.dof_indices(elem, dofs);
      const unsigned int n = dofs.size();
      Ke.resize(n, n);
      Fe.resize(n);
      fe->reinit(elem);
      for (const auto qp : index_range(JxW))
        {
          const RealGradient f = source(xyz[qp]);
          for (const auto i : make_range(n))
            {
              Fe(i) += JxW[qp] * (f * phi[i][qp]);
              for (const auto j : make_range(n))
                Ke(i, j) += JxW[qp] * (curl_phi[i][qp] * curl_phi[j][qp] + phi[i][qp] * phi[j][qp]);
            }
        }
      // natural boundary condition: - int (n x curl E) . v
      for (const auto s : elem->side_index_range())
        if (!elem->neighbor_ptr(s))
          {
            fe_face->reinit(elem, s);
            for (const auto qp : index_range(JxW_face))
              {
                const RealGradient g = normals[qp].cross(exact_curl(xyz_face[qp]));
                for (const auto i : make_range(n))
                  Fe(i) -= JxW_face[qp] * (g * phi_face[i][qp]);
              }
          }
      sys.matrix->add_matrix(Ke, dofs);
      sys.rhs->add_vector(Fe, dofs);
    }
}

/// L2 errors of E and curl E
std::pair<Real, Real> errors(EquationSystems & es)
{
  const MeshBase & mesh = es.get_mesh();
  LinearImplicitSystem & sys = es.get_system<LinearImplicitSystem>("curl_curl");
  const DofMap & dof_map = sys.get_dof_map();
  const FEType fe_type = dof_map.variable_type(0);

  std::unique_ptr<FEVectorBase> fe = FEVectorBase::build(3, fe_type);
  QGauss qrule(3, static_cast<Order>(2*fe_type.order + 5));
  fe->attach_quadrature_rule(&qrule);
  const auto & phi = fe->get_phi();
  const auto & curl_phi = fe->get_curl_phi();
  const auto & JxW = fe->get_JxW();
  const auto & xyz = fe->get_xyz();

  Real e0 = 0, e1 = 0;
  std::vector<dof_id_type> dofs;
  for (const auto * elem : mesh.active_local_element_ptr_range())
    {
      dof_map.dof_indices(elem, dofs);
      fe->reinit(elem);
      for (const auto qp : index_range(JxW))
        {
          RealGradient u, cu;
          for (const auto i : index_range(dofs))
            {
              const Number c = sys.current_solution(dofs[i]);
              u += c * phi[i][qp];
              cu += c * curl_phi[i][qp];
            }
          e0 += JxW[qp] * (u - exact(xyz[qp])).norm_sq();
          e1 += JxW[qp] * (cu - exact_curl(xyz[qp])).norm_sq();
        }
    }
  mesh.comm().sum(e0);
  mesh.comm().sum(e1);
  return {std::sqrt(e0), std::sqrt(e1)};
}

void test_convergence()
{
  std::cout << "convergence: curl curl E + E = f, L2 errors of E and curl E (rate)" << std::endl;
  for (const Order order : {FIRST, SECOND})
    for (const ElemType t : types(order))
      {
        std::vector<Real> e0, e1;
        std::cout << "  " << std::setw(18) << std::left << name(t, order);
        for (const unsigned int n : {2u, 4u, 8u})
          {
            Mesh mesh(*comm_world);
            MeshTools::Generation::build_cube(mesh, n, n, n, 0., 1., 0., 1., 0., 1., t);
            skew(mesh);
            mesh.allow_renumbering(false);
            permute(mesh, n);

            EquationSystems es(mesh);
            LinearImplicitSystem & sys = es.add_system<LinearImplicitSystem>("curl_curl");
            sys.add_variable("E", order, NEDELEC_ONE);
            sys.attach_assemble_function(assemble_curl_curl);
            es.init();
            sys.solve();
            const auto [a, b] = errors(es);
            e0.push_back(a);
            e1.push_back(b);
            std::cout << std::scientific << std::setprecision(2) << "  n=" << n << " [" << sys.n_dofs()
                      << " dofs] " << a << " " << b;
            if (e0.size() > 1)
              std::cout << std::fixed << std::setprecision(2) << " ("
                        << std::log2(Real(e0[e0.size()-2]/a)) << ", " << std::log2(Real(e1[e1.size()-2]/b)) << ")";
          }
        std::cout << std::endl;
        const Real r0 = std::log2(e0[1]/e0[2]), r1 = std::log2(e1[1]/e1[2]);
        const Real expected = static_cast<int>(order);
        check(r0 > expected - 0.25 && r1 > expected - 0.25,
              name(t, order) + ": convergence rates " + std::to_string(r0) + ", " + std::to_string(r1));
      }
}
} // anonymous namespace

int main(int argc, char ** argv)
{
  LibMeshInit init(argc, argv);
  comm_world = &init.comm();

  std::vector<std::string> tests;
  for (int a = 1; a < argc; ++a)
    if (argv[a][0] != '-' && (a == 1 || argv[a-1][0] != '-'))
      tests.push_back(argv[a]);
  if (tests.empty())
    tests = {"derivatives", "conformity", "convergence"};

  for (const auto & t : tests)
    {
      if (t == "derivatives")
        test_derivatives();
      else if (t == "conformity")
        test_conformity();
      else if (t == "convergence")
        test_convergence();
      else
        libmesh_error_msg("Unknown test " << t);
    }

  std::cout << (n_failures ? "FAILED: " : "PASSED: ") << n_failures << " failures in "
            << n_checks << " checks" << std::endl;
  return n_failures ? 1 : 0;
}
