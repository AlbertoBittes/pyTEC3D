// shape_second_deriv() of NEDELEC_ONE on 3D cells against central differences of shape_deriv().
//
// libMesh is normally built without second derivatives (--disable-second), so this test compiles
// the patched fe_nedelec_one_shape_3D.C itself, with LIBMESH_ENABLE_SECOND_DERIVATIVES defined:
//
//   mpicxx -DLIBMESH_ENABLE_SECOND_DERIVATIVES -DNEDELEC_ONE_SHAPE_3D=\"<libMesh src>/src/fe/fe_nedelec_one_shape_3D.C\" ...
//
// Its shape functions take precedence over those of the shared libMesh library. Only the static
// shape functions are used here, so the other effects of the define on libMesh classes are moot.

#ifndef LIBMESH_ENABLE_SECOND_DERIVATIVES
#error "compile with -DLIBMESH_ENABLE_SECOND_DERIVATIVES"
#endif

#include NEDELEC_ONE_SHAPE_3D

#include "libmesh/libmesh.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/fe_interface.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/quadrature_gauss.h"

#include <iomanip>
#include <iostream>

using namespace libMesh;

int main(int argc, char ** argv)
{
  LibMeshInit init(argc, argv);

  // j of shape_second_deriv() -> the two derivatives
  const unsigned int ab[6][2] = {{0, 0}, {0, 1}, {1, 1}, {0, 2}, {1, 2}, {2, 2}};
  const Real h = 1.e-6;
  unsigned int n_failures = 0, n_checks = 0;

  const std::vector<std::pair<Order, std::vector<ElemType>>> cases =
    {{FIRST, {TET10, TET14, HEX20, HEX27, PRISM15, PRISM18, PRISM20, PRISM21,
              PYRAMID13, PYRAMID14, PYRAMID18}},
     {SECOND, {TET14, HEX27, PRISM20, PRISM21, PYRAMID18}}};

  std::cout << "second derivatives: shape_second_deriv against central differences of shape_deriv"
            << std::endl;
  for (const auto & [order, types] : cases)
    for (const ElemType t : types)
      {
        Mesh mesh(init.comm());
        MeshTools::Generation::build_cube(mesh, 1, 1, 1, 0., 1., 0., 1., 0., 1., t);
        Elem * elem = *mesh.elements_begin();
        QGauss qrule(3, FOURTH);
        qrule.init(*elem);
        Real worst = 0;

        for (const auto perm : make_range(elem->n_permutations()))
          {
            elem->permute(perm);
            const unsigned int n_dofs = FEInterface::n_dofs(FEType(order, NEDELEC_ONE), elem);
            for (const Point & p : qrule.get_points())
              for (const auto i : make_range(n_dofs))
                for (const auto j : make_range(6u))
                  {
                    Point dp;
                    dp(ab[j][1]) = h;
                    const RealGradient fd =
                      (FE<3,NEDELEC_ONE>::shape_deriv(elem, order, i, ab[j][0], p + dp) -
                       FE<3,NEDELEC_ONE>::shape_deriv(elem, order, i, ab[j][0], p - dp)) / (2*h);
                    const RealGradient d2 = FE<3,NEDELEC_ONE>::shape_second_deriv(elem, order, i, j, p);
                    const Real err = (fd - d2).norm() / std::max(1., d2.norm());
                    worst = std::max(worst, err);
                    ++n_checks;
                    if (err > 1.e-6)
                      {
                        if (++n_failures < 20)
                          std::cout << "  FAIL: " << Utility::enum_to_string(t) << " "
                                    << Utility::enum_to_string(order) << " i=" << i << " j=" << j
                                    << " error " << err << std::endl;
                      }
                  }
          }
        std::cout << "  " << std::setw(18) << std::left
                  << (Utility::enum_to_string(t) + " " + Utility::enum_to_string(order))
                  << " worst relative error " << worst << std::endl;
      }

  std::cout << (n_failures ? "FAILED: " : "PASSED: ") << n_failures << " failures in " << n_checks
            << " checks" << std::endl;
  return n_failures ? 1 : 0;
}
