"""FIRST order NEDELEC_ONE shape functions of the prism and the pyramid (the pyramid: MFEM's
lowest order Fuentes basis, mapped to libMesh's pyramid, with libMesh's sign)."""
from nedelec_ref import *

h = R(1, 2)
prism1 = [vec(-h*(1-eta)*(1-zeta), -h*xi*(1-zeta), 0), vec(h*eta*(1-zeta), -h*xi*(1-zeta), 0),
          vec(-h*eta*(1-zeta), -h*(1-xi)*(1-zeta), 0), vec(0, 0, -h*(1-xi-eta)),
          vec(0, 0, -h*xi), vec(0, 0, -h*eta),
          vec(-h*(1-eta)*(1+zeta), -h*xi*(1+zeta), 0), vec(h*eta*(1+zeta), -h*xi*(1+zeta), 0),
          vec(-h*eta*(1+zeta), -h*(1-xi)*(1+zeta), 0)]

# pyramid: in the coordinates x, y, z of MFEM's pyramid, with libMesh's sign (-J^T phi_MFEM)
x = h*(1+xi-zeta); y = h*(1+eta-zeta); ox = h*(1-xi-zeta); oy = h*(1-eta-zeta)
ozi = 1/(1-zeta); q = 1-2*zeta
def pv(vx, vy, vz):
    return vec(-h*vx, -h*vy, -vz + h*(vx+vy))
pyr1 = [vec(-h*oy, 0, -x*oy*ozi + h*oy), vec(0, -h*x, -x*y*ozi + h*x),
        vec(h*y, 0, x*y*ozi - h*y), vec(0, -h*ox, -ox*y*ozi + h*ox),
        pv(oy*zeta*ozi, ox*zeta*ozi, 1-x-y+x*y*q*ozi**2),
        pv(-oy*zeta*ozi, x*zeta*ozi, x*(1-y*q*ozi**2)),
        pv(-y*zeta*ozi, -x*zeta*ozi, x*y*q*ozi**2),
        pv(y*zeta*ozi, -ox*zeta*ozi, y*(1-x*q*ozi**2))]
