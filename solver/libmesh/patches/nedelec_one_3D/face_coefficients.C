// Anonymous namespace for local helper functions
namespace
{
// The coefficients of the second order shape function for degree of
// freedom k on side s of elem, in terms of the face shape functions
// written below for the side node order of s.
//
// Those face shape functions are dual to moments of the tangential
// trace over the face, parametrized from its first vertex along the
// second and last ones: on a triangle against the vectors from the
// first vertex to the other two, on a quadrilateral against each of
// the two vectors times the linear Lagrange basis along it.  The
// degrees of freedom of the element are instead those moments for the
// parametrization from the face's oriented vertices, which every
// element sharing the face agrees on (Elem::oriented_face_vertices()):
// on a triangle they mix the face shape functions with coefficients
// -1, 0 or 1, on a quadrilateral they permute them and change their
// sign.
std::array<Real, 4> nedelec_one_face_coefficients(const Elem & elem,
                                                  const unsigned int s,
                                                  const unsigned int k)
{
  // The corners of the unit triangle or square that the side node
  // order maps to the face
  static const int tri_corner[3][2]  = {{0, 0}, {1, 0}, {0, 1}};
  static const int quad_corner[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};

  const std::vector<unsigned int> vertices = elem.oriented_face_vertices(s);
  const unsigned int n_vertices = vertices.size();
  libmesh_assert_less (k, 2*n_vertices - 4);

  // The corners of the oriented vertices
  int x[4][2];
  for (unsigned int v = 0; v != n_vertices; ++v)
    {
      unsigned int m = 0;
      while (elem.local_side_node(s, m) != vertices[v])
        ++m;
      libmesh_assert_less (m, n_vertices);

      for (unsigned int d = 0; d != 2; ++d)
        x[v][d] = (n_vertices == 3) ? tri_corner[m][d] : quad_corner[m][d];
    }

  std::array<Real, 4> c {};

  if (n_vertices == 3)
    {
      // The moments against x1-x0 and x2-x0 are those of the side
      // node order combined by the rows of [[a, b], [e, d]], so the
      // shape functions combine by its inverse transpose.  Its
      // determinant is 1 or -1, its own inverse.
      const int a = x[1][0] - x[0][0], b = x[1][1] - x[0][1],
                e = x[2][0] - x[0][0], d = x[2][1] - x[0][1];
      const int det = a*d - b*e;
      libmesh_assert(det == 1 || det == -1);

      if (k == 0)
        {
          c[0] =  d*det;
          c[1] = -e*det;
        }
      else
        {
          c[0] = -b*det;
          c[1] =  a*det;
        }

      return c;
    }

  // On a quadrilateral, degrees of freedom 0, 1 are the moments
  // against T (1-S) and T S, with T the vector from the least vertex
  // to its lesser neighbor and S the coordinate along T; 2, 3 likewise
  // along the vector to its other neighbor.  T is plus or minus the
  // vector along s (face shape functions 0, 1) or along t (face shape
  // functions 2, 3) of the side node order, and S is s (or t) or one
  // minus it.
  const unsigned int v = (k < 2) ? 1 : 3;
  const int ds = x[v][0] - x[0][0],
            dt = x[v][1] - x[0][1];
  const unsigned int dir = (ds != 0) ? 0 : 1;
  const int sgn = ds + dt;
  libmesh_assert(sgn == 1 || sgn == -1);

  const unsigned int lagrange = k % 2;
  c[2*dir + ((sgn > 0) ? lagrange : 1 - lagrange)] = sgn;

  return c;
}
} // anonymous namespace


