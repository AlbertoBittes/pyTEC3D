      // quadratic Nedelec (first kind) shape functions
    case SECOND:
      {
        switch (elem->type())
          {
          case TET14:
            {
              switch(ii)
                {
                case 0:
                  return sign * RealGradient( 2.0*(4.0*eta*eta+4.0*eta*xi+8.0*eta*zeta-6.0*eta+4.0*xi*zeta-3.0*xi+4.0*zeta*zeta-6.0*zeta+2.0), -2.0*xi*(4.0*eta+4.0*xi+4.0*zeta-3.0), -2.0*xi*(4.0*eta+4.0*xi+4.0*zeta-3.0) );
                case 1:
                  return sign * RealGradient( -2.0*(4.0*eta*xi-eta+4.0*xi*zeta-3.0*xi-zeta+1.0), 4.0*xi*(2.0*xi-1.0), 4.0*xi*(2.0*xi-1.0) );
                case 2:
                  return sign * RealGradient( -2.0*eta*(4.0*xi-1.0), 4.0*xi*(2.0*xi-1.0), 0.0 );
                case 3:
                  return sign * RealGradient( -4.0*eta*(2.0*eta-1.0), 2.0*xi*(4.0*eta-1.0), 0.0 );
                case 4:
                  return sign * RealGradient( -2.0*eta*(4.0*eta+4.0*xi+4.0*zeta-3.0), 2.0*(4.0*eta*xi+4.0*eta*zeta-3.0*eta+4.0*xi*xi+8.0*xi*zeta-6.0*xi+4.0*zeta*zeta-6.0*zeta+2.0), -2.0*eta*(4.0*eta+4.0*xi+4.0*zeta-3.0) );
                case 5:
                  return sign * RealGradient( 4.0*eta*(2.0*eta-1.0), -2.0*(4.0*eta*xi+4.0*eta*zeta-3.0*eta-xi-zeta+1.0), 4.0*eta*(2.0*eta-1.0) );
                case 6:
                  return sign * RealGradient( -2.0*zeta*(4.0*eta+4.0*xi+4.0*zeta-3.0), -2.0*zeta*(4.0*eta+4.0*xi+4.0*zeta-3.0), 2.0*(4.0*eta*eta+8.0*eta*xi+4.0*eta*zeta-6.0*eta+4.0*xi*xi+4.0*xi*zeta-6.0*xi-3.0*zeta+2.0) );
                case 7:
                  return sign * RealGradient( 4.0*zeta*(2.0*zeta-1.0), 4.0*zeta*(2.0*zeta-1.0), -2.0*(4.0*eta*zeta-eta+4.0*xi*zeta-xi-3.0*zeta+1.0) );
                case 8:
                  return sign * RealGradient( -2.0*zeta*(4.0*xi-1.0), 0.0, 4.0*xi*(2.0*xi-1.0) );
                case 9:
                  return sign * RealGradient( -4.0*zeta*(2.0*zeta-1.0), 0.0, 2.0*xi*(4.0*zeta-1.0) );
                case 10:
                  return sign * RealGradient( 0.0, -2.0*zeta*(4.0*eta-1.0), 4.0*eta*(2.0*eta-1.0) );
                case 11:
                  return sign * RealGradient( 0.0, -4.0*zeta*(2.0*zeta-1.0), 2.0*eta*(4.0*zeta-1.0) );
                case 12:
                case 13:
                  {
                    const std::array<Real, 4> c =
                      nedelec_one_face_coefficients(*elem, 0, i - 12);
                    return c[0] * RealGradient( 8.0*eta*(eta+2.0*xi+zeta-1.0), -8.0*xi*(eta+2.0*xi+2.0*zeta-2.0), 8.0*eta*xi ) +
                           c[1] * RealGradient( -8.0*eta*(2.0*eta+xi+2.0*zeta-2.0), 8.0*xi*(2.0*eta+xi+zeta-1.0), 8.0*eta*xi );
                  }
                case 14:
                case 15:
                  {
                    const std::array<Real, 4> c =
                      nedelec_one_face_coefficients(*elem, 1, i - 14);
                    return c[0] * RealGradient( -8.0*zeta*(2.0*eta+xi+2.0*zeta-2.0), 8.0*xi*zeta, 8.0*xi*(eta+xi+2.0*zeta-1.0) ) +
                           c[1] * RealGradient( 8.0*zeta*(eta+2.0*xi+zeta-1.0), 8.0*xi*zeta, -8.0*xi*(2.0*eta+2.0*xi+zeta-2.0) );
                  }
                case 16:
                case 17:
                  {
                    const std::array<Real, 4> c =
                      nedelec_one_face_coefficients(*elem, 2, i - 16);
                    return c[0] * RealGradient( -8.0*eta*zeta, 16.0*xi*zeta, -8.0*eta*xi ) +
                           c[1] * RealGradient( -8.0*eta*zeta, -8.0*xi*zeta, 16.0*eta*xi );
                  }
                case 18:
                case 19:
                  {
                    const std::array<Real, 4> c =
                      nedelec_one_face_coefficients(*elem, 3, i - 18);
                    return c[0] * RealGradient( -16.0*eta*zeta, -8.0*zeta*(eta-xi-zeta+1.0), 8.0*eta*(eta+xi-zeta-1.0) ) +
                           c[1] * RealGradient( 8.0*eta*zeta, 8.0*zeta*(2.0*eta+xi+zeta-1.0), -8.0*eta*(2.0*eta+2.0*xi+zeta-2.0) );
                  }

                default:
                  libmesh_error_msg("Invalid i = " << i);
                } // switch(ii)
            }

          case HEX27:
            {
              // Even with a loose inverse_map tolerance we ought to
              // be nearly on the element interior in master
              // coordinates
              libmesh_assert_less_equal ( std::fabs(xi),   1.0+10*TOLERANCE );
              libmesh_assert_less_equal ( std::fabs(eta),  1.0+10*TOLERANCE );
              libmesh_assert_less_equal ( std::fabs(zeta), 1.0+10*TOLERANCE );
              switch(ii)
                {
                case 0:
                  return sign * RealGradient( -0.03125*(3.0*eta+1.0)*(3.0*zeta+1.0)*(eta-1.0)*(zeta-1.0)*(3.0*xi-1.0), 0.0, 0.0 );
                case 1:
                  return sign * RealGradient( 0.03125*(3.0*eta+1.0)*(3.0*xi+1.0)*(3.0*zeta+1.0)*(eta-1.0)*(zeta-1.0), 0.0, 0.0 );
                case 2:
                  return sign * RealGradient( 0.0, -0.03125*(xi+1.0)*(3.0*zeta+1.0)*(zeta-1.0)*(3.0*eta-1.0)*(3.0*xi-1.0), 0.0 );
                case 3:
                  return sign * RealGradient( 0.0, 0.03125*(xi+1.0)*(3.0*eta+1.0)*(3.0*zeta+1.0)*(zeta-1.0)*(3.0*xi-1.0), 0.0 );
                case 4:
                  return sign * RealGradient( -0.03125*(eta+1.0)*(3.0*xi+1.0)*(3.0*zeta+1.0)*(zeta-1.0)*(3.0*eta-1.0), 0.0, 0.0 );
                case 5:
                  return sign * RealGradient( 0.03125*(eta+1.0)*(3.0*zeta+1.0)*(zeta-1.0)*(3.0*eta-1.0)*(3.0*xi-1.0), 0.0, 0.0 );
                case 6:
                  return sign * RealGradient( 0.0, -0.03125*(3.0*xi+1.0)*(3.0*zeta+1.0)*(xi-1.0)*(zeta-1.0)*(3.0*eta-1.0), 0.0 );
                case 7:
                  return sign * RealGradient( 0.0, 0.03125*(3.0*eta+1.0)*(3.0*xi+1.0)*(3.0*zeta+1.0)*(xi-1.0)*(zeta-1.0), 0.0 );
                case 8:
                  return sign * RealGradient( 0.0, 0.0, -0.03125*(3.0*eta+1.0)*(3.0*xi+1.0)*(eta-1.0)*(xi-1.0)*(3.0*zeta-1.0) );
                case 9:
                  return sign * RealGradient( 0.0, 0.0, 0.03125*(3.0*eta+1.0)*(3.0*xi+1.0)*(3.0*zeta+1.0)*(eta-1.0)*(xi-1.0) );
                case 10:
                  return sign * RealGradient( 0.0, 0.0, -0.03125*(xi+1.0)*(3.0*eta+1.0)*(eta-1.0)*(3.0*xi-1.0)*(3.0*zeta-1.0) );
                case 11:
                  return sign * RealGradient( 0.0, 0.0, 0.03125*(xi+1.0)*(3.0*eta+1.0)*(3.0*zeta+1.0)*(eta-1.0)*(3.0*xi-1.0) );
                case 12:
                  return sign * RealGradient( 0.0, 0.0, -0.03125*(eta+1.0)*(xi+1.0)*(3.0*eta-1.0)*(3.0*xi-1.0)*(3.0*zeta-1.0) );
                case 13:
                  return sign * RealGradient( 0.0, 0.0, 0.03125*(eta+1.0)*(xi+1.0)*(3.0*zeta+1.0)*(3.0*eta-1.0)*(3.0*xi-1.0) );
                case 14:
                  return sign * RealGradient( 0.0, 0.0, -0.03125*(eta+1.0)*(3.0*xi+1.0)*(xi-1.0)*(3.0*eta-1.0)*(3.0*zeta-1.0) );
                case 15:
                  return sign * RealGradient( 0.0, 0.0, 0.03125*(eta+1.0)*(3.0*xi+1.0)*(3.0*zeta+1.0)*(xi-1.0)*(3.0*eta-1.0) );
                case 16:
                  return sign * RealGradient( -0.03125*(zeta+1.0)*(3.0*eta+1.0)*(eta-1.0)*(3.0*xi-1.0)*(3.0*zeta-1.0), 0.0, 0.0 );
                case 17:
                  return sign * RealGradient( 0.03125*(zeta+1.0)*(3.0*eta+1.0)*(3.0*xi+1.0)*(eta-1.0)*(3.0*zeta-1.0), 0.0, 0.0 );
                case 18:
                  return sign * RealGradient( 0.0, -0.03125*(xi+1.0)*(zeta+1.0)*(3.0*eta-1.0)*(3.0*xi-1.0)*(3.0*zeta-1.0), 0.0 );
                case 19:
                  return sign * RealGradient( 0.0, 0.03125*(xi+1.0)*(zeta+1.0)*(3.0*eta+1.0)*(3.0*xi-1.0)*(3.0*zeta-1.0), 0.0 );
                case 20:
                  return sign * RealGradient( -0.03125*(eta+1.0)*(zeta+1.0)*(3.0*xi+1.0)*(3.0*eta-1.0)*(3.0*zeta-1.0), 0.0, 0.0 );
                case 21:
                  return sign * RealGradient( 0.03125*(eta+1.0)*(zeta+1.0)*(3.0*eta-1.0)*(3.0*xi-1.0)*(3.0*zeta-1.0), 0.0, 0.0 );
                case 22:
                  return sign * RealGradient( 0.0, -0.03125*(zeta+1.0)*(3.0*xi+1.0)*(xi-1.0)*(3.0*eta-1.0)*(3.0*zeta-1.0), 0.0 );
                case 23:
                  return sign * RealGradient( 0.0, 0.03125*(zeta+1.0)*(3.0*eta+1.0)*(3.0*xi+1.0)*(xi-1.0)*(3.0*zeta-1.0), 0.0 );
                case 24:
                case 25:
                case 26:
                case 27:
                  {
                    const std::array<Real, 4> c =
                      nedelec_one_face_coefficients(*elem, 0, i - 24);
                    return c[0] * RealGradient( 0.0, 0.1875*(xi+1.0)*(3.0*zeta+1.0)*(xi-1.0)*(zeta-1.0)*(3.0*eta-1.0), 0.0 ) +
                           c[1] * RealGradient( 0.0, -0.1875*(xi+1.0)*(3.0*eta+1.0)*(3.0*zeta+1.0)*(xi-1.0)*(zeta-1.0), 0.0 ) +
                           c[2] * RealGradient( 0.1875*(eta+1.0)*(3.0*zeta+1.0)*(eta-1.0)*(zeta-1.0)*(3.0*xi-1.0), 0.0, 0.0 ) +
                           c[3] * RealGradient( -0.1875*(eta+1.0)*(3.0*xi+1.0)*(3.0*zeta+1.0)*(eta-1.0)*(zeta-1.0), 0.0, 0.0 );
                  }
                case 28:
                case 29:
                case 30:
                case 31:
                  {
                    const std::array<Real, 4> c =
                      nedelec_one_face_coefficients(*elem, 1, i - 28);
                    return c[0] * RealGradient( 0.1875*(zeta+1.0)*(3.0*eta+1.0)*(eta-1.0)*(zeta-1.0)*(3.0*xi-1.0), 0.0, 0.0 ) +
                           c[1] * RealGradient( -0.1875*(zeta+1.0)*(3.0*eta+1.0)*(3.0*xi+1.0)*(eta-1.0)*(zeta-1.0), 0.0, 0.0 ) +
                           c[2] * RealGradient( 0.0, 0.0, 0.1875*(xi+1.0)*(3.0*eta+1.0)*(eta-1.0)*(xi-1.0)*(3.0*zeta-1.0) ) +
                           c[3] * RealGradient( 0.0, 0.0, -0.1875*(xi+1.0)*(3.0*eta+1.0)*(3.0*zeta+1.0)*(eta-1.0)*(xi-1.0) );
                  }
                case 32:
                case 33:
                case 34:
                case 35:
                  {
                    const std::array<Real, 4> c =
                      nedelec_one_face_coefficients(*elem, 2, i - 32);
                    return c[0] * RealGradient( 0.0, 0.1875*(xi+1.0)*(zeta+1.0)*(zeta-1.0)*(3.0*eta-1.0)*(3.0*xi-1.0), 0.0 ) +
                           c[1] * RealGradient( 0.0, -0.1875*(xi+1.0)*(zeta+1.0)*(3.0*eta+1.0)*(zeta-1.0)*(3.0*xi-1.0), 0.0 ) +
                           c[2] * RealGradient( 0.0, 0.0, 0.1875*(eta+1.0)*(xi+1.0)*(eta-1.0)*(3.0*xi-1.0)*(3.0*zeta-1.0) ) +
                           c[3] * RealGradient( 0.0, 0.0, -0.1875*(eta+1.0)*(xi+1.0)*(3.0*zeta+1.0)*(eta-1.0)*(3.0*xi-1.0) );
                  }
                case 36:
                case 37:
                case 38:
                case 39:
                  {
                    const std::array<Real, 4> c =
                      nedelec_one_face_coefficients(*elem, 3, i - 36);
                    return c[0] * RealGradient( 0.1875*(eta+1.0)*(zeta+1.0)*(3.0*xi+1.0)*(zeta-1.0)*(3.0*eta-1.0), 0.0, 0.0 ) +
                           c[1] * RealGradient( -0.1875*(eta+1.0)*(zeta+1.0)*(zeta-1.0)*(3.0*eta-1.0)*(3.0*xi-1.0), 0.0, 0.0 ) +
                           c[2] * RealGradient( 0.0, 0.0, 0.1875*(eta+1.0)*(xi+1.0)*(xi-1.0)*(3.0*eta-1.0)*(3.0*zeta-1.0) ) +
                           c[3] * RealGradient( 0.0, 0.0, -0.1875*(eta+1.0)*(xi+1.0)*(3.0*zeta+1.0)*(xi-1.0)*(3.0*eta-1.0) );
                  }
                case 40:
                case 41:
                case 42:
                case 43:
                  {
                    const std::array<Real, 4> c =
                      nedelec_one_face_coefficients(*elem, 4, i - 40);
                    return c[0] * RealGradient( 0.0, 0.1875*(zeta+1.0)*(3.0*eta+1.0)*(3.0*xi+1.0)*(xi-1.0)*(zeta-1.0), 0.0 ) +
                           c[1] * RealGradient( 0.0, -0.1875*(zeta+1.0)*(3.0*xi+1.0)*(xi-1.0)*(zeta-1.0)*(3.0*eta-1.0), 0.0 ) +
                           c[2] * RealGradient( 0.0, 0.0, 0.1875*(eta+1.0)*(3.0*xi+1.0)*(eta-1.0)*(xi-1.0)*(3.0*zeta-1.0) ) +
                           c[3] * RealGradient( 0.0, 0.0, -0.1875*(eta+1.0)*(3.0*xi+1.0)*(3.0*zeta+1.0)*(eta-1.0)*(xi-1.0) );
                  }
                case 44:
                case 45:
                case 46:
                case 47:
                  {
                    const std::array<Real, 4> c =
                      nedelec_one_face_coefficients(*elem, 5, i - 44);
                    return c[0] * RealGradient( 0.1875*(eta+1.0)*(zeta+1.0)*(eta-1.0)*(3.0*xi-1.0)*(3.0*zeta-1.0), 0.0, 0.0 ) +
                           c[1] * RealGradient( -0.1875*(eta+1.0)*(zeta+1.0)*(3.0*xi+1.0)*(eta-1.0)*(3.0*zeta-1.0), 0.0, 0.0 ) +
                           c[2] * RealGradient( 0.0, 0.1875*(xi+1.0)*(zeta+1.0)*(xi-1.0)*(3.0*eta-1.0)*(3.0*zeta-1.0), 0.0 ) +
                           c[3] * RealGradient( 0.0, -0.1875*(xi+1.0)*(zeta+1.0)*(3.0*eta+1.0)*(xi-1.0)*(3.0*zeta-1.0), 0.0 );
                  }
                case 48:
                  return RealGradient( -0.28125*(eta+1.0)*(zeta+1.0)*(eta-1.0)*(zeta-1.0)*(3.0*xi-1.0), 0.0, 0.0 );
                case 49:
                  return RealGradient( 0.28125*(eta+1.0)*(zeta+1.0)*(3.0*xi+1.0)*(eta-1.0)*(zeta-1.0), 0.0, 0.0 );
                case 50:
                  return RealGradient( 0.0, -0.28125*(xi+1.0)*(zeta+1.0)*(xi-1.0)*(zeta-1.0)*(3.0*eta-1.0), 0.0 );
                case 51:
                  return RealGradient( 0.0, 0.28125*(xi+1.0)*(zeta+1.0)*(3.0*eta+1.0)*(xi-1.0)*(zeta-1.0), 0.0 );
                case 52:
                  return RealGradient( 0.0, 0.0, -0.28125*(eta+1.0)*(xi+1.0)*(eta-1.0)*(xi-1.0)*(3.0*zeta-1.0) );
                case 53:
                  return RealGradient( 0.0, 0.0, 0.28125*(eta+1.0)*(xi+1.0)*(3.0*zeta+1.0)*(eta-1.0)*(xi-1.0) );

                default:
                  libmesh_error_msg("Invalid i = " << i);
                } // switch(ii)
            }

          case PRISM20:
          case PRISM21:
            {
              // Even with a loose inverse_map tolerance we ought to
              // be nearly on the element interior in master
              // coordinates
              libmesh_assert_greater_equal ( xi,  -10*TOLERANCE );
              libmesh_assert_greater_equal ( eta, -10*TOLERANCE );
              libmesh_assert_less_equal ( xi+eta, 1.0+10*TOLERANCE );
              libmesh_assert_less_equal ( std::fabs(zeta), 1.0+10*TOLERANCE );
              switch(ii)
                {
                case 0:
                  return sign * RealGradient( 0.5*(3.0*zeta+1.0)*(zeta-1.0)*(4.0*eta*eta+4.0*eta*xi-6.0*eta-3.0*xi+2.0), -0.5*xi*(3.0*zeta+1.0)*(zeta-1.0)*(4.0*eta+4.0*xi-3.0), 0.0 );
                case 1:
                  return sign * RealGradient( -0.5*(3.0*zeta+1.0)*(zeta-1.0)*(4.0*eta*xi-eta-3.0*xi+1.0), xi*(3.0*zeta+1.0)*(zeta-1.0)*(2.0*xi-1.0), 0.0 );
                case 2:
                  return sign * RealGradient( -0.5*eta*(3.0*zeta+1.0)*(zeta-1.0)*(4.0*xi-1.0), xi*(3.0*zeta+1.0)*(zeta-1.0)*(2.0*xi-1.0), 0.0 );
                case 3:
                  return sign * RealGradient( -eta*(3.0*zeta+1.0)*(zeta-1.0)*(2.0*eta-1.0), 0.5*xi*(3.0*zeta+1.0)*(zeta-1.0)*(4.0*eta-1.0), 0.0 );
                case 4:
                  return sign * RealGradient( -0.5*eta*(3.0*zeta+1.0)*(zeta-1.0)*(4.0*eta+4.0*xi-3.0), 0.5*(3.0*zeta+1.0)*(zeta-1.0)*(4.0*eta*xi-3.0*eta+4.0*xi*xi-6.0*xi+2.0), 0.0 );
                case 5:
                  return sign * RealGradient( eta*(3.0*zeta+1.0)*(zeta-1.0)*(2.0*eta-1.0), -0.5*(3.0*zeta+1.0)*(zeta-1.0)*(4.0*eta*xi-3.0*eta-xi+1.0), 0.0 );
                case 6:
                  return sign * RealGradient( 0.0, 0.0, -0.5*(3.0*zeta-1.0)*(eta+xi-1.0)*(3.0*eta+3.0*xi-1.0) );
                case 7:
                  return sign * RealGradient( 0.0, 0.0, 0.5*(3.0*zeta+1.0)*(eta+xi-1.0)*(3.0*eta+3.0*xi-1.0) );
                case 8:
                  return sign * RealGradient( 0.0, 0.0, -0.5*xi*(3.0*zeta-1.0)*(3.0*xi-2.0) );
                case 9:
                  return sign * RealGradient( 0.0, 0.0, 0.5*xi*(3.0*zeta+1.0)*(3.0*xi-2.0) );
                case 10:
                  return sign * RealGradient( 0.0, 0.0, -0.5*eta*(3.0*zeta-1.0)*(3.0*eta-2.0) );
                case 11:
                  return sign * RealGradient( 0.0, 0.0, 0.5*eta*(3.0*zeta+1.0)*(3.0*eta-2.0) );
                case 12:
                  return sign * RealGradient( 0.5*(zeta+1.0)*(3.0*zeta-1.0)*(4.0*eta*eta+4.0*eta*xi-6.0*eta-3.0*xi+2.0), -0.5*xi*(zeta+1.0)*(3.0*zeta-1.0)*(4.0*eta+4.0*xi-3.0), 0.0 );
                case 13:
                  return sign * RealGradient( -0.5*(zeta+1.0)*(3.0*zeta-1.0)*(4.0*eta*xi-eta-3.0*xi+1.0), xi*(zeta+1.0)*(2.0*xi-1.0)*(3.0*zeta-1.0), 0.0 );
                case 14:
                  return sign * RealGradient( -0.5*eta*(zeta+1.0)*(3.0*zeta-1.0)*(4.0*xi-1.0), xi*(zeta+1.0)*(2.0*xi-1.0)*(3.0*zeta-1.0), 0.0 );
                case 15:
                  return sign * RealGradient( -eta*(zeta+1.0)*(2.0*eta-1.0)*(3.0*zeta-1.0), 0.5*xi*(zeta+1.0)*(3.0*zeta-1.0)*(4.0*eta-1.0), 0.0 );
                case 16:
                  return sign * RealGradient( -0.5*eta*(zeta+1.0)*(3.0*zeta-1.0)*(4.0*eta+4.0*xi-3.0), 0.5*(zeta+1.0)*(3.0*zeta-1.0)*(4.0*eta*xi-3.0*eta+4.0*xi*xi-6.0*xi+2.0), 0.0 );
                case 17:
                  return sign * RealGradient( eta*(zeta+1.0)*(2.0*eta-1.0)*(3.0*zeta-1.0), -0.5*(zeta+1.0)*(3.0*zeta-1.0)*(4.0*eta*xi-3.0*eta-xi+1.0), 0.0 );
                case 18:
                case 19:
                case 20:
                case 21:
                  {
                    const std::array<Real, 4> c =
                      nedelec_one_face_coefficients(*elem, 1, i - 18);
                    return c[0] * RealGradient( -3.0*(zeta+1.0)*(zeta-1.0)*(4.0*eta*eta+4.0*eta*xi-6.0*eta-3.0*xi+2.0), 3.0*xi*(zeta+1.0)*(zeta-1.0)*(4.0*eta+4.0*xi-3.0), 0.0 ) +
                           c[1] * RealGradient( 3.0*(zeta+1.0)*(zeta-1.0)*(4.0*eta*xi-eta-3.0*xi+1.0), -6.0*xi*(zeta+1.0)*(zeta-1.0)*(2.0*xi-1.0), 0.0 ) +
                           c[2] * RealGradient( 0.0, 0.0, 3.0*xi*(3.0*zeta-1.0)*(eta+xi-1.0) ) +
                           c[3] * RealGradient( 0.0, 0.0, -3.0*xi*(3.0*zeta+1.0)*(eta+xi-1.0) );
                  }
                case 22:
                case 23:
                case 24:
                case 25:
                  {
                    const std::array<Real, 4> c =
                      nedelec_one_face_coefficients(*elem, 2, i - 22);
                    return c[0] * RealGradient( 3.0*eta*(zeta+1.0)*(zeta-1.0)*(4.0*xi-1.0), -6.0*xi*(zeta+1.0)*(zeta-1.0)*(2.0*xi-1.0), 0.0 ) +
                           c[1] * RealGradient( 6.0*eta*(zeta+1.0)*(zeta-1.0)*(2.0*eta-1.0), -3.0*xi*(zeta+1.0)*(zeta-1.0)*(4.0*eta-1.0), 0.0 ) +
                           c[2] * RealGradient( 0.0, 0.0, -3.0*eta*xi*(3.0*zeta-1.0) ) +
                           c[3] * RealGradient( 0.0, 0.0, 3.0*eta*xi*(3.0*zeta+1.0) );
                  }
                case 26:
                case 27:
                case 28:
                case 29:
                  {
                    const std::array<Real, 4> c =
                      nedelec_one_face_coefficients(*elem, 3, i - 26);
                    return c[0] * RealGradient( 6.0*eta*(zeta+1.0)*(zeta-1.0)*(2.0*eta-1.0), -3.0*(zeta+1.0)*(zeta-1.0)*(4.0*eta*xi-3.0*eta-xi+1.0), 0.0 ) +
                           c[1] * RealGradient( -3.0*eta*(zeta+1.0)*(zeta-1.0)*(4.0*eta+4.0*xi-3.0), 3.0*(zeta+1.0)*(zeta-1.0)*(4.0*eta*xi-3.0*eta+4.0*xi*xi-6.0*xi+2.0), 0.0 ) +
                           c[2] * RealGradient( 0.0, 0.0, 3.0*eta*(3.0*zeta-1.0)*(eta+xi-1.0) ) +
                           c[3] * RealGradient( 0.0, 0.0, -3.0*eta*(3.0*zeta+1.0)*(eta+xi-1.0) );
                  }
                case 30:
                case 31:
                  {
                    const std::array<Real, 4> c =
                      nedelec_one_face_coefficients(*elem, 0, i - 30);
                    return c[0] * RealGradient( 2.0*eta*(3.0*zeta+1.0)*(zeta-1.0)*(eta+2.0*xi-1.0), -2.0*xi*(3.0*zeta+1.0)*(zeta-1.0)*(eta+2.0*xi-2.0), 0.0 ) +
                           c[1] * RealGradient( -2.0*eta*(3.0*zeta+1.0)*(zeta-1.0)*(2.0*eta+xi-2.0), 2.0*xi*(3.0*zeta+1.0)*(zeta-1.0)*(2.0*eta+xi-1.0), 0.0 );
                  }
                case 32:
                case 33:
                  {
                    const std::array<Real, 4> c =
                      nedelec_one_face_coefficients(*elem, 4, i - 32);
                    return c[0] * RealGradient( -2.0*eta*(zeta+1.0)*(3.0*zeta-1.0)*(2.0*eta+xi-2.0), 2.0*xi*(zeta+1.0)*(3.0*zeta-1.0)*(2.0*eta+xi-1.0), 0.0 ) +
                           c[1] * RealGradient( 2.0*eta*(zeta+1.0)*(3.0*zeta-1.0)*(eta+2.0*xi-1.0), -2.0*xi*(zeta+1.0)*(3.0*zeta-1.0)*(eta+2.0*xi-2.0), 0.0 );
                  }
                case 34:
                  return RealGradient( 6.0*eta*(zeta+1.0)*(zeta-1.0)*(2.0*eta+xi-2.0), -6.0*xi*(zeta+1.0)*(zeta-1.0)*(2.0*eta+xi-1.0), 0.0 );
                case 35:
                  return RealGradient( -6.0*eta*(zeta+1.0)*(zeta-1.0)*(eta+2.0*xi-1.0), 6.0*xi*(zeta+1.0)*(zeta-1.0)*(eta+2.0*xi-2.0), 0.0 );

                default:
                  libmesh_error_msg("Invalid i = " << i);
                } // switch(ii)
            }

          case PYRAMID18:
            {
              // Even with a loose inverse_map tolerance we ought to
              // be nearly on the element interior in master
              // coordinates
              libmesh_assert_greater_equal ( zeta, -10*TOLERANCE );
              libmesh_assert_less_equal ( zeta, 1.0+10*TOLERANCE );
              libmesh_assert_less_equal ( std::fabs(xi),  1.0-zeta+10*TOLERANCE );
              libmesh_assert_less_equal ( std::fabs(eta), 1.0-zeta+10*TOLERANCE );

              // The shape functions are rational in den = 1-zeta. At the apex,
              // where den vanishes, they only have limits along each direction:
              // take those along the pyramid axis.
              if (1.0 - zeta < 1.e-6)
                {
                  switch(ii)
                    {
                    case 8:
                      return sign * RealGradient( -0.5, -0.5, 0.125 );
                    case 9:
                      return sign * RealGradient( 1.0, 1.0, -0.25 );
                    case 10:
                      return sign * RealGradient( 0.5, -0.5, 0.125 );
                    case 11:
                      return sign * RealGradient( -1.0, 1.0, -0.25 );
                    case 12:
                      return sign * RealGradient( 0.5, 0.5, 0.125 );
                    case 13:
                      return sign * RealGradient( -1.0, -1.0, -0.25 );
                    case 14:
                      return sign * RealGradient( -0.5, 0.5, 0.125 );
                    case 15:
                      return sign * RealGradient( 1.0, -1.0, -0.25 );
                    case 16:
                    case 17:
                    case 18:
                    case 19:
                      {
                        const std::array<Real, 4> c =
                          nedelec_one_face_coefficients(*elem, 4, i - 16);
                        return c[0] * RealGradient( 0.0, 0.0, -1.5 ) +
                               c[1] * RealGradient( 0.0, 0.0, 1.5 ) +
                               c[2] * RealGradient( 0.0, 0.0, -1.5 ) +
                               c[3] * RealGradient( 0.0, 0.0, 1.5 );
                      }
                    case 30:
                      return RealGradient( 0.0, 0.0, -3.75 );
                    case 31:
                      return RealGradient( 0.0, 0.0, -3.75 );
                    case 32:
                      return RealGradient( 0.0, 0.0, -3.75 );
                    case 33:
                      return RealGradient( 0.0, 0.0, 18.75 );
                    case 0:
                    case 1:
                    case 2:
                    case 3:
                    case 4:
                    case 5:
                    case 6:
                    case 7:
                    case 20:
                    case 21:
                    case 22:
                    case 23:
                    case 24:
                    case 25:
                    case 26:
                    case 27:
                    case 28:
                    case 29:
                      return RealGradient();

                    default:
                      libmesh_error_msg("Invalid i = " << i);
                    } // switch(ii)
                }

              const Real den = 1.0 - zeta;

              switch(ii)
                {
                case 0:
                  return sign * RealGradient( 0.125*(3.0*eta-zeta+1.0)*(eta+zeta-1.0)*(4.0*xi*zeta-3.0*xi+4.0*zeta*zeta-5.0*zeta+1.0)/(den*den), 0.0, -0.125*(3.0*eta-zeta+1.0)*(eta+zeta-1.0)*(4.0*xi*xi+4.0*xi*zeta-xi+zeta-1.0)/(den*den) );
                case 1:
                  return sign * RealGradient( -0.125*(3.0*eta-zeta+1.0)*(eta+zeta-1.0)*(4.0*xi*zeta-3.0*xi-4.0*zeta*zeta+5.0*zeta-1.0)/(den*den), 0.0, 0.125*(3.0*eta-zeta+1.0)*(eta+zeta-1.0)*(4.0*xi*xi-4.0*xi*zeta+xi+zeta-1.0)/(den*den) );
                case 2:
                  return sign * RealGradient( 0.0, 0.125*(xi-zeta+1.0)*(3.0*xi+zeta-1.0)*(4.0*eta*zeta-3.0*eta+4.0*zeta*zeta-5.0*zeta+1.0)/(den*den), -0.125*(xi-zeta+1.0)*(3.0*xi+zeta-1.0)*(4.0*eta*eta+4.0*eta*zeta-eta+zeta-1.0)/(den*den) );
                case 3:
                  return sign * RealGradient( 0.0, -0.125*(xi-zeta+1.0)*(3.0*xi+zeta-1.0)*(4.0*eta*zeta-3.0*eta-4.0*zeta*zeta+5.0*zeta-1.0)/(den*den), 0.125*(xi-zeta+1.0)*(3.0*xi+zeta-1.0)*(4.0*eta*eta-4.0*eta*zeta+eta+zeta-1.0)/(den*den) );
                case 4:
                  return sign * RealGradient( 0.125*(eta-zeta+1.0)*(3.0*eta+zeta-1.0)*(4.0*xi*zeta-3.0*xi-4.0*zeta*zeta+5.0*zeta-1.0)/(den*den), 0.0, -0.125*(eta-zeta+1.0)*(3.0*eta+zeta-1.0)*(4.0*xi*xi-4.0*xi*zeta+xi+zeta-1.0)/(den*den) );
                case 5:
                  return sign * RealGradient( -0.125*(eta-zeta+1.0)*(3.0*eta+zeta-1.0)*(4.0*xi*zeta-3.0*xi+4.0*zeta*zeta-5.0*zeta+1.0)/(den*den), 0.0, 0.125*(eta-zeta+1.0)*(3.0*eta+zeta-1.0)*(4.0*xi*xi+4.0*xi*zeta-xi+zeta-1.0)/(den*den) );
                case 6:
                  return sign * RealGradient( 0.0, 0.125*(3.0*xi-zeta+1.0)*(xi+zeta-1.0)*(4.0*eta*zeta-3.0*eta+4.0*zeta*zeta-5.0*zeta+1.0)/(den*den), -0.125*(3.0*xi-zeta+1.0)*(xi+zeta-1.0)*(4.0*eta*eta+4.0*eta*zeta-eta+zeta-1.0)/(den*den) );
                case 7:
                  return sign * RealGradient( 0.0, -0.125*(3.0*xi-zeta+1.0)*(xi+zeta-1.0)*(4.0*eta*zeta-3.0*eta-4.0*zeta*zeta+5.0*zeta-1.0)/(den*den), 0.125*(3.0*xi-zeta+1.0)*(xi+zeta-1.0)*(4.0*eta*eta-4.0*eta*zeta+eta+zeta-1.0)/(den*den) );
                case 8:
                  return sign * RealGradient( -0.125*zeta*(eta+zeta-1.0)*(3.0*eta*xi*zeta+7.0*eta*xi-2.0*eta*zeta*zeta+4.0*eta*zeta-2.0*eta-11.0*xi*zeta*zeta+12.0*xi*zeta-xi-6.0*zeta*zeta*zeta+14.0*zeta*zeta-10.0*zeta+2.0)/(den*den*den), -0.125*zeta*(xi+zeta-1.0)*(3.0*eta*xi*zeta+7.0*eta*xi-11.0*eta*zeta*zeta+12.0*eta*zeta-eta-2.0*xi*zeta*zeta+4.0*xi*zeta-2.0*xi-6.0*zeta*zeta*zeta+14.0*zeta*zeta-10.0*zeta+2.0)/(den*den*den), 0.125*(11.0*eta*eta+16.0*eta*xi+6.0*eta*zeta+11.0*xi*xi+6.0*xi*zeta+2.0*zeta-1.0)-0.125*(6.0*eta*eta*xi*xi*zeta*zeta+16.0*zeta*eta*eta*xi*xi-7.0*eta*eta*xi*xi-10.0*xi*eta*eta*zeta*zeta*zeta+28.0*xi*eta*eta*zeta*zeta-26.0*xi*zeta*eta*eta+8.0*xi*eta*eta-15.0*eta*eta*zeta*zeta*zeta+40.0*eta*eta*zeta*zeta-35.0*zeta*eta*eta+10.0*eta*eta-10.0*eta*xi*xi*zeta*zeta*zeta+28.0*eta*xi*xi*zeta*zeta-26.0*eta*zeta*xi*xi+8.0*eta*xi*xi-12.0*eta*xi*zeta*zeta*zeta+32.0*eta*xi*zeta*zeta-28.0*eta*xi*zeta+8.0*eta*xi-15.0*xi*xi*zeta*zeta*zeta+40.0*xi*xi*zeta*zeta-35.0*zeta*xi*xi+10.0*xi*xi)/(den*den*den*den) );
                case 9:
                  return sign * RealGradient( 0.25*zeta*(eta+zeta-1.0)*(15.0*eta*xi*zeta-5.0*eta*xi-2.0*eta*zeta*zeta+4.0*eta*zeta-2.0*eta-15.0*xi*zeta*zeta+20.0*xi*zeta-5.0*xi-6.0*zeta*zeta*zeta+14.0*zeta*zeta-10.0*zeta+2.0)/(den*den*den), 0.25*zeta*(xi+zeta-1.0)*(15.0*eta*xi*zeta-5.0*eta*xi-15.0*eta*zeta*zeta+20.0*eta*zeta-5.0*eta-2.0*xi*zeta*zeta+4.0*xi*zeta-2.0*xi-6.0*zeta*zeta*zeta+14.0*zeta*zeta-10.0*zeta+2.0)/(den*den*den), -0.125*(30.0*eta*eta+32.0*eta*xi+12.0*eta*zeta-4.0*eta+30.0*xi*xi+12.0*xi*zeta-4.0*xi+5.0*zeta-3.0)+0.125*(60.0*eta*eta*xi*xi*zeta*zeta-37.0*zeta*eta*eta*xi*xi+7.0*eta*eta*xi*xi-4.0*xi*eta*eta*zeta*zeta*zeta+8.0*xi*eta*eta*zeta*zeta-4.0*xi*zeta*eta*eta-33.0*eta*eta*zeta*zeta*zeta+89.0*eta*eta*zeta*zeta-79.0*zeta*eta*eta+23.0*eta*eta-4.0*eta*xi*xi*zeta*zeta*zeta+8.0*eta*xi*xi*zeta*zeta-4.0*eta*zeta*xi*xi-36.0*eta*xi*zeta*zeta*zeta+100.0*eta*xi*zeta*zeta-92.0*eta*xi*zeta+28.0*eta*xi-33.0*xi*xi*zeta*zeta*zeta+89.0*xi*xi*zeta*zeta-79.0*zeta*xi*xi+23.0*xi*xi)/(den*den*den*den) );
                case 10:
                  return sign * RealGradient( -0.125*zeta*(eta+zeta-1.0)*(3.0*eta*xi*zeta+7.0*eta*xi+2.0*eta*zeta*zeta-4.0*eta*zeta+2.0*eta-11.0*xi*zeta*zeta+12.0*xi*zeta-xi+6.0*zeta*zeta*zeta-14.0*zeta*zeta+10.0*zeta-2.0)/(den*den*den), -0.125*zeta*(xi-zeta+1.0)*(3.0*eta*xi*zeta+7.0*eta*xi+11.0*eta*zeta*zeta-12.0*eta*zeta+eta-2.0*xi*zeta*zeta+4.0*xi*zeta-2.0*xi+6.0*zeta*zeta*zeta-14.0*zeta*zeta+10.0*zeta-2.0)/(den*den*den), 0.125*(11.0*eta*eta-16.0*eta*xi+6.0*eta*zeta+11.0*xi*xi-6.0*xi*zeta+2.0*zeta-1.0)-0.125*(6.0*eta*eta*xi*xi*zeta*zeta+16.0*zeta*eta*eta*xi*xi-7.0*eta*eta*xi*xi+10.0*xi*eta*eta*zeta*zeta*zeta-28.0*xi*eta*eta*zeta*zeta+26.0*xi*zeta*eta*eta-8.0*xi*eta*eta-15.0*eta*eta*zeta*zeta*zeta+40.0*eta*eta*zeta*zeta-35.0*zeta*eta*eta+10.0*eta*eta-10.0*eta*xi*xi*zeta*zeta*zeta+28.0*eta*xi*xi*zeta*zeta-26.0*eta*zeta*xi*xi+8.0*eta*xi*xi+12.0*eta*xi*zeta*zeta*zeta-32.0*eta*xi*zeta*zeta+28.0*eta*xi*zeta-8.0*eta*xi-15.0*xi*xi*zeta*zeta*zeta+40.0*xi*xi*zeta*zeta-35.0*zeta*xi*xi+10.0*xi*xi)/(den*den*den*den) );
                case 11:
                  return sign * RealGradient( 0.25*zeta*(eta+zeta-1.0)*(15.0*eta*xi*zeta-5.0*eta*xi+2.0*eta*zeta*zeta-4.0*eta*zeta+2.0*eta-15.0*xi*zeta*zeta+20.0*xi*zeta-5.0*xi+6.0*zeta*zeta*zeta-14.0*zeta*zeta+10.0*zeta-2.0)/(den*den*den), 0.25*zeta*(xi-zeta+1.0)*(15.0*eta*xi*zeta-5.0*eta*xi+15.0*eta*zeta*zeta-20.0*eta*zeta+5.0*eta-2.0*xi*zeta*zeta+4.0*xi*zeta-2.0*xi+6.0*zeta*zeta*zeta-14.0*zeta*zeta+10.0*zeta-2.0)/(den*den*den), -0.125*(30.0*eta*eta-32.0*eta*xi+12.0*eta*zeta-4.0*eta+30.0*xi*xi-12.0*xi*zeta+4.0*xi+5.0*zeta-3.0)+0.125*(60.0*eta*eta*xi*xi*zeta*zeta-37.0*zeta*eta*eta*xi*xi+7.0*eta*eta*xi*xi+4.0*xi*eta*eta*zeta*zeta*zeta-8.0*xi*eta*eta*zeta*zeta+4.0*xi*zeta*eta*eta-33.0*eta*eta*zeta*zeta*zeta+89.0*eta*eta*zeta*zeta-79.0*zeta*eta*eta+23.0*eta*eta-4.0*eta*xi*xi*zeta*zeta*zeta+8.0*eta*xi*xi*zeta*zeta-4.0*eta*zeta*xi*xi+36.0*eta*xi*zeta*zeta*zeta-100.0*eta*xi*zeta*zeta+92.0*eta*xi*zeta-28.0*eta*xi-33.0*xi*xi*zeta*zeta*zeta+89.0*xi*xi*zeta*zeta-79.0*zeta*xi*xi+23.0*xi*xi)/(den*den*den*den) );
                case 12:
                  return sign * RealGradient( -0.125*zeta*(eta-zeta+1.0)*(3.0*eta*xi*zeta+7.0*eta*xi+2.0*eta*zeta*zeta-4.0*eta*zeta+2.0*eta+11.0*xi*zeta*zeta-12.0*xi*zeta+xi-6.0*zeta*zeta*zeta+14.0*zeta*zeta-10.0*zeta+2.0)/(den*den*den), -0.125*zeta*(xi-zeta+1.0)*(3.0*eta*xi*zeta+7.0*eta*xi+11.0*eta*zeta*zeta-12.0*eta*zeta+eta+2.0*xi*zeta*zeta-4.0*xi*zeta+2.0*xi-6.0*zeta*zeta*zeta+14.0*zeta*zeta-10.0*zeta+2.0)/(den*den*den), 0.125*(11.0*eta*eta+16.0*eta*xi-6.0*eta*zeta+11.0*xi*xi-6.0*xi*zeta+2.0*zeta-1.0)-0.125*(6.0*eta*eta*xi*xi*zeta*zeta+16.0*zeta*eta*eta*xi*xi-7.0*eta*eta*xi*xi+10.0*xi*eta*eta*zeta*zeta*zeta-28.0*xi*eta*eta*zeta*zeta+26.0*xi*zeta*eta*eta-8.0*xi*eta*eta-15.0*eta*eta*zeta*zeta*zeta+40.0*eta*eta*zeta*zeta-35.0*zeta*eta*eta+10.0*eta*eta+10.0*eta*xi*xi*zeta*zeta*zeta-28.0*eta*xi*xi*zeta*zeta+26.0*eta*zeta*xi*xi-8.0*eta*xi*xi-12.0*eta*xi*zeta*zeta*zeta+32.0*eta*xi*zeta*zeta-28.0*eta*xi*zeta+8.0*eta*xi-15.0*xi*xi*zeta*zeta*zeta+40.0*xi*xi*zeta*zeta-35.0*zeta*xi*xi+10.0*xi*xi)/(den*den*den*den) );
                case 13:
                  return sign * RealGradient( 0.25*zeta*(eta-zeta+1.0)*(15.0*eta*xi*zeta-5.0*eta*xi+2.0*eta*zeta*zeta-4.0*eta*zeta+2.0*eta+15.0*xi*zeta*zeta-20.0*xi*zeta+5.0*xi-6.0*zeta*zeta*zeta+14.0*zeta*zeta-10.0*zeta+2.0)/(den*den*den), 0.25*zeta*(xi-zeta+1.0)*(15.0*eta*xi*zeta-5.0*eta*xi+15.0*eta*zeta*zeta-20.0*eta*zeta+5.0*eta+2.0*xi*zeta*zeta-4.0*xi*zeta+2.0*xi-6.0*zeta*zeta*zeta+14.0*zeta*zeta-10.0*zeta+2.0)/(den*den*den), -0.125*(30.0*eta*eta+32.0*eta*xi-12.0*eta*zeta+4.0*eta+30.0*xi*xi-12.0*xi*zeta+4.0*xi+5.0*zeta-3.0)+0.125*(60.0*eta*eta*xi*xi*zeta*zeta-37.0*zeta*eta*eta*xi*xi+7.0*eta*eta*xi*xi+4.0*xi*eta*eta*zeta*zeta*zeta-8.0*xi*eta*eta*zeta*zeta+4.0*xi*zeta*eta*eta-33.0*eta*eta*zeta*zeta*zeta+89.0*eta*eta*zeta*zeta-79.0*zeta*eta*eta+23.0*eta*eta+4.0*eta*xi*xi*zeta*zeta*zeta-8.0*eta*xi*xi*zeta*zeta+4.0*eta*zeta*xi*xi-36.0*eta*xi*zeta*zeta*zeta+100.0*eta*xi*zeta*zeta-92.0*eta*xi*zeta+28.0*eta*xi-33.0*xi*xi*zeta*zeta*zeta+89.0*xi*xi*zeta*zeta-79.0*zeta*xi*xi+23.0*xi*xi)/(den*den*den*den) );
                case 14:
                  return sign * RealGradient( -0.125*zeta*(eta-zeta+1.0)*(3.0*eta*xi*zeta+7.0*eta*xi-2.0*eta*zeta*zeta+4.0*eta*zeta-2.0*eta+11.0*xi*zeta*zeta-12.0*xi*zeta+xi+6.0*zeta*zeta*zeta-14.0*zeta*zeta+10.0*zeta-2.0)/(den*den*den), -0.125*zeta*(xi+zeta-1.0)*(3.0*eta*xi*zeta+7.0*eta*xi-11.0*eta*zeta*zeta+12.0*eta*zeta-eta+2.0*xi*zeta*zeta-4.0*xi*zeta+2.0*xi+6.0*zeta*zeta*zeta-14.0*zeta*zeta+10.0*zeta-2.0)/(den*den*den), 0.125*(11.0*eta*eta-16.0*eta*xi-6.0*eta*zeta+11.0*xi*xi+6.0*xi*zeta+2.0*zeta-1.0)-0.125*(6.0*eta*eta*xi*xi*zeta*zeta+16.0*zeta*eta*eta*xi*xi-7.0*eta*eta*xi*xi-10.0*xi*eta*eta*zeta*zeta*zeta+28.0*xi*eta*eta*zeta*zeta-26.0*xi*zeta*eta*eta+8.0*xi*eta*eta-15.0*eta*eta*zeta*zeta*zeta+40.0*eta*eta*zeta*zeta-35.0*zeta*eta*eta+10.0*eta*eta+10.0*eta*xi*xi*zeta*zeta*zeta-28.0*eta*xi*xi*zeta*zeta+26.0*eta*zeta*xi*xi-8.0*eta*xi*xi+12.0*eta*xi*zeta*zeta*zeta-32.0*eta*xi*zeta*zeta+28.0*eta*xi*zeta-8.0*eta*xi-15.0*xi*xi*zeta*zeta*zeta+40.0*xi*xi*zeta*zeta-35.0*zeta*xi*xi+10.0*xi*xi)/(den*den*den*den) );
                case 15:
                  return sign * RealGradient( 0.25*zeta*(eta-zeta+1.0)*(15.0*eta*xi*zeta-5.0*eta*xi-2.0*eta*zeta*zeta+4.0*eta*zeta-2.0*eta+15.0*xi*zeta*zeta-20.0*xi*zeta+5.0*xi+6.0*zeta*zeta*zeta-14.0*zeta*zeta+10.0*zeta-2.0)/(den*den*den), 0.25*zeta*(xi+zeta-1.0)*(15.0*eta*xi*zeta-5.0*eta*xi-15.0*eta*zeta*zeta+20.0*eta*zeta-5.0*eta+2.0*xi*zeta*zeta-4.0*xi*zeta+2.0*xi+6.0*zeta*zeta*zeta-14.0*zeta*zeta+10.0*zeta-2.0)/(den*den*den), -0.125*(30.0*eta*eta-32.0*eta*xi-12.0*eta*zeta+4.0*eta+30.0*xi*xi+12.0*xi*zeta-4.0*xi+5.0*zeta-3.0)+0.125*(60.0*eta*eta*xi*xi*zeta*zeta-37.0*zeta*eta*eta*xi*xi+7.0*eta*eta*xi*xi-4.0*xi*eta*eta*zeta*zeta*zeta+8.0*xi*eta*eta*zeta*zeta-4.0*xi*zeta*eta*eta-33.0*eta*eta*zeta*zeta*zeta+89.0*eta*eta*zeta*zeta-79.0*zeta*eta*eta+23.0*eta*eta+4.0*eta*xi*xi*zeta*zeta*zeta-8.0*eta*xi*xi*zeta*zeta+4.0*eta*zeta*xi*xi+36.0*eta*xi*zeta*zeta*zeta-100.0*eta*xi*zeta*zeta+92.0*eta*xi*zeta-28.0*eta*xi-33.0*xi*xi*zeta*zeta*zeta+89.0*xi*xi*zeta*zeta-79.0*zeta*xi*xi+23.0*xi*xi)/(den*den*den*den) );
                case 16:
                case 17:
                case 18:
                case 19:
                  {
                    const std::array<Real, 4> c =
                      nedelec_one_face_coefficients(*elem, 4, i - 16);
                    return c[0] * RealGradient( 1.5*xi*zeta*(3.0*zeta-1.0)*(eta-zeta+1.0)*(eta+zeta-1.0)/(den*den*den), 0.75*(xi-zeta+1.0)*(xi+zeta-1.0)*(24.0*eta*zeta*zeta-23.0*eta*zeta+3.0*eta+5.0*zeta*zeta*zeta-11.0*zeta*zeta+7.0*zeta-1.0)/(den*den*den), -0.375*(48.0*eta*eta+10.0*eta*zeta-2.0*eta+12.0*xi*xi+5.0*zeta-1.0)+0.375*(60.0*eta*eta*xi*xi*zeta*zeta-55.0*zeta*eta*eta*xi*xi+7.0*eta*eta*xi*xi-45.0*eta*eta*zeta*zeta*zeta+131.0*eta*eta*zeta*zeta-127.0*zeta*eta*eta+41.0*eta*eta+10.0*eta*xi*xi*zeta*zeta*zeta-22.0*eta*xi*xi*zeta*zeta+14.0*eta*zeta*xi*xi-2.0*eta*xi*xi-15.0*xi*xi*zeta*zeta*zeta+41.0*xi*xi*zeta*zeta-37.0*zeta*xi*xi+11.0*xi*xi)/(den*den*den*den) ) +
                           c[1] * RealGradient( -1.5*xi*zeta*(3.0*zeta-1.0)*(eta-zeta+1.0)*(eta+zeta-1.0)/(den*den*den), -0.75*(xi-zeta+1.0)*(xi+zeta-1.0)*(24.0*eta*zeta*zeta-23.0*eta*zeta+3.0*eta-5.0*zeta*zeta*zeta+11.0*zeta*zeta-7.0*zeta+1.0)/(den*den*den), 0.375*(48.0*eta*eta-10.0*eta*zeta+2.0*eta+12.0*xi*xi+5.0*zeta-1.0)-0.375*(60.0*eta*eta*xi*xi*zeta*zeta-55.0*zeta*eta*eta*xi*xi+7.0*eta*eta*xi*xi-45.0*eta*eta*zeta*zeta*zeta+131.0*eta*eta*zeta*zeta-127.0*zeta*eta*eta+41.0*eta*eta-10.0*eta*xi*xi*zeta*zeta*zeta+22.0*eta*xi*xi*zeta*zeta-14.0*eta*zeta*xi*xi+2.0*eta*xi*xi-15.0*xi*xi*zeta*zeta*zeta+41.0*xi*xi*zeta*zeta-37.0*zeta*xi*xi+11.0*xi*xi)/(den*den*den*den) ) +
                           c[2] * RealGradient( 0.75*(eta-zeta+1.0)*(eta+zeta-1.0)*(24.0*xi*zeta*zeta-23.0*xi*zeta+3.0*xi+5.0*zeta*zeta*zeta-11.0*zeta*zeta+7.0*zeta-1.0)/(den*den*den), 1.5*eta*zeta*(3.0*zeta-1.0)*(xi-zeta+1.0)*(xi+zeta-1.0)/(den*den*den), -0.375*(12.0*eta*eta+48.0*xi*xi+10.0*xi*zeta-2.0*xi+5.0*zeta-1.0)+0.375*(60.0*eta*eta*xi*xi*zeta*zeta-55.0*zeta*eta*eta*xi*xi+7.0*eta*eta*xi*xi+10.0*xi*eta*eta*zeta*zeta*zeta-22.0*xi*eta*eta*zeta*zeta+14.0*xi*zeta*eta*eta-2.0*xi*eta*eta-15.0*eta*eta*zeta*zeta*zeta+41.0*eta*eta*zeta*zeta-37.0*zeta*eta*eta+11.0*eta*eta-45.0*xi*xi*zeta*zeta*zeta+131.0*xi*xi*zeta*zeta-127.0*zeta*xi*xi+41.0*xi*xi)/(den*den*den*den) ) +
                           c[3] * RealGradient( -0.75*(eta-zeta+1.0)*(eta+zeta-1.0)*(24.0*xi*zeta*zeta-23.0*xi*zeta+3.0*xi-5.0*zeta*zeta*zeta+11.0*zeta*zeta-7.0*zeta+1.0)/(den*den*den), -1.5*eta*zeta*(3.0*zeta-1.0)*(xi-zeta+1.0)*(xi+zeta-1.0)/(den*den*den), 0.375*(12.0*eta*eta+48.0*xi*xi-10.0*xi*zeta+2.0*xi+5.0*zeta-1.0)-0.375*(60.0*eta*eta*xi*xi*zeta*zeta-55.0*zeta*eta*eta*xi*xi+7.0*eta*eta*xi*xi-10.0*xi*eta*eta*zeta*zeta*zeta+22.0*xi*eta*eta*zeta*zeta-14.0*xi*zeta*eta*eta+2.0*xi*eta*eta-15.0*eta*eta*zeta*zeta*zeta+41.0*eta*eta*zeta*zeta-37.0*zeta*eta*eta+11.0*eta*eta-45.0*xi*xi*zeta*zeta*zeta+131.0*xi*xi*zeta*zeta-127.0*zeta*xi*xi+41.0*xi*xi)/(den*den*den*den) );
                  }
                case 20:
                case 21:
                  {
                    const std::array<Real, 4> c =
                      nedelec_one_face_coefficients(*elem, 0, i - 20);
                    return c[0] * RealGradient( -0.5*zeta*(3.0*eta-zeta+1.0)*(eta+zeta-1.0)*(xi+3.0*zeta-3.0)/(den*den), 0.0, 0.5*(3.0*eta-zeta+1.0)*(eta+zeta-1.0)*(xi*xi+3.0*xi*zeta+zeta-1.0)/(den*den) ) +
                           c[1] * RealGradient( xi*zeta*(3.0*eta-zeta+1.0)*(eta+zeta-1.0)/(den*den), 0.0, -(3.0*eta-zeta+1.0)*(eta+zeta-1.0)*(xi*xi+zeta-1.0)/(den*den) );
                  }
                case 22:
                case 23:
                  {
                    const std::array<Real, 4> c =
                      nedelec_one_face_coefficients(*elem, 1, i - 22);
                    return c[0] * RealGradient( 0.0, -0.5*zeta*(xi-zeta+1.0)*(3.0*xi+zeta-1.0)*(eta+3.0*zeta-3.0)/(den*den), 0.5*(xi-zeta+1.0)*(3.0*xi+zeta-1.0)*(eta*eta+3.0*eta*zeta+zeta-1.0)/(den*den) ) +
                           c[1] * RealGradient( 0.0, eta*zeta*(xi-zeta+1.0)*(3.0*xi+zeta-1.0)/(den*den), -(xi-zeta+1.0)*(eta*eta+zeta-1.0)*(3.0*xi+zeta-1.0)/(den*den) );
                  }
                case 24:
                case 25:
                  {
                    const std::array<Real, 4> c =
                      nedelec_one_face_coefficients(*elem, 2, i - 24);
                    return c[0] * RealGradient( -0.5*zeta*(eta-zeta+1.0)*(3.0*eta+zeta-1.0)*(xi-3.0*zeta+3.0)/(den*den), 0.0, 0.5*(eta-zeta+1.0)*(3.0*eta+zeta-1.0)*(xi*xi-3.0*xi*zeta+zeta-1.0)/(den*den) ) +
                           c[1] * RealGradient( xi*zeta*(eta-zeta+1.0)*(3.0*eta+zeta-1.0)/(den*den), 0.0, -(eta-zeta+1.0)*(xi*xi+zeta-1.0)*(3.0*eta+zeta-1.0)/(den*den) );
                  }
                case 26:
                case 27:
                  {
                    const std::array<Real, 4> c =
                      nedelec_one_face_coefficients(*elem, 3, i - 26);
                    return c[0] * RealGradient( 0.0, -0.5*zeta*(3.0*xi-zeta+1.0)*(xi+zeta-1.0)*(eta-3.0*zeta+3.0)/(den*den), 0.5*(3.0*xi-zeta+1.0)*(xi+zeta-1.0)*(eta*eta-3.0*eta*zeta+zeta-1.0)/(den*den) ) +
                           c[1] * RealGradient( 0.0, eta*zeta*(3.0*xi-zeta+1.0)*(xi+zeta-1.0)/(den*den), -(3.0*xi-zeta+1.0)*(xi+zeta-1.0)*(eta*eta+zeta-1.0)/(den*den) );
                  }
                case 28:
                  return RealGradient( -7.5*zeta*(eta-zeta+1.0)*(eta+zeta-1.0)/den, 0.0, 7.5*xi*zeta-7.5*xi*zeta*eta*eta/(den*den) );
                case 29:
                  return RealGradient( 0.0, 7.5*(xi*xi-zeta*zeta+zeta)-7.5*xi*xi/den, 7.5*eta*zeta-7.5*eta*zeta*xi*xi/(den*den) );
                case 30:
                  return RealGradient( 3.75*xi*zeta*(3.0*zeta-1.0)*(eta-zeta+1.0)*(eta+zeta-1.0)/(den*den*den), 3.75*eta*zeta*(3.0*zeta-1.0)*(xi-zeta+1.0)*(xi+zeta-1.0)/(den*den*den), -3.75*(3.0*eta*eta+3.0*xi*xi+2.0*zeta-1.0)+3.75*(6.0*eta*eta*xi*xi*zeta*zeta-4.0*zeta*eta*eta*xi*xi+eta*eta*xi*xi-3.0*eta*eta*zeta*zeta*zeta+8.0*eta*eta*zeta*zeta-7.0*zeta*eta*eta+2.0*eta*eta-3.0*xi*xi*zeta*zeta*zeta+8.0*xi*xi*zeta*zeta-7.0*zeta*xi*xi+2.0*xi*xi)/(den*den*den*den) );
                case 31:
                  return RealGradient( 7.5*xi*zeta*(6.0*zeta-5.0)*(eta-zeta+1.0)*(eta+zeta-1.0)/(den*den*den), 3.75*eta*zeta*(3.0*zeta-1.0)*(xi-zeta+1.0)*(xi+zeta-1.0)/(den*den*den), -1.875*(6.0*eta*eta+24.0*xi*xi+zeta+1.0)+1.875*(30.0*eta*eta*xi*xi*zeta*zeta-23.0*zeta*eta*eta*xi*xi-eta*eta*xi*xi-9.0*eta*eta*zeta*zeta*zeta+25.0*eta*eta*zeta*zeta-23.0*zeta*eta*eta+7.0*eta*eta-27.0*xi*xi*zeta*zeta*zeta+79.0*xi*xi*zeta*zeta-77.0*zeta*xi*xi+25.0*xi*xi)/(den*den*den*den) );
                case 32:
                  return RealGradient( 3.75*xi*zeta*(3.0*zeta-1.0)*(eta-zeta+1.0)*(eta+zeta-1.0)/(den*den*den), 7.5*eta*zeta*(6.0*zeta-5.0)*(xi-zeta+1.0)*(xi+zeta-1.0)/(den*den*den), -1.875*(24.0*eta*eta+6.0*xi*xi+zeta+1.0)+1.875*(30.0*eta*eta*xi*xi*zeta*zeta-23.0*zeta*eta*eta*xi*xi-eta*eta*xi*xi-27.0*eta*eta*zeta*zeta*zeta+79.0*eta*eta*zeta*zeta-77.0*zeta*eta*eta+25.0*eta*eta-9.0*xi*xi*zeta*zeta*zeta+25.0*xi*xi*zeta*zeta-23.0*zeta*xi*xi+7.0*xi*xi)/(den*den*den*den) );
                case 33:
                  return RealGradient( -18.75*xi*zeta*(3.0*zeta-1.0)*(eta-zeta+1.0)*(eta+zeta-1.0)/(den*den*den), -18.75*eta*zeta*(3.0*zeta-1.0)*(xi-zeta+1.0)*(xi+zeta-1.0)/(den*den*den), 3.75*(15.0*eta*eta+15.0*xi*xi+7.0*zeta-2.0)-3.75*(30.0*eta*eta*xi*xi*zeta*zeta-17.0*zeta*eta*eta*xi*xi+2.0*eta*eta*xi*xi-18.0*eta*eta*zeta*zeta*zeta+49.0*eta*eta*zeta*zeta-44.0*zeta*eta*eta+13.0*eta*eta-18.0*xi*xi*zeta*zeta*zeta+49.0*xi*xi*zeta*zeta-44.0*zeta*xi*xi+13.0*xi*xi)/(den*den*den*den) );

                default:
                  libmesh_error_msg("Invalid i = " << i);
                } // switch(ii)
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 3D element type!: " << Utility::enum_to_string(elem->type()));
          }
      }

