      // quadratic Nedelec (first kind) shape function first derivatives
    case SECOND:
      {
        switch (elem->type())
          {
          case TET14:
            {
              switch (j)
                {
                  // d()/dxi
                case 0:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient( 2.0*(4.0*eta+4.0*zeta-3.0), -2.0*(4.0*eta+8.0*xi+4.0*zeta-3.0), -2.0*(4.0*eta+8.0*xi+4.0*zeta-3.0) );
                      case 1:
                        return sign * RealGradient( -2.0*(4.0*eta+4.0*zeta-3.0), 4.0*(4.0*xi-1.0), 4.0*(4.0*xi-1.0) );
                      case 2:
                        return sign * RealGradient( -8.0*eta, 4.0*(4.0*xi-1.0), 0.0 );
                      case 3:
                        return sign * RealGradient( 0.0, 2.0*(4.0*eta-1.0), 0.0 );
                      case 4:
                        return sign * RealGradient( -8.0*eta, 4.0*(2.0*eta+4.0*xi+4.0*zeta-3.0), -8.0*eta );
                      case 5:
                        return sign * RealGradient( 0.0, -2.0*(4.0*eta-1.0), 0.0 );
                      case 6:
                        return sign * RealGradient( -8.0*zeta, -8.0*zeta, 4.0*(4.0*eta+4.0*xi+2.0*zeta-3.0) );
                      case 7:
                        return sign * RealGradient( 0.0, 0.0, -2.0*(4.0*zeta-1.0) );
                      case 8:
                        return sign * RealGradient( -8.0*zeta, 0.0, 4.0*(4.0*xi-1.0) );
                      case 9:
                        return sign * RealGradient( 0.0, 0.0, 2.0*(4.0*zeta-1.0) );
                      case 12:
                      case 13:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 0, i - 12);
                          return c[0] * RealGradient( 16.0*eta, -8.0*(eta+4.0*xi+2.0*zeta-2.0), 8.0*eta ) +
                                 c[1] * RealGradient( -8.0*eta, 8.0*(2.0*eta+2.0*xi+zeta-1.0), 8.0*eta );
                        }
                      case 14:
                      case 15:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 1, i - 14);
                          return c[0] * RealGradient( -8.0*zeta, 8.0*zeta, 8.0*(eta+2.0*xi+2.0*zeta-1.0) ) +
                                 c[1] * RealGradient( 16.0*zeta, 8.0*zeta, -8.0*(2.0*eta+4.0*xi+zeta-2.0) );
                        }
                      case 16:
                      case 17:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 2, i - 16);
                          return c[0] * RealGradient( 0.0, 16.0*zeta, -8.0*eta ) +
                                 c[1] * RealGradient( 0.0, -8.0*zeta, 16.0*eta );
                        }
                      case 18:
                      case 19:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 3, i - 18);
                          return c[0] * RealGradient( 0.0, 8.0*zeta, 8.0*eta ) +
                                 c[1] * RealGradient( 0.0, 8.0*zeta, -16.0*eta );
                        }
                      case 10:
                      case 11:
                        return RealGradient();

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      } // switch(ii)
                  } // j = 0

                  // d()/deta
                case 1:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient( 4.0*(4.0*eta+2.0*xi+4.0*zeta-3.0), -8.0*xi, -8.0*xi );
                      case 1:
                        return sign * RealGradient( -2.0*(4.0*xi-1.0), 0.0, 0.0 );
                      case 2:
                        return sign * RealGradient( -2.0*(4.0*xi-1.0), 0.0, 0.0 );
                      case 3:
                        return sign * RealGradient( -4.0*(4.0*eta-1.0), 8.0*xi, 0.0 );
                      case 4:
                        return sign * RealGradient( -2.0*(8.0*eta+4.0*xi+4.0*zeta-3.0), 2.0*(4.0*xi+4.0*zeta-3.0), -2.0*(8.0*eta+4.0*xi+4.0*zeta-3.0) );
                      case 5:
                        return sign * RealGradient( 4.0*(4.0*eta-1.0), -2.0*(4.0*xi+4.0*zeta-3.0), 4.0*(4.0*eta-1.0) );
                      case 6:
                        return sign * RealGradient( -8.0*zeta, -8.0*zeta, 4.0*(4.0*eta+4.0*xi+2.0*zeta-3.0) );
                      case 7:
                        return sign * RealGradient( 0.0, 0.0, -2.0*(4.0*zeta-1.0) );
                      case 10:
                        return sign * RealGradient( 0.0, -8.0*zeta, 4.0*(4.0*eta-1.0) );
                      case 11:
                        return sign * RealGradient( 0.0, 0.0, 2.0*(4.0*zeta-1.0) );
                      case 12:
                      case 13:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 0, i - 12);
                          return c[0] * RealGradient( 8.0*(2.0*eta+2.0*xi+zeta-1.0), -8.0*xi, 8.0*xi ) +
                                 c[1] * RealGradient( -8.0*(4.0*eta+xi+2.0*zeta-2.0), 16.0*xi, 8.0*xi );
                        }
                      case 14:
                      case 15:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 1, i - 14);
                          return c[0] * RealGradient( -16.0*zeta, 0.0, 8.0*xi ) +
                                 c[1] * RealGradient( 8.0*zeta, 0.0, -16.0*xi );
                        }
                      case 16:
                      case 17:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 2, i - 16);
                          return c[0] * RealGradient( -8.0*zeta, 0.0, -8.0*xi ) +
                                 c[1] * RealGradient( -8.0*zeta, 0.0, 16.0*xi );
                        }
                      case 18:
                      case 19:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 3, i - 18);
                          return c[0] * RealGradient( -16.0*zeta, -8.0*zeta, 8.0*(2.0*eta+xi-zeta-1.0) ) +
                                 c[1] * RealGradient( 8.0*zeta, 16.0*zeta, -8.0*(4.0*eta+2.0*xi+zeta-2.0) );
                        }
                      case 8:
                      case 9:
                        return RealGradient();

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      } // switch(ii)
                  } // j = 1

                  // d()/dzeta
                case 2:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient( 4.0*(4.0*eta+2.0*xi+4.0*zeta-3.0), -8.0*xi, -8.0*xi );
                      case 1:
                        return sign * RealGradient( -2.0*(4.0*xi-1.0), 0.0, 0.0 );
                      case 4:
                        return sign * RealGradient( -8.0*eta, 4.0*(2.0*eta+4.0*xi+4.0*zeta-3.0), -8.0*eta );
                      case 5:
                        return sign * RealGradient( 0.0, -2.0*(4.0*eta-1.0), 0.0 );
                      case 6:
                        return sign * RealGradient( -2.0*(4.0*eta+4.0*xi+8.0*zeta-3.0), -2.0*(4.0*eta+4.0*xi+8.0*zeta-3.0), 2.0*(4.0*eta+4.0*xi-3.0) );
                      case 7:
                        return sign * RealGradient( 4.0*(4.0*zeta-1.0), 4.0*(4.0*zeta-1.0), -2.0*(4.0*eta+4.0*xi-3.0) );
                      case 8:
                        return sign * RealGradient( -2.0*(4.0*xi-1.0), 0.0, 0.0 );
                      case 9:
                        return sign * RealGradient( -4.0*(4.0*zeta-1.0), 0.0, 8.0*xi );
                      case 10:
                        return sign * RealGradient( 0.0, -2.0*(4.0*eta-1.0), 0.0 );
                      case 11:
                        return sign * RealGradient( 0.0, -4.0*(4.0*zeta-1.0), 8.0*eta );
                      case 12:
                      case 13:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 0, i - 12);
                          return c[0] * RealGradient( 8.0*eta, -16.0*xi, 0.0 ) +
                                 c[1] * RealGradient( -16.0*eta, 8.0*xi, 0.0 );
                        }
                      case 14:
                      case 15:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 1, i - 14);
                          return c[0] * RealGradient( -8.0*(2.0*eta+xi+4.0*zeta-2.0), 8.0*xi, 16.0*xi ) +
                                 c[1] * RealGradient( 8.0*(eta+2.0*xi+2.0*zeta-1.0), 8.0*xi, -8.0*xi );
                        }
                      case 16:
                      case 17:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 2, i - 16);
                          return c[0] * RealGradient( -8.0*eta, 16.0*xi, 0.0 ) +
                                 c[1] * RealGradient( -8.0*eta, -8.0*xi, 0.0 );
                        }
                      case 18:
                      case 19:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 3, i - 18);
                          return c[0] * RealGradient( -16.0*eta, -8.0*(eta-xi-2.0*zeta+1.0), -8.0*eta ) +
                                 c[1] * RealGradient( 8.0*eta, 8.0*(2.0*eta+xi+2.0*zeta-1.0), -8.0*eta );
                        }
                      case 2:
                      case 3:
                        return RealGradient();

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      } // switch(ii)
                  } // j = 2

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          case HEX27:
            {
              // Even with a loose inverse_map tolerance we ought to
              // be nearly on the element interior in master
              // coordinates
              libmesh_assert_less_equal ( std::fabs(xi),   1.0+10*TOLERANCE );
              libmesh_assert_less_equal ( std::fabs(eta),  1.0+10*TOLERANCE );
              libmesh_assert_less_equal ( std::fabs(zeta), 1.0+10*TOLERANCE );
              switch (j)
                {
                  // d()/dxi
                case 0:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient( -0.09375*(3.0*eta+1.0)*(3.0*zeta+1.0)*(eta-1.0)*(zeta-1.0), 0.0, 0.0 );
                      case 1:
                        return sign * RealGradient( 0.09375*(3.0*eta+1.0)*(3.0*zeta+1.0)*(eta-1.0)*(zeta-1.0), 0.0, 0.0 );
                      case 2:
                        return sign * RealGradient( 0.0, -0.0625*(3.0*xi+1.0)*(3.0*zeta+1.0)*(zeta-1.0)*(3.0*eta-1.0), 0.0 );
                      case 3:
                        return sign * RealGradient( 0.0, 0.0625*(3.0*eta+1.0)*(3.0*xi+1.0)*(3.0*zeta+1.0)*(zeta-1.0), 0.0 );
                      case 4:
                        return sign * RealGradient( -0.09375*(eta+1.0)*(3.0*zeta+1.0)*(zeta-1.0)*(3.0*eta-1.0), 0.0, 0.0 );
                      case 5:
                        return sign * RealGradient( 0.09375*(eta+1.0)*(3.0*zeta+1.0)*(zeta-1.0)*(3.0*eta-1.0), 0.0, 0.0 );
                      case 6:
                        return sign * RealGradient( 0.0, -0.0625*(3.0*zeta+1.0)*(zeta-1.0)*(3.0*eta-1.0)*(3.0*xi-1.0), 0.0 );
                      case 7:
                        return sign * RealGradient( 0.0, 0.0625*(3.0*eta+1.0)*(3.0*zeta+1.0)*(zeta-1.0)*(3.0*xi-1.0), 0.0 );
                      case 8:
                        return sign * RealGradient( 0.0, 0.0, -0.0625*(3.0*eta+1.0)*(eta-1.0)*(3.0*xi-1.0)*(3.0*zeta-1.0) );
                      case 9:
                        return sign * RealGradient( 0.0, 0.0, 0.0625*(3.0*eta+1.0)*(3.0*zeta+1.0)*(eta-1.0)*(3.0*xi-1.0) );
                      case 10:
                        return sign * RealGradient( 0.0, 0.0, -0.0625*(3.0*eta+1.0)*(3.0*xi+1.0)*(eta-1.0)*(3.0*zeta-1.0) );
                      case 11:
                        return sign * RealGradient( 0.0, 0.0, 0.0625*(3.0*eta+1.0)*(3.0*xi+1.0)*(3.0*zeta+1.0)*(eta-1.0) );
                      case 12:
                        return sign * RealGradient( 0.0, 0.0, -0.0625*(eta+1.0)*(3.0*xi+1.0)*(3.0*eta-1.0)*(3.0*zeta-1.0) );
                      case 13:
                        return sign * RealGradient( 0.0, 0.0, 0.0625*(eta+1.0)*(3.0*xi+1.0)*(3.0*zeta+1.0)*(3.0*eta-1.0) );
                      case 14:
                        return sign * RealGradient( 0.0, 0.0, -0.0625*(eta+1.0)*(3.0*eta-1.0)*(3.0*xi-1.0)*(3.0*zeta-1.0) );
                      case 15:
                        return sign * RealGradient( 0.0, 0.0, 0.0625*(eta+1.0)*(3.0*zeta+1.0)*(3.0*eta-1.0)*(3.0*xi-1.0) );
                      case 16:
                        return sign * RealGradient( -0.09375*(zeta+1.0)*(3.0*eta+1.0)*(eta-1.0)*(3.0*zeta-1.0), 0.0, 0.0 );
                      case 17:
                        return sign * RealGradient( 0.09375*(zeta+1.0)*(3.0*eta+1.0)*(eta-1.0)*(3.0*zeta-1.0), 0.0, 0.0 );
                      case 18:
                        return sign * RealGradient( 0.0, -0.0625*(zeta+1.0)*(3.0*xi+1.0)*(3.0*eta-1.0)*(3.0*zeta-1.0), 0.0 );
                      case 19:
                        return sign * RealGradient( 0.0, 0.0625*(zeta+1.0)*(3.0*eta+1.0)*(3.0*xi+1.0)*(3.0*zeta-1.0), 0.0 );
                      case 20:
                        return sign * RealGradient( -0.09375*(eta+1.0)*(zeta+1.0)*(3.0*eta-1.0)*(3.0*zeta-1.0), 0.0, 0.0 );
                      case 21:
                        return sign * RealGradient( 0.09375*(eta+1.0)*(zeta+1.0)*(3.0*eta-1.0)*(3.0*zeta-1.0), 0.0, 0.0 );
                      case 22:
                        return sign * RealGradient( 0.0, -0.0625*(zeta+1.0)*(3.0*eta-1.0)*(3.0*xi-1.0)*(3.0*zeta-1.0), 0.0 );
                      case 23:
                        return sign * RealGradient( 0.0, 0.0625*(zeta+1.0)*(3.0*eta+1.0)*(3.0*xi-1.0)*(3.0*zeta-1.0), 0.0 );
                      case 24:
                      case 25:
                      case 26:
                      case 27:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 0, i - 24);
                          return c[0] * RealGradient( 0.0, 0.375*xi*(3.0*zeta+1.0)*(zeta-1.0)*(3.0*eta-1.0), 0.0 ) +
                                 c[1] * RealGradient( 0.0, -0.375*xi*(3.0*eta+1.0)*(3.0*zeta+1.0)*(zeta-1.0), 0.0 ) +
                                 c[2] * RealGradient( 0.5625*(eta+1.0)*(3.0*zeta+1.0)*(eta-1.0)*(zeta-1.0), 0.0, 0.0 ) +
                                 c[3] * RealGradient( -0.5625*(eta+1.0)*(3.0*zeta+1.0)*(eta-1.0)*(zeta-1.0), 0.0, 0.0 );
                        }
                      case 28:
                      case 29:
                      case 30:
                      case 31:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 1, i - 28);
                          return c[0] * RealGradient( 0.5625*(zeta+1.0)*(3.0*eta+1.0)*(eta-1.0)*(zeta-1.0), 0.0, 0.0 ) +
                                 c[1] * RealGradient( -0.5625*(zeta+1.0)*(3.0*eta+1.0)*(eta-1.0)*(zeta-1.0), 0.0, 0.0 ) +
                                 c[2] * RealGradient( 0.0, 0.0, 0.375*xi*(3.0*eta+1.0)*(eta-1.0)*(3.0*zeta-1.0) ) +
                                 c[3] * RealGradient( 0.0, 0.0, -0.375*xi*(3.0*eta+1.0)*(3.0*zeta+1.0)*(eta-1.0) );
                        }
                      case 32:
                      case 33:
                      case 34:
                      case 35:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 2, i - 32);
                          return c[0] * RealGradient( 0.0, 0.375*(zeta+1.0)*(3.0*xi+1.0)*(zeta-1.0)*(3.0*eta-1.0), 0.0 ) +
                                 c[1] * RealGradient( 0.0, -0.375*(zeta+1.0)*(3.0*eta+1.0)*(3.0*xi+1.0)*(zeta-1.0), 0.0 ) +
                                 c[2] * RealGradient( 0.0, 0.0, 0.375*(eta+1.0)*(3.0*xi+1.0)*(eta-1.0)*(3.0*zeta-1.0) ) +
                                 c[3] * RealGradient( 0.0, 0.0, -0.375*(eta+1.0)*(3.0*xi+1.0)*(3.0*zeta+1.0)*(eta-1.0) );
                        }
                      case 36:
                      case 37:
                      case 38:
                      case 39:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 3, i - 36);
                          return c[0] * RealGradient( 0.5625*(eta+1.0)*(zeta+1.0)*(zeta-1.0)*(3.0*eta-1.0), 0.0, 0.0 ) +
                                 c[1] * RealGradient( -0.5625*(eta+1.0)*(zeta+1.0)*(zeta-1.0)*(3.0*eta-1.0), 0.0, 0.0 ) +
                                 c[2] * RealGradient( 0.0, 0.0, 0.375*xi*(eta+1.0)*(3.0*eta-1.0)*(3.0*zeta-1.0) ) +
                                 c[3] * RealGradient( 0.0, 0.0, -0.375*xi*(eta+1.0)*(3.0*zeta+1.0)*(3.0*eta-1.0) );
                        }
                      case 40:
                      case 41:
                      case 42:
                      case 43:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 4, i - 40);
                          return c[0] * RealGradient( 0.0, 0.375*(zeta+1.0)*(3.0*eta+1.0)*(zeta-1.0)*(3.0*xi-1.0), 0.0 ) +
                                 c[1] * RealGradient( 0.0, -0.375*(zeta+1.0)*(zeta-1.0)*(3.0*eta-1.0)*(3.0*xi-1.0), 0.0 ) +
                                 c[2] * RealGradient( 0.0, 0.0, 0.375*(eta+1.0)*(eta-1.0)*(3.0*xi-1.0)*(3.0*zeta-1.0) ) +
                                 c[3] * RealGradient( 0.0, 0.0, -0.375*(eta+1.0)*(3.0*zeta+1.0)*(eta-1.0)*(3.0*xi-1.0) );
                        }
                      case 44:
                      case 45:
                      case 46:
                      case 47:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 5, i - 44);
                          return c[0] * RealGradient( 0.5625*(eta+1.0)*(zeta+1.0)*(eta-1.0)*(3.0*zeta-1.0), 0.0, 0.0 ) +
                                 c[1] * RealGradient( -0.5625*(eta+1.0)*(zeta+1.0)*(eta-1.0)*(3.0*zeta-1.0), 0.0, 0.0 ) +
                                 c[2] * RealGradient( 0.0, 0.375*xi*(zeta+1.0)*(3.0*eta-1.0)*(3.0*zeta-1.0), 0.0 ) +
                                 c[3] * RealGradient( 0.0, -0.375*xi*(zeta+1.0)*(3.0*eta+1.0)*(3.0*zeta-1.0), 0.0 );
                        }
                      case 48:
                        return RealGradient( -0.84375*(eta+1.0)*(zeta+1.0)*(eta-1.0)*(zeta-1.0), 0.0, 0.0 );
                      case 49:
                        return RealGradient( 0.84375*(eta+1.0)*(zeta+1.0)*(eta-1.0)*(zeta-1.0), 0.0, 0.0 );
                      case 50:
                        return RealGradient( 0.0, -0.5625*xi*(zeta+1.0)*(zeta-1.0)*(3.0*eta-1.0), 0.0 );
                      case 51:
                        return RealGradient( 0.0, 0.5625*xi*(zeta+1.0)*(3.0*eta+1.0)*(zeta-1.0), 0.0 );
                      case 52:
                        return RealGradient( 0.0, 0.0, -0.5625*xi*(eta+1.0)*(eta-1.0)*(3.0*zeta-1.0) );
                      case 53:
                        return RealGradient( 0.0, 0.0, 0.5625*xi*(eta+1.0)*(3.0*zeta+1.0)*(eta-1.0) );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      } // switch(ii)
                  } // j = 0

                  // d()/deta
                case 1:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient( -0.0625*(3.0*zeta+1.0)*(zeta-1.0)*(3.0*eta-1.0)*(3.0*xi-1.0), 0.0, 0.0 );
                      case 1:
                        return sign * RealGradient( 0.0625*(3.0*xi+1.0)*(3.0*zeta+1.0)*(zeta-1.0)*(3.0*eta-1.0), 0.0, 0.0 );
                      case 2:
                        return sign * RealGradient( 0.0, -0.09375*(xi+1.0)*(3.0*zeta+1.0)*(zeta-1.0)*(3.0*xi-1.0), 0.0 );
                      case 3:
                        return sign * RealGradient( 0.0, 0.09375*(xi+1.0)*(3.0*zeta+1.0)*(zeta-1.0)*(3.0*xi-1.0), 0.0 );
                      case 4:
                        return sign * RealGradient( -0.0625*(3.0*eta+1.0)*(3.0*xi+1.0)*(3.0*zeta+1.0)*(zeta-1.0), 0.0, 0.0 );
                      case 5:
                        return sign * RealGradient( 0.0625*(3.0*eta+1.0)*(3.0*zeta+1.0)*(zeta-1.0)*(3.0*xi-1.0), 0.0, 0.0 );
                      case 6:
                        return sign * RealGradient( 0.0, -0.09375*(3.0*xi+1.0)*(3.0*zeta+1.0)*(xi-1.0)*(zeta-1.0), 0.0 );
                      case 7:
                        return sign * RealGradient( 0.0, 0.09375*(3.0*xi+1.0)*(3.0*zeta+1.0)*(xi-1.0)*(zeta-1.0), 0.0 );
                      case 8:
                        return sign * RealGradient( 0.0, 0.0, -0.0625*(3.0*xi+1.0)*(xi-1.0)*(3.0*eta-1.0)*(3.0*zeta-1.0) );
                      case 9:
                        return sign * RealGradient( 0.0, 0.0, 0.0625*(3.0*xi+1.0)*(3.0*zeta+1.0)*(xi-1.0)*(3.0*eta-1.0) );
                      case 10:
                        return sign * RealGradient( 0.0, 0.0, -0.0625*(xi+1.0)*(3.0*eta-1.0)*(3.0*xi-1.0)*(3.0*zeta-1.0) );
                      case 11:
                        return sign * RealGradient( 0.0, 0.0, 0.0625*(xi+1.0)*(3.0*zeta+1.0)*(3.0*eta-1.0)*(3.0*xi-1.0) );
                      case 12:
                        return sign * RealGradient( 0.0, 0.0, -0.0625*(xi+1.0)*(3.0*eta+1.0)*(3.0*xi-1.0)*(3.0*zeta-1.0) );
                      case 13:
                        return sign * RealGradient( 0.0, 0.0, 0.0625*(xi+1.0)*(3.0*eta+1.0)*(3.0*zeta+1.0)*(3.0*xi-1.0) );
                      case 14:
                        return sign * RealGradient( 0.0, 0.0, -0.0625*(3.0*eta+1.0)*(3.0*xi+1.0)*(xi-1.0)*(3.0*zeta-1.0) );
                      case 15:
                        return sign * RealGradient( 0.0, 0.0, 0.0625*(3.0*eta+1.0)*(3.0*xi+1.0)*(3.0*zeta+1.0)*(xi-1.0) );
                      case 16:
                        return sign * RealGradient( -0.0625*(zeta+1.0)*(3.0*eta-1.0)*(3.0*xi-1.0)*(3.0*zeta-1.0), 0.0, 0.0 );
                      case 17:
                        return sign * RealGradient( 0.0625*(zeta+1.0)*(3.0*xi+1.0)*(3.0*eta-1.0)*(3.0*zeta-1.0), 0.0, 0.0 );
                      case 18:
                        return sign * RealGradient( 0.0, -0.09375*(xi+1.0)*(zeta+1.0)*(3.0*xi-1.0)*(3.0*zeta-1.0), 0.0 );
                      case 19:
                        return sign * RealGradient( 0.0, 0.09375*(xi+1.0)*(zeta+1.0)*(3.0*xi-1.0)*(3.0*zeta-1.0), 0.0 );
                      case 20:
                        return sign * RealGradient( -0.0625*(zeta+1.0)*(3.0*eta+1.0)*(3.0*xi+1.0)*(3.0*zeta-1.0), 0.0, 0.0 );
                      case 21:
                        return sign * RealGradient( 0.0625*(zeta+1.0)*(3.0*eta+1.0)*(3.0*xi-1.0)*(3.0*zeta-1.0), 0.0, 0.0 );
                      case 22:
                        return sign * RealGradient( 0.0, -0.09375*(zeta+1.0)*(3.0*xi+1.0)*(xi-1.0)*(3.0*zeta-1.0), 0.0 );
                      case 23:
                        return sign * RealGradient( 0.0, 0.09375*(zeta+1.0)*(3.0*xi+1.0)*(xi-1.0)*(3.0*zeta-1.0), 0.0 );
                      case 24:
                      case 25:
                      case 26:
                      case 27:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 0, i - 24);
                          return c[0] * RealGradient( 0.0, 0.5625*(xi+1.0)*(3.0*zeta+1.0)*(xi-1.0)*(zeta-1.0), 0.0 ) +
                                 c[1] * RealGradient( 0.0, -0.5625*(xi+1.0)*(3.0*zeta+1.0)*(xi-1.0)*(zeta-1.0), 0.0 ) +
                                 c[2] * RealGradient( 0.375*eta*(3.0*zeta+1.0)*(zeta-1.0)*(3.0*xi-1.0), 0.0, 0.0 ) +
                                 c[3] * RealGradient( -0.375*eta*(3.0*xi+1.0)*(3.0*zeta+1.0)*(zeta-1.0), 0.0, 0.0 );
                        }
                      case 28:
                      case 29:
                      case 30:
                      case 31:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 1, i - 28);
                          return c[0] * RealGradient( 0.375*(zeta+1.0)*(zeta-1.0)*(3.0*eta-1.0)*(3.0*xi-1.0), 0.0, 0.0 ) +
                                 c[1] * RealGradient( -0.375*(zeta+1.0)*(3.0*xi+1.0)*(zeta-1.0)*(3.0*eta-1.0), 0.0, 0.0 ) +
                                 c[2] * RealGradient( 0.0, 0.0, 0.375*(xi+1.0)*(xi-1.0)*(3.0*eta-1.0)*(3.0*zeta-1.0) ) +
                                 c[3] * RealGradient( 0.0, 0.0, -0.375*(xi+1.0)*(3.0*zeta+1.0)*(xi-1.0)*(3.0*eta-1.0) );
                        }
                      case 32:
                      case 33:
                      case 34:
                      case 35:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 2, i - 32);
                          return c[0] * RealGradient( 0.0, 0.5625*(xi+1.0)*(zeta+1.0)*(zeta-1.0)*(3.0*xi-1.0), 0.0 ) +
                                 c[1] * RealGradient( 0.0, -0.5625*(xi+1.0)*(zeta+1.0)*(zeta-1.0)*(3.0*xi-1.0), 0.0 ) +
                                 c[2] * RealGradient( 0.0, 0.0, 0.375*eta*(xi+1.0)*(3.0*xi-1.0)*(3.0*zeta-1.0) ) +
                                 c[3] * RealGradient( 0.0, 0.0, -0.375*eta*(xi+1.0)*(3.0*zeta+1.0)*(3.0*xi-1.0) );
                        }
                      case 36:
                      case 37:
                      case 38:
                      case 39:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 3, i - 36);
                          return c[0] * RealGradient( 0.375*(zeta+1.0)*(3.0*eta+1.0)*(3.0*xi+1.0)*(zeta-1.0), 0.0, 0.0 ) +
                                 c[1] * RealGradient( -0.375*(zeta+1.0)*(3.0*eta+1.0)*(zeta-1.0)*(3.0*xi-1.0), 0.0, 0.0 ) +
                                 c[2] * RealGradient( 0.0, 0.0, 0.375*(xi+1.0)*(3.0*eta+1.0)*(xi-1.0)*(3.0*zeta-1.0) ) +
                                 c[3] * RealGradient( 0.0, 0.0, -0.375*(xi+1.0)*(3.0*eta+1.0)*(3.0*zeta+1.0)*(xi-1.0) );
                        }
                      case 40:
                      case 41:
                      case 42:
                      case 43:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 4, i - 40);
                          return c[0] * RealGradient( 0.0, 0.5625*(zeta+1.0)*(3.0*xi+1.0)*(xi-1.0)*(zeta-1.0), 0.0 ) +
                                 c[1] * RealGradient( 0.0, -0.5625*(zeta+1.0)*(3.0*xi+1.0)*(xi-1.0)*(zeta-1.0), 0.0 ) +
                                 c[2] * RealGradient( 0.0, 0.0, 0.375*eta*(3.0*xi+1.0)*(xi-1.0)*(3.0*zeta-1.0) ) +
                                 c[3] * RealGradient( 0.0, 0.0, -0.375*eta*(3.0*xi+1.0)*(3.0*zeta+1.0)*(xi-1.0) );
                        }
                      case 44:
                      case 45:
                      case 46:
                      case 47:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 5, i - 44);
                          return c[0] * RealGradient( 0.375*eta*(zeta+1.0)*(3.0*xi-1.0)*(3.0*zeta-1.0), 0.0, 0.0 ) +
                                 c[1] * RealGradient( -0.375*eta*(zeta+1.0)*(3.0*xi+1.0)*(3.0*zeta-1.0), 0.0, 0.0 ) +
                                 c[2] * RealGradient( 0.0, 0.5625*(xi+1.0)*(zeta+1.0)*(xi-1.0)*(3.0*zeta-1.0), 0.0 ) +
                                 c[3] * RealGradient( 0.0, -0.5625*(xi+1.0)*(zeta+1.0)*(xi-1.0)*(3.0*zeta-1.0), 0.0 );
                        }
                      case 48:
                        return RealGradient( -0.5625*eta*(zeta+1.0)*(zeta-1.0)*(3.0*xi-1.0), 0.0, 0.0 );
                      case 49:
                        return RealGradient( 0.5625*eta*(zeta+1.0)*(3.0*xi+1.0)*(zeta-1.0), 0.0, 0.0 );
                      case 50:
                        return RealGradient( 0.0, -0.84375*(xi+1.0)*(zeta+1.0)*(xi-1.0)*(zeta-1.0), 0.0 );
                      case 51:
                        return RealGradient( 0.0, 0.84375*(xi+1.0)*(zeta+1.0)*(xi-1.0)*(zeta-1.0), 0.0 );
                      case 52:
                        return RealGradient( 0.0, 0.0, -0.5625*eta*(xi+1.0)*(xi-1.0)*(3.0*zeta-1.0) );
                      case 53:
                        return RealGradient( 0.0, 0.0, 0.5625*eta*(xi+1.0)*(3.0*zeta+1.0)*(xi-1.0) );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      } // switch(ii)
                  } // j = 1

                  // d()/dzeta
                case 2:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient( -0.0625*(3.0*eta+1.0)*(eta-1.0)*(3.0*xi-1.0)*(3.0*zeta-1.0), 0.0, 0.0 );
                      case 1:
                        return sign * RealGradient( 0.0625*(3.0*eta+1.0)*(3.0*xi+1.0)*(eta-1.0)*(3.0*zeta-1.0), 0.0, 0.0 );
                      case 2:
                        return sign * RealGradient( 0.0, -0.0625*(xi+1.0)*(3.0*eta-1.0)*(3.0*xi-1.0)*(3.0*zeta-1.0), 0.0 );
                      case 3:
                        return sign * RealGradient( 0.0, 0.0625*(xi+1.0)*(3.0*eta+1.0)*(3.0*xi-1.0)*(3.0*zeta-1.0), 0.0 );
                      case 4:
                        return sign * RealGradient( -0.0625*(eta+1.0)*(3.0*xi+1.0)*(3.0*eta-1.0)*(3.0*zeta-1.0), 0.0, 0.0 );
                      case 5:
                        return sign * RealGradient( 0.0625*(eta+1.0)*(3.0*eta-1.0)*(3.0*xi-1.0)*(3.0*zeta-1.0), 0.0, 0.0 );
                      case 6:
                        return sign * RealGradient( 0.0, -0.0625*(3.0*xi+1.0)*(xi-1.0)*(3.0*eta-1.0)*(3.0*zeta-1.0), 0.0 );
                      case 7:
                        return sign * RealGradient( 0.0, 0.0625*(3.0*eta+1.0)*(3.0*xi+1.0)*(xi-1.0)*(3.0*zeta-1.0), 0.0 );
                      case 8:
                        return sign * RealGradient( 0.0, 0.0, -0.09375*(3.0*eta+1.0)*(3.0*xi+1.0)*(eta-1.0)*(xi-1.0) );
                      case 9:
                        return sign * RealGradient( 0.0, 0.0, 0.09375*(3.0*eta+1.0)*(3.0*xi+1.0)*(eta-1.0)*(xi-1.0) );
                      case 10:
                        return sign * RealGradient( 0.0, 0.0, -0.09375*(xi+1.0)*(3.0*eta+1.0)*(eta-1.0)*(3.0*xi-1.0) );
                      case 11:
                        return sign * RealGradient( 0.0, 0.0, 0.09375*(xi+1.0)*(3.0*eta+1.0)*(eta-1.0)*(3.0*xi-1.0) );
                      case 12:
                        return sign * RealGradient( 0.0, 0.0, -0.09375*(eta+1.0)*(xi+1.0)*(3.0*eta-1.0)*(3.0*xi-1.0) );
                      case 13:
                        return sign * RealGradient( 0.0, 0.0, 0.09375*(eta+1.0)*(xi+1.0)*(3.0*eta-1.0)*(3.0*xi-1.0) );
                      case 14:
                        return sign * RealGradient( 0.0, 0.0, -0.09375*(eta+1.0)*(3.0*xi+1.0)*(xi-1.0)*(3.0*eta-1.0) );
                      case 15:
                        return sign * RealGradient( 0.0, 0.0, 0.09375*(eta+1.0)*(3.0*xi+1.0)*(xi-1.0)*(3.0*eta-1.0) );
                      case 16:
                        return sign * RealGradient( -0.0625*(3.0*eta+1.0)*(3.0*zeta+1.0)*(eta-1.0)*(3.0*xi-1.0), 0.0, 0.0 );
                      case 17:
                        return sign * RealGradient( 0.0625*(3.0*eta+1.0)*(3.0*xi+1.0)*(3.0*zeta+1.0)*(eta-1.0), 0.0, 0.0 );
                      case 18:
                        return sign * RealGradient( 0.0, -0.0625*(xi+1.0)*(3.0*zeta+1.0)*(3.0*eta-1.0)*(3.0*xi-1.0), 0.0 );
                      case 19:
                        return sign * RealGradient( 0.0, 0.0625*(xi+1.0)*(3.0*eta+1.0)*(3.0*zeta+1.0)*(3.0*xi-1.0), 0.0 );
                      case 20:
                        return sign * RealGradient( -0.0625*(eta+1.0)*(3.0*xi+1.0)*(3.0*zeta+1.0)*(3.0*eta-1.0), 0.0, 0.0 );
                      case 21:
                        return sign * RealGradient( 0.0625*(eta+1.0)*(3.0*zeta+1.0)*(3.0*eta-1.0)*(3.0*xi-1.0), 0.0, 0.0 );
                      case 22:
                        return sign * RealGradient( 0.0, -0.0625*(3.0*xi+1.0)*(3.0*zeta+1.0)*(xi-1.0)*(3.0*eta-1.0), 0.0 );
                      case 23:
                        return sign * RealGradient( 0.0, 0.0625*(3.0*eta+1.0)*(3.0*xi+1.0)*(3.0*zeta+1.0)*(xi-1.0), 0.0 );
                      case 24:
                      case 25:
                      case 26:
                      case 27:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 0, i - 24);
                          return c[0] * RealGradient( 0.0, 0.375*(xi+1.0)*(xi-1.0)*(3.0*eta-1.0)*(3.0*zeta-1.0), 0.0 ) +
                                 c[1] * RealGradient( 0.0, -0.375*(xi+1.0)*(3.0*eta+1.0)*(xi-1.0)*(3.0*zeta-1.0), 0.0 ) +
                                 c[2] * RealGradient( 0.375*(eta+1.0)*(eta-1.0)*(3.0*xi-1.0)*(3.0*zeta-1.0), 0.0, 0.0 ) +
                                 c[3] * RealGradient( -0.375*(eta+1.0)*(3.0*xi+1.0)*(eta-1.0)*(3.0*zeta-1.0), 0.0, 0.0 );
                        }
                      case 28:
                      case 29:
                      case 30:
                      case 31:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 1, i - 28);
                          return c[0] * RealGradient( 0.375*zeta*(3.0*eta+1.0)*(eta-1.0)*(3.0*xi-1.0), 0.0, 0.0 ) +
                                 c[1] * RealGradient( -0.375*zeta*(3.0*eta+1.0)*(3.0*xi+1.0)*(eta-1.0), 0.0, 0.0 ) +
                                 c[2] * RealGradient( 0.0, 0.0, 0.5625*(xi+1.0)*(3.0*eta+1.0)*(eta-1.0)*(xi-1.0) ) +
                                 c[3] * RealGradient( 0.0, 0.0, -0.5625*(xi+1.0)*(3.0*eta+1.0)*(eta-1.0)*(xi-1.0) );
                        }
                      case 32:
                      case 33:
                      case 34:
                      case 35:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 2, i - 32);
                          return c[0] * RealGradient( 0.0, 0.375*zeta*(xi+1.0)*(3.0*eta-1.0)*(3.0*xi-1.0), 0.0 ) +
                                 c[1] * RealGradient( 0.0, -0.375*zeta*(xi+1.0)*(3.0*eta+1.0)*(3.0*xi-1.0), 0.0 ) +
                                 c[2] * RealGradient( 0.0, 0.0, 0.5625*(eta+1.0)*(xi+1.0)*(eta-1.0)*(3.0*xi-1.0) ) +
                                 c[3] * RealGradient( 0.0, 0.0, -0.5625*(eta+1.0)*(xi+1.0)*(eta-1.0)*(3.0*xi-1.0) );
                        }
                      case 36:
                      case 37:
                      case 38:
                      case 39:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 3, i - 36);
                          return c[0] * RealGradient( 0.375*zeta*(eta+1.0)*(3.0*xi+1.0)*(3.0*eta-1.0), 0.0, 0.0 ) +
                                 c[1] * RealGradient( -0.375*zeta*(eta+1.0)*(3.0*eta-1.0)*(3.0*xi-1.0), 0.0, 0.0 ) +
                                 c[2] * RealGradient( 0.0, 0.0, 0.5625*(eta+1.0)*(xi+1.0)*(xi-1.0)*(3.0*eta-1.0) ) +
                                 c[3] * RealGradient( 0.0, 0.0, -0.5625*(eta+1.0)*(xi+1.0)*(xi-1.0)*(3.0*eta-1.0) );
                        }
                      case 40:
                      case 41:
                      case 42:
                      case 43:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 4, i - 40);
                          return c[0] * RealGradient( 0.0, 0.375*zeta*(3.0*eta+1.0)*(3.0*xi+1.0)*(xi-1.0), 0.0 ) +
                                 c[1] * RealGradient( 0.0, -0.375*zeta*(3.0*xi+1.0)*(xi-1.0)*(3.0*eta-1.0), 0.0 ) +
                                 c[2] * RealGradient( 0.0, 0.0, 0.5625*(eta+1.0)*(3.0*xi+1.0)*(eta-1.0)*(xi-1.0) ) +
                                 c[3] * RealGradient( 0.0, 0.0, -0.5625*(eta+1.0)*(3.0*xi+1.0)*(eta-1.0)*(xi-1.0) );
                        }
                      case 44:
                      case 45:
                      case 46:
                      case 47:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 5, i - 44);
                          return c[0] * RealGradient( 0.375*(eta+1.0)*(3.0*zeta+1.0)*(eta-1.0)*(3.0*xi-1.0), 0.0, 0.0 ) +
                                 c[1] * RealGradient( -0.375*(eta+1.0)*(3.0*xi+1.0)*(3.0*zeta+1.0)*(eta-1.0), 0.0, 0.0 ) +
                                 c[2] * RealGradient( 0.0, 0.375*(xi+1.0)*(3.0*zeta+1.0)*(xi-1.0)*(3.0*eta-1.0), 0.0 ) +
                                 c[3] * RealGradient( 0.0, -0.375*(xi+1.0)*(3.0*eta+1.0)*(3.0*zeta+1.0)*(xi-1.0), 0.0 );
                        }
                      case 48:
                        return RealGradient( -0.5625*zeta*(eta+1.0)*(eta-1.0)*(3.0*xi-1.0), 0.0, 0.0 );
                      case 49:
                        return RealGradient( 0.5625*zeta*(eta+1.0)*(3.0*xi+1.0)*(eta-1.0), 0.0, 0.0 );
                      case 50:
                        return RealGradient( 0.0, -0.5625*zeta*(xi+1.0)*(xi-1.0)*(3.0*eta-1.0), 0.0 );
                      case 51:
                        return RealGradient( 0.0, 0.5625*zeta*(xi+1.0)*(3.0*eta+1.0)*(xi-1.0), 0.0 );
                      case 52:
                        return RealGradient( 0.0, 0.0, -0.84375*(eta+1.0)*(xi+1.0)*(eta-1.0)*(xi-1.0) );
                      case 53:
                        return RealGradient( 0.0, 0.0, 0.84375*(eta+1.0)*(xi+1.0)*(eta-1.0)*(xi-1.0) );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      } // switch(ii)
                  } // j = 2

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
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
              switch (j)
                {
                  // d()/dxi
                case 0:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient( 0.5*(3.0*zeta+1.0)*(zeta-1.0)*(4.0*eta-3.0), -0.5*(3.0*zeta+1.0)*(zeta-1.0)*(4.0*eta+8.0*xi-3.0), 0.0 );
                      case 1:
                        return sign * RealGradient( -0.5*(3.0*zeta+1.0)*(zeta-1.0)*(4.0*eta-3.0), (3.0*zeta+1.0)*(zeta-1.0)*(4.0*xi-1.0), 0.0 );
                      case 2:
                        return sign * RealGradient( -2.0*eta*(3.0*zeta+1.0)*(zeta-1.0), (3.0*zeta+1.0)*(zeta-1.0)*(4.0*xi-1.0), 0.0 );
                      case 3:
                        return sign * RealGradient( 0.0, 0.5*(3.0*zeta+1.0)*(zeta-1.0)*(4.0*eta-1.0), 0.0 );
                      case 4:
                        return sign * RealGradient( -2.0*eta*(3.0*zeta+1.0)*(zeta-1.0), (3.0*zeta+1.0)*(zeta-1.0)*(2.0*eta+4.0*xi-3.0), 0.0 );
                      case 5:
                        return sign * RealGradient( 0.0, -0.5*(3.0*zeta+1.0)*(zeta-1.0)*(4.0*eta-1.0), 0.0 );
                      case 6:
                        return sign * RealGradient( 0.0, 0.0, -(3.0*zeta-1.0)*(3.0*eta+3.0*xi-2.0) );
                      case 7:
                        return sign * RealGradient( 0.0, 0.0, (3.0*zeta+1.0)*(3.0*eta+3.0*xi-2.0) );
                      case 8:
                        return sign * RealGradient( 0.0, 0.0, -(3.0*xi-1.0)*(3.0*zeta-1.0) );
                      case 9:
                        return sign * RealGradient( 0.0, 0.0, (3.0*zeta+1.0)*(3.0*xi-1.0) );
                      case 12:
                        return sign * RealGradient( 0.5*(zeta+1.0)*(3.0*zeta-1.0)*(4.0*eta-3.0), -0.5*(zeta+1.0)*(3.0*zeta-1.0)*(4.0*eta+8.0*xi-3.0), 0.0 );
                      case 13:
                        return sign * RealGradient( -0.5*(zeta+1.0)*(3.0*zeta-1.0)*(4.0*eta-3.0), (zeta+1.0)*(3.0*zeta-1.0)*(4.0*xi-1.0), 0.0 );
                      case 14:
                        return sign * RealGradient( -2.0*eta*(zeta+1.0)*(3.0*zeta-1.0), (zeta+1.0)*(3.0*zeta-1.0)*(4.0*xi-1.0), 0.0 );
                      case 15:
                        return sign * RealGradient( 0.0, 0.5*(zeta+1.0)*(3.0*zeta-1.0)*(4.0*eta-1.0), 0.0 );
                      case 16:
                        return sign * RealGradient( -2.0*eta*(zeta+1.0)*(3.0*zeta-1.0), (zeta+1.0)*(3.0*zeta-1.0)*(2.0*eta+4.0*xi-3.0), 0.0 );
                      case 17:
                        return sign * RealGradient( 0.0, -0.5*(zeta+1.0)*(3.0*zeta-1.0)*(4.0*eta-1.0), 0.0 );
                      case 18:
                      case 19:
                      case 20:
                      case 21:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 1, i - 18);
                          return c[0] * RealGradient( -3.0*(zeta+1.0)*(zeta-1.0)*(4.0*eta-3.0), 3.0*(zeta+1.0)*(zeta-1.0)*(4.0*eta+8.0*xi-3.0), 0.0 ) +
                                 c[1] * RealGradient( 3.0*(zeta+1.0)*(zeta-1.0)*(4.0*eta-3.0), -6.0*(zeta+1.0)*(zeta-1.0)*(4.0*xi-1.0), 0.0 ) +
                                 c[2] * RealGradient( 0.0, 0.0, 3.0*(3.0*zeta-1.0)*(eta+2.0*xi-1.0) ) +
                                 c[3] * RealGradient( 0.0, 0.0, -3.0*(3.0*zeta+1.0)*(eta+2.0*xi-1.0) );
                        }
                      case 22:
                      case 23:
                      case 24:
                      case 25:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 2, i - 22);
                          return c[0] * RealGradient( 12.0*eta*(zeta+1.0)*(zeta-1.0), -6.0*(zeta+1.0)*(zeta-1.0)*(4.0*xi-1.0), 0.0 ) +
                                 c[1] * RealGradient( 0.0, -3.0*(zeta+1.0)*(zeta-1.0)*(4.0*eta-1.0), 0.0 ) +
                                 c[2] * RealGradient( 0.0, 0.0, -3.0*eta*(3.0*zeta-1.0) ) +
                                 c[3] * RealGradient( 0.0, 0.0, 3.0*eta*(3.0*zeta+1.0) );
                        }
                      case 26:
                      case 27:
                      case 28:
                      case 29:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 3, i - 26);
                          return c[0] * RealGradient( 0.0, -3.0*(zeta+1.0)*(zeta-1.0)*(4.0*eta-1.0), 0.0 ) +
                                 c[1] * RealGradient( -12.0*eta*(zeta+1.0)*(zeta-1.0), 6.0*(zeta+1.0)*(zeta-1.0)*(2.0*eta+4.0*xi-3.0), 0.0 ) +
                                 c[2] * RealGradient( 0.0, 0.0, 3.0*eta*(3.0*zeta-1.0) ) +
                                 c[3] * RealGradient( 0.0, 0.0, -3.0*eta*(3.0*zeta+1.0) );
                        }
                      case 30:
                      case 31:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 0, i - 30);
                          return c[0] * RealGradient( 4.0*eta*(3.0*zeta+1.0)*(zeta-1.0), -2.0*(3.0*zeta+1.0)*(zeta-1.0)*(eta+4.0*xi-2.0), 0.0 ) +
                                 c[1] * RealGradient( -2.0*eta*(3.0*zeta+1.0)*(zeta-1.0), 2.0*(3.0*zeta+1.0)*(zeta-1.0)*(2.0*eta+2.0*xi-1.0), 0.0 );
                        }
                      case 32:
                      case 33:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 4, i - 32);
                          return c[0] * RealGradient( -2.0*eta*(zeta+1.0)*(3.0*zeta-1.0), 2.0*(zeta+1.0)*(3.0*zeta-1.0)*(2.0*eta+2.0*xi-1.0), 0.0 ) +
                                 c[1] * RealGradient( 4.0*eta*(zeta+1.0)*(3.0*zeta-1.0), -2.0*(zeta+1.0)*(3.0*zeta-1.0)*(eta+4.0*xi-2.0), 0.0 );
                        }
                      case 34:
                        return RealGradient( 6.0*eta*(zeta+1.0)*(zeta-1.0), -6.0*(zeta+1.0)*(zeta-1.0)*(2.0*eta+2.0*xi-1.0), 0.0 );
                      case 35:
                        return RealGradient( -12.0*eta*(zeta+1.0)*(zeta-1.0), 6.0*(zeta+1.0)*(zeta-1.0)*(eta+4.0*xi-2.0), 0.0 );
                      case 10:
                      case 11:
                        return RealGradient();

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      } // switch(ii)
                  } // j = 0

                  // d()/deta
                case 1:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient( (3.0*zeta+1.0)*(zeta-1.0)*(4.0*eta+2.0*xi-3.0), -2.0*xi*(3.0*zeta+1.0)*(zeta-1.0), 0.0 );
                      case 1:
                        return sign * RealGradient( -0.5*(3.0*zeta+1.0)*(zeta-1.0)*(4.0*xi-1.0), 0.0, 0.0 );
                      case 2:
                        return sign * RealGradient( -0.5*(3.0*zeta+1.0)*(zeta-1.0)*(4.0*xi-1.0), 0.0, 0.0 );
                      case 3:
                        return sign * RealGradient( -(3.0*zeta+1.0)*(zeta-1.0)*(4.0*eta-1.0), 2.0*xi*(3.0*zeta+1.0)*(zeta-1.0), 0.0 );
                      case 4:
                        return sign * RealGradient( -0.5*(3.0*zeta+1.0)*(zeta-1.0)*(8.0*eta+4.0*xi-3.0), 0.5*(3.0*zeta+1.0)*(zeta-1.0)*(4.0*xi-3.0), 0.0 );
                      case 5:
                        return sign * RealGradient( (3.0*zeta+1.0)*(zeta-1.0)*(4.0*eta-1.0), -0.5*(3.0*zeta+1.0)*(zeta-1.0)*(4.0*xi-3.0), 0.0 );
                      case 6:
                        return sign * RealGradient( 0.0, 0.0, -(3.0*zeta-1.0)*(3.0*eta+3.0*xi-2.0) );
                      case 7:
                        return sign * RealGradient( 0.0, 0.0, (3.0*zeta+1.0)*(3.0*eta+3.0*xi-2.0) );
                      case 10:
                        return sign * RealGradient( 0.0, 0.0, -(3.0*eta-1.0)*(3.0*zeta-1.0) );
                      case 11:
                        return sign * RealGradient( 0.0, 0.0, (3.0*zeta+1.0)*(3.0*eta-1.0) );
                      case 12:
                        return sign * RealGradient( (zeta+1.0)*(3.0*zeta-1.0)*(4.0*eta+2.0*xi-3.0), -2.0*xi*(zeta+1.0)*(3.0*zeta-1.0), 0.0 );
                      case 13:
                        return sign * RealGradient( -0.5*(zeta+1.0)*(3.0*zeta-1.0)*(4.0*xi-1.0), 0.0, 0.0 );
                      case 14:
                        return sign * RealGradient( -0.5*(zeta+1.0)*(3.0*zeta-1.0)*(4.0*xi-1.0), 0.0, 0.0 );
                      case 15:
                        return sign * RealGradient( -(zeta+1.0)*(3.0*zeta-1.0)*(4.0*eta-1.0), 2.0*xi*(zeta+1.0)*(3.0*zeta-1.0), 0.0 );
                      case 16:
                        return sign * RealGradient( -0.5*(zeta+1.0)*(3.0*zeta-1.0)*(8.0*eta+4.0*xi-3.0), 0.5*(zeta+1.0)*(3.0*zeta-1.0)*(4.0*xi-3.0), 0.0 );
                      case 17:
                        return sign * RealGradient( (zeta+1.0)*(3.0*zeta-1.0)*(4.0*eta-1.0), -0.5*(zeta+1.0)*(3.0*zeta-1.0)*(4.0*xi-3.0), 0.0 );
                      case 18:
                      case 19:
                      case 20:
                      case 21:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 1, i - 18);
                          return c[0] * RealGradient( -6.0*(zeta+1.0)*(zeta-1.0)*(4.0*eta+2.0*xi-3.0), 12.0*xi*(zeta+1.0)*(zeta-1.0), 0.0 ) +
                                 c[1] * RealGradient( 3.0*(zeta+1.0)*(zeta-1.0)*(4.0*xi-1.0), 0.0, 0.0 ) +
                                 c[2] * RealGradient( 0.0, 0.0, 3.0*xi*(3.0*zeta-1.0) ) +
                                 c[3] * RealGradient( 0.0, 0.0, -3.0*xi*(3.0*zeta+1.0) );
                        }
                      case 22:
                      case 23:
                      case 24:
                      case 25:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 2, i - 22);
                          return c[0] * RealGradient( 3.0*(zeta+1.0)*(zeta-1.0)*(4.0*xi-1.0), 0.0, 0.0 ) +
                                 c[1] * RealGradient( 6.0*(zeta+1.0)*(zeta-1.0)*(4.0*eta-1.0), -12.0*xi*(zeta+1.0)*(zeta-1.0), 0.0 ) +
                                 c[2] * RealGradient( 0.0, 0.0, -3.0*xi*(3.0*zeta-1.0) ) +
                                 c[3] * RealGradient( 0.0, 0.0, 3.0*xi*(3.0*zeta+1.0) );
                        }
                      case 26:
                      case 27:
                      case 28:
                      case 29:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 3, i - 26);
                          return c[0] * RealGradient( 6.0*(zeta+1.0)*(zeta-1.0)*(4.0*eta-1.0), -3.0*(zeta+1.0)*(zeta-1.0)*(4.0*xi-3.0), 0.0 ) +
                                 c[1] * RealGradient( -3.0*(zeta+1.0)*(zeta-1.0)*(8.0*eta+4.0*xi-3.0), 3.0*(zeta+1.0)*(zeta-1.0)*(4.0*xi-3.0), 0.0 ) +
                                 c[2] * RealGradient( 0.0, 0.0, 3.0*(3.0*zeta-1.0)*(2.0*eta+xi-1.0) ) +
                                 c[3] * RealGradient( 0.0, 0.0, -3.0*(3.0*zeta+1.0)*(2.0*eta+xi-1.0) );
                        }
                      case 30:
                      case 31:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 0, i - 30);
                          return c[0] * RealGradient( 2.0*(3.0*zeta+1.0)*(zeta-1.0)*(2.0*eta+2.0*xi-1.0), -2.0*xi*(3.0*zeta+1.0)*(zeta-1.0), 0.0 ) +
                                 c[1] * RealGradient( -2.0*(3.0*zeta+1.0)*(zeta-1.0)*(4.0*eta+xi-2.0), 4.0*xi*(3.0*zeta+1.0)*(zeta-1.0), 0.0 );
                        }
                      case 32:
                      case 33:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 4, i - 32);
                          return c[0] * RealGradient( -2.0*(zeta+1.0)*(3.0*zeta-1.0)*(4.0*eta+xi-2.0), 4.0*xi*(zeta+1.0)*(3.0*zeta-1.0), 0.0 ) +
                                 c[1] * RealGradient( 2.0*(zeta+1.0)*(3.0*zeta-1.0)*(2.0*eta+2.0*xi-1.0), -2.0*xi*(zeta+1.0)*(3.0*zeta-1.0), 0.0 );
                        }
                      case 34:
                        return RealGradient( 6.0*(zeta+1.0)*(zeta-1.0)*(4.0*eta+xi-2.0), -12.0*xi*(zeta+1.0)*(zeta-1.0), 0.0 );
                      case 35:
                        return RealGradient( -6.0*(zeta+1.0)*(zeta-1.0)*(2.0*eta+2.0*xi-1.0), 6.0*xi*(zeta+1.0)*(zeta-1.0), 0.0 );
                      case 8:
                      case 9:
                        return RealGradient();

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      } // switch(ii)
                  } // j = 1

                  // d()/dzeta
                case 2:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient( (3.0*zeta-1.0)*(4.0*eta*eta+4.0*eta*xi-6.0*eta-3.0*xi+2.0), -xi*(3.0*zeta-1.0)*(4.0*eta+4.0*xi-3.0), 0.0 );
                      case 1:
                        return sign * RealGradient( -(3.0*zeta-1.0)*(4.0*eta*xi-eta-3.0*xi+1.0), 2.0*xi*(2.0*xi-1.0)*(3.0*zeta-1.0), 0.0 );
                      case 2:
                        return sign * RealGradient( -eta*(3.0*zeta-1.0)*(4.0*xi-1.0), 2.0*xi*(2.0*xi-1.0)*(3.0*zeta-1.0), 0.0 );
                      case 3:
                        return sign * RealGradient( -2.0*eta*(2.0*eta-1.0)*(3.0*zeta-1.0), xi*(3.0*zeta-1.0)*(4.0*eta-1.0), 0.0 );
                      case 4:
                        return sign * RealGradient( -eta*(3.0*zeta-1.0)*(4.0*eta+4.0*xi-3.0), (3.0*zeta-1.0)*(4.0*eta*xi-3.0*eta+4.0*xi*xi-6.0*xi+2.0), 0.0 );
                      case 5:
                        return sign * RealGradient( 2.0*eta*(2.0*eta-1.0)*(3.0*zeta-1.0), -(3.0*zeta-1.0)*(4.0*eta*xi-3.0*eta-xi+1.0), 0.0 );
                      case 6:
                        return sign * RealGradient( 0.0, 0.0, -1.5*(eta+xi-1.0)*(3.0*eta+3.0*xi-1.0) );
                      case 7:
                        return sign * RealGradient( 0.0, 0.0, 1.5*(eta+xi-1.0)*(3.0*eta+3.0*xi-1.0) );
                      case 8:
                        return sign * RealGradient( 0.0, 0.0, -1.5*xi*(3.0*xi-2.0) );
                      case 9:
                        return sign * RealGradient( 0.0, 0.0, 1.5*xi*(3.0*xi-2.0) );
                      case 10:
                        return sign * RealGradient( 0.0, 0.0, -1.5*eta*(3.0*eta-2.0) );
                      case 11:
                        return sign * RealGradient( 0.0, 0.0, 1.5*eta*(3.0*eta-2.0) );
                      case 12:
                        return sign * RealGradient( (3.0*zeta+1.0)*(4.0*eta*eta+4.0*eta*xi-6.0*eta-3.0*xi+2.0), -xi*(3.0*zeta+1.0)*(4.0*eta+4.0*xi-3.0), 0.0 );
                      case 13:
                        return sign * RealGradient( -(3.0*zeta+1.0)*(4.0*eta*xi-eta-3.0*xi+1.0), 2.0*xi*(3.0*zeta+1.0)*(2.0*xi-1.0), 0.0 );
                      case 14:
                        return sign * RealGradient( -eta*(3.0*zeta+1.0)*(4.0*xi-1.0), 2.0*xi*(3.0*zeta+1.0)*(2.0*xi-1.0), 0.0 );
                      case 15:
                        return sign * RealGradient( -2.0*eta*(3.0*zeta+1.0)*(2.0*eta-1.0), xi*(3.0*zeta+1.0)*(4.0*eta-1.0), 0.0 );
                      case 16:
                        return sign * RealGradient( -eta*(3.0*zeta+1.0)*(4.0*eta+4.0*xi-3.0), (3.0*zeta+1.0)*(4.0*eta*xi-3.0*eta+4.0*xi*xi-6.0*xi+2.0), 0.0 );
                      case 17:
                        return sign * RealGradient( 2.0*eta*(3.0*zeta+1.0)*(2.0*eta-1.0), -(3.0*zeta+1.0)*(4.0*eta*xi-3.0*eta-xi+1.0), 0.0 );
                      case 18:
                      case 19:
                      case 20:
                      case 21:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 1, i - 18);
                          return c[0] * RealGradient( -6.0*zeta*(4.0*eta*eta+4.0*eta*xi-6.0*eta-3.0*xi+2.0), 6.0*xi*zeta*(4.0*eta+4.0*xi-3.0), 0.0 ) +
                                 c[1] * RealGradient( 6.0*zeta*(4.0*eta*xi-eta-3.0*xi+1.0), -12.0*xi*zeta*(2.0*xi-1.0), 0.0 ) +
                                 c[2] * RealGradient( 0.0, 0.0, 9.0*xi*(eta+xi-1.0) ) +
                                 c[3] * RealGradient( 0.0, 0.0, -9.0*xi*(eta+xi-1.0) );
                        }
                      case 22:
                      case 23:
                      case 24:
                      case 25:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 2, i - 22);
                          return c[0] * RealGradient( 6.0*eta*zeta*(4.0*xi-1.0), -12.0*xi*zeta*(2.0*xi-1.0), 0.0 ) +
                                 c[1] * RealGradient( 12.0*eta*zeta*(2.0*eta-1.0), -6.0*xi*zeta*(4.0*eta-1.0), 0.0 ) +
                                 c[2] * RealGradient( 0.0, 0.0, -9.0*eta*xi ) +
                                 c[3] * RealGradient( 0.0, 0.0, 9.0*eta*xi );
                        }
                      case 26:
                      case 27:
                      case 28:
                      case 29:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 3, i - 26);
                          return c[0] * RealGradient( 12.0*eta*zeta*(2.0*eta-1.0), -6.0*zeta*(4.0*eta*xi-3.0*eta-xi+1.0), 0.0 ) +
                                 c[1] * RealGradient( -6.0*eta*zeta*(4.0*eta+4.0*xi-3.0), 6.0*zeta*(4.0*eta*xi-3.0*eta+4.0*xi*xi-6.0*xi+2.0), 0.0 ) +
                                 c[2] * RealGradient( 0.0, 0.0, 9.0*eta*(eta+xi-1.0) ) +
                                 c[3] * RealGradient( 0.0, 0.0, -9.0*eta*(eta+xi-1.0) );
                        }
                      case 30:
                      case 31:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 0, i - 30);
                          return c[0] * RealGradient( 4.0*eta*(3.0*zeta-1.0)*(eta+2.0*xi-1.0), -4.0*xi*(3.0*zeta-1.0)*(eta+2.0*xi-2.0), 0.0 ) +
                                 c[1] * RealGradient( -4.0*eta*(3.0*zeta-1.0)*(2.0*eta+xi-2.0), 4.0*xi*(3.0*zeta-1.0)*(2.0*eta+xi-1.0), 0.0 );
                        }
                      case 32:
                      case 33:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 4, i - 32);
                          return c[0] * RealGradient( -4.0*eta*(3.0*zeta+1.0)*(2.0*eta+xi-2.0), 4.0*xi*(3.0*zeta+1.0)*(2.0*eta+xi-1.0), 0.0 ) +
                                 c[1] * RealGradient( 4.0*eta*(3.0*zeta+1.0)*(eta+2.0*xi-1.0), -4.0*xi*(3.0*zeta+1.0)*(eta+2.0*xi-2.0), 0.0 );
                        }
                      case 34:
                        return RealGradient( 12.0*eta*zeta*(2.0*eta+xi-2.0), -12.0*xi*zeta*(2.0*eta+xi-1.0), 0.0 );
                      case 35:
                        return RealGradient( -12.0*eta*zeta*(eta+2.0*xi-1.0), 12.0*xi*zeta*(eta+2.0*xi-2.0), 0.0 );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      } // switch(ii)
                  } // j = 2

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
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

              // The derivatives of the rational shape functions have no limit at
              // the apex: only keep them finite there.
              const Real den = std::max(1.0 - zeta, 1.e-6);

              switch (j)
                {
                  // d()/dxi
                case 0:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient( 0.125*(4.0*zeta-3.0)*(3.0*eta-zeta+1.0)*(eta+zeta-1.0)/(den*den), 0.0, -0.125*(3.0*eta-zeta+1.0)*(eta+zeta-1.0)*(8.0*xi+4.0*zeta-1.0)/(den*den) );
                      case 1:
                        return sign * RealGradient( -0.125*(4.0*zeta-3.0)*(3.0*eta-zeta+1.0)*(eta+zeta-1.0)/(den*den), 0.0, 0.125*(3.0*eta-zeta+1.0)*(8.0*xi-4.0*zeta+1.0)*(eta+zeta-1.0)/(den*den) );
                      case 2:
                        return sign * RealGradient( 0.0, 0.25*(3.0*xi-zeta+1.0)*(4.0*eta*zeta-3.0*eta+4.0*zeta*zeta-5.0*zeta+1.0)/(den*den), -0.25*(3.0*xi-zeta+1.0)*(4.0*eta*eta+4.0*eta*zeta-eta+zeta-1.0)/(den*den) );
                      case 3:
                        return sign * RealGradient( 0.0, -0.25*(3.0*xi-zeta+1.0)*(4.0*eta*zeta-3.0*eta-4.0*zeta*zeta+5.0*zeta-1.0)/(den*den), 0.25*(3.0*xi-zeta+1.0)*(4.0*eta*eta-4.0*eta*zeta+eta+zeta-1.0)/(den*den) );
                      case 4:
                        return sign * RealGradient( 0.125*(4.0*zeta-3.0)*(eta-zeta+1.0)*(3.0*eta+zeta-1.0)/(den*den), 0.0, -0.125*(eta-zeta+1.0)*(8.0*xi-4.0*zeta+1.0)*(3.0*eta+zeta-1.0)/(den*den) );
                      case 5:
                        return sign * RealGradient( -0.125*(4.0*zeta-3.0)*(eta-zeta+1.0)*(3.0*eta+zeta-1.0)/(den*den), 0.0, 0.125*(eta-zeta+1.0)*(3.0*eta+zeta-1.0)*(8.0*xi+4.0*zeta-1.0)/(den*den) );
                      case 6:
                        return sign * RealGradient( 0.0, 0.25*(3.0*xi+zeta-1.0)*(4.0*eta*zeta-3.0*eta+4.0*zeta*zeta-5.0*zeta+1.0)/(den*den), -0.25*(3.0*xi+zeta-1.0)*(4.0*eta*eta+4.0*eta*zeta-eta+zeta-1.0)/(den*den) );
                      case 7:
                        return sign * RealGradient( 0.0, -0.25*(3.0*xi+zeta-1.0)*(4.0*eta*zeta-3.0*eta-4.0*zeta*zeta+5.0*zeta-1.0)/(den*den), 0.25*(3.0*xi+zeta-1.0)*(4.0*eta*eta-4.0*eta*zeta+eta+zeta-1.0)/(den*den) );
                      case 8:
                        return sign * RealGradient( -0.125*zeta*(eta+zeta-1.0)*(3.0*eta*zeta+7.0*eta-11.0*zeta*zeta+12.0*zeta-1.0)/(den*den*den), -0.25*zeta*(3.0*eta*xi*zeta+7.0*eta*xi-4.0*eta*zeta*zeta+8.0*eta*zeta-4.0*eta-2.0*xi*zeta*zeta+4.0*xi*zeta-2.0*xi-4.0*zeta*zeta*zeta+10.0*zeta*zeta-8.0*zeta+2.0)/(den*den*den), 0.25*(8.0*eta+11.0*xi+3.0*zeta)-0.25*(6.0*xi*eta*eta*zeta*zeta+16.0*xi*zeta*eta*eta-7.0*xi*eta*eta-5.0*eta*eta*zeta*zeta*zeta+14.0*eta*eta*zeta*zeta-13.0*zeta*eta*eta+4.0*eta*eta-10.0*eta*xi*zeta*zeta*zeta+28.0*eta*xi*zeta*zeta-26.0*eta*xi*zeta+8.0*eta*xi-6.0*eta*zeta*zeta*zeta+16.0*eta*zeta*zeta-14.0*eta*zeta+4.0*eta-15.0*xi*zeta*zeta*zeta+40.0*xi*zeta*zeta-35.0*xi*zeta+10.0*xi)/(den*den*den*den) );
                      case 9:
                        return sign * RealGradient( 1.25*zeta*(3.0*zeta-1.0)*(eta-zeta+1.0)*(eta+zeta-1.0)/(den*den*den), 0.5*zeta*(15.0*eta*xi*zeta-5.0*eta*xi-2.0*xi*zeta*zeta+4.0*xi*zeta-2.0*xi-4.0*zeta*zeta*zeta+10.0*zeta*zeta-8.0*zeta+2.0)/(den*den*den), -0.5*(8.0*eta+15.0*xi+3.0*zeta-1.0)+0.25*(60.0*xi*eta*eta*zeta*zeta-37.0*xi*zeta*eta*eta+7.0*xi*eta*eta-2.0*eta*eta*zeta*zeta*zeta+4.0*eta*eta*zeta*zeta-2.0*zeta*eta*eta-4.0*eta*xi*zeta*zeta*zeta+8.0*eta*xi*zeta*zeta-4.0*eta*xi*zeta-18.0*eta*zeta*zeta*zeta+50.0*eta*zeta*zeta-46.0*eta*zeta+14.0*eta-33.0*xi*zeta*zeta*zeta+89.0*xi*zeta*zeta-79.0*xi*zeta+23.0*xi)/(den*den*den*den) );
                      case 10:
                        return sign * RealGradient( -0.125*zeta*(eta+zeta-1.0)*(3.0*eta*zeta+7.0*eta-11.0*zeta*zeta+12.0*zeta-1.0)/(den*den*den), -0.25*zeta*(3.0*eta*xi*zeta+7.0*eta*xi+4.0*eta*zeta*zeta-8.0*eta*zeta+4.0*eta-2.0*xi*zeta*zeta+4.0*xi*zeta-2.0*xi+4.0*zeta*zeta*zeta-10.0*zeta*zeta+8.0*zeta-2.0)/(den*den*den), -0.25*(8.0*eta-11.0*xi+3.0*zeta)-0.25*(6.0*xi*eta*eta*zeta*zeta+16.0*xi*zeta*eta*eta-7.0*xi*eta*eta+5.0*eta*eta*zeta*zeta*zeta-14.0*eta*eta*zeta*zeta+13.0*zeta*eta*eta-4.0*eta*eta-10.0*eta*xi*zeta*zeta*zeta+28.0*eta*xi*zeta*zeta-26.0*eta*xi*zeta+8.0*eta*xi+6.0*eta*zeta*zeta*zeta-16.0*eta*zeta*zeta+14.0*eta*zeta-4.0*eta-15.0*xi*zeta*zeta*zeta+40.0*xi*zeta*zeta-35.0*xi*zeta+10.0*xi)/(den*den*den*den) );
                      case 11:
                        return sign * RealGradient( 1.25*zeta*(3.0*zeta-1.0)*(eta-zeta+1.0)*(eta+zeta-1.0)/(den*den*den), 0.5*zeta*(15.0*eta*xi*zeta-5.0*eta*xi-2.0*xi*zeta*zeta+4.0*xi*zeta-2.0*xi+4.0*zeta*zeta*zeta-10.0*zeta*zeta+8.0*zeta-2.0)/(den*den*den), 0.5*(8.0*eta-15.0*xi+3.0*zeta-1.0)+0.25*(60.0*xi*eta*eta*zeta*zeta-37.0*xi*zeta*eta*eta+7.0*xi*eta*eta+2.0*eta*eta*zeta*zeta*zeta-4.0*eta*eta*zeta*zeta+2.0*zeta*eta*eta-4.0*eta*xi*zeta*zeta*zeta+8.0*eta*xi*zeta*zeta-4.0*eta*xi*zeta+18.0*eta*zeta*zeta*zeta-50.0*eta*zeta*zeta+46.0*eta*zeta-14.0*eta-33.0*xi*zeta*zeta*zeta+89.0*xi*zeta*zeta-79.0*xi*zeta+23.0*xi)/(den*den*den*den) );
                      case 12:
                        return sign * RealGradient( -0.125*zeta*(eta-zeta+1.0)*(3.0*eta*zeta+7.0*eta+11.0*zeta*zeta-12.0*zeta+1.0)/(den*den*den), -0.25*zeta*(3.0*eta*xi*zeta+7.0*eta*xi+4.0*eta*zeta*zeta-8.0*eta*zeta+4.0*eta+2.0*xi*zeta*zeta-4.0*xi*zeta+2.0*xi-4.0*zeta*zeta*zeta+10.0*zeta*zeta-8.0*zeta+2.0)/(den*den*den), 0.25*(8.0*eta+11.0*xi-3.0*zeta)-0.25*(6.0*xi*eta*eta*zeta*zeta+16.0*xi*zeta*eta*eta-7.0*xi*eta*eta+5.0*eta*eta*zeta*zeta*zeta-14.0*eta*eta*zeta*zeta+13.0*zeta*eta*eta-4.0*eta*eta+10.0*eta*xi*zeta*zeta*zeta-28.0*eta*xi*zeta*zeta+26.0*eta*xi*zeta-8.0*eta*xi-6.0*eta*zeta*zeta*zeta+16.0*eta*zeta*zeta-14.0*eta*zeta+4.0*eta-15.0*xi*zeta*zeta*zeta+40.0*xi*zeta*zeta-35.0*xi*zeta+10.0*xi)/(den*den*den*den) );
                      case 13:
                        return sign * RealGradient( 1.25*zeta*(3.0*zeta-1.0)*(eta-zeta+1.0)*(eta+zeta-1.0)/(den*den*den), 0.5*zeta*(15.0*eta*xi*zeta-5.0*eta*xi+2.0*xi*zeta*zeta-4.0*xi*zeta+2.0*xi-4.0*zeta*zeta*zeta+10.0*zeta*zeta-8.0*zeta+2.0)/(den*den*den), -0.5*(8.0*eta+15.0*xi-3.0*zeta+1.0)+0.25*(60.0*xi*eta*eta*zeta*zeta-37.0*xi*zeta*eta*eta+7.0*xi*eta*eta+2.0*eta*eta*zeta*zeta*zeta-4.0*eta*eta*zeta*zeta+2.0*zeta*eta*eta+4.0*eta*xi*zeta*zeta*zeta-8.0*eta*xi*zeta*zeta+4.0*eta*xi*zeta-18.0*eta*zeta*zeta*zeta+50.0*eta*zeta*zeta-46.0*eta*zeta+14.0*eta-33.0*xi*zeta*zeta*zeta+89.0*xi*zeta*zeta-79.0*xi*zeta+23.0*xi)/(den*den*den*den) );
                      case 14:
                        return sign * RealGradient( -0.125*zeta*(eta-zeta+1.0)*(3.0*eta*zeta+7.0*eta+11.0*zeta*zeta-12.0*zeta+1.0)/(den*den*den), -0.25*zeta*(3.0*eta*xi*zeta+7.0*eta*xi-4.0*eta*zeta*zeta+8.0*eta*zeta-4.0*eta+2.0*xi*zeta*zeta-4.0*xi*zeta+2.0*xi+4.0*zeta*zeta*zeta-10.0*zeta*zeta+8.0*zeta-2.0)/(den*den*den), -0.25*(8.0*eta-11.0*xi-3.0*zeta)-0.25*(6.0*xi*eta*eta*zeta*zeta+16.0*xi*zeta*eta*eta-7.0*xi*eta*eta-5.0*eta*eta*zeta*zeta*zeta+14.0*eta*eta*zeta*zeta-13.0*zeta*eta*eta+4.0*eta*eta+10.0*eta*xi*zeta*zeta*zeta-28.0*eta*xi*zeta*zeta+26.0*eta*xi*zeta-8.0*eta*xi+6.0*eta*zeta*zeta*zeta-16.0*eta*zeta*zeta+14.0*eta*zeta-4.0*eta-15.0*xi*zeta*zeta*zeta+40.0*xi*zeta*zeta-35.0*xi*zeta+10.0*xi)/(den*den*den*den) );
                      case 15:
                        return sign * RealGradient( 1.25*zeta*(3.0*zeta-1.0)*(eta-zeta+1.0)*(eta+zeta-1.0)/(den*den*den), 0.5*zeta*(15.0*eta*xi*zeta-5.0*eta*xi+2.0*xi*zeta*zeta-4.0*xi*zeta+2.0*xi+4.0*zeta*zeta*zeta-10.0*zeta*zeta+8.0*zeta-2.0)/(den*den*den), 0.5*(8.0*eta-15.0*xi-3.0*zeta+1.0)+0.25*(60.0*xi*eta*eta*zeta*zeta-37.0*xi*zeta*eta*eta+7.0*xi*eta*eta-2.0*eta*eta*zeta*zeta*zeta+4.0*eta*eta*zeta*zeta-2.0*zeta*eta*eta+4.0*eta*xi*zeta*zeta*zeta-8.0*eta*xi*zeta*zeta+4.0*eta*xi*zeta+18.0*eta*zeta*zeta*zeta-50.0*eta*zeta*zeta+46.0*eta*zeta-14.0*eta-33.0*xi*zeta*zeta*zeta+89.0*xi*zeta*zeta-79.0*xi*zeta+23.0*xi)/(den*den*den*den) );
                      case 16:
                      case 17:
                      case 18:
                      case 19:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 4, i - 16);
                          return c[0] * RealGradient( 1.5*zeta*(3.0*zeta-1.0)*(eta-zeta+1.0)*(eta+zeta-1.0)/(den*den*den), -7.5*xi+1.5*xi*(24.0*eta*zeta*zeta-23.0*eta*zeta+3.0*eta+4.0*zeta*zeta-8.0*zeta+4.0)/(den*den*den), -9.0*xi+0.75*xi*(60.0*eta*eta*zeta*zeta-55.0*zeta*eta*eta+7.0*eta*eta+10.0*eta*zeta*zeta*zeta-22.0*eta*zeta*zeta+14.0*eta*zeta-2.0*eta-15.0*zeta*zeta*zeta+41.0*zeta*zeta-37.0*zeta+11.0)/(den*den*den*den) ) +
                                 c[1] * RealGradient( -1.5*zeta*(3.0*zeta-1.0)*(eta-zeta+1.0)*(eta+zeta-1.0)/(den*den*den), -7.5*xi-1.5*xi*(24.0*eta*zeta*zeta-23.0*eta*zeta+3.0*eta-4.0*zeta*zeta+8.0*zeta-4.0)/(den*den*den), 9.0*xi-0.75*xi*(60.0*eta*eta*zeta*zeta-55.0*zeta*eta*eta+7.0*eta*eta-10.0*eta*zeta*zeta*zeta+22.0*eta*zeta*zeta-14.0*eta*zeta+2.0*eta-15.0*zeta*zeta*zeta+41.0*zeta*zeta-37.0*zeta+11.0)/(den*den*den*den) ) +
                                 c[2] * RealGradient( 0.75*(eta-zeta+1.0)*(eta+zeta-1.0)*(24.0*zeta*zeta-23.0*zeta+3.0)/(den*den*den), 3.0*eta*xi*zeta*(3.0*zeta-1.0)/(den*den*den), -0.75*(48.0*xi+5.0*zeta-1.0)+0.75*(60.0*xi*eta*eta*zeta*zeta-55.0*xi*zeta*eta*eta+7.0*xi*eta*eta+5.0*eta*eta*zeta*zeta*zeta-11.0*eta*eta*zeta*zeta+7.0*zeta*eta*eta-eta*eta-45.0*xi*zeta*zeta*zeta+131.0*xi*zeta*zeta-127.0*xi*zeta+41.0*xi)/(den*den*den*den) ) +
                                 c[3] * RealGradient( -0.75*(eta-zeta+1.0)*(eta+zeta-1.0)*(24.0*zeta*zeta-23.0*zeta+3.0)/(den*den*den), -3.0*eta*xi*zeta*(3.0*zeta-1.0)/(den*den*den), 0.75*(48.0*xi-5.0*zeta+1.0)-0.75*(60.0*xi*eta*eta*zeta*zeta-55.0*xi*zeta*eta*eta+7.0*xi*eta*eta-5.0*eta*eta*zeta*zeta*zeta+11.0*eta*eta*zeta*zeta-7.0*zeta*eta*eta+eta*eta-45.0*xi*zeta*zeta*zeta+131.0*xi*zeta*zeta-127.0*xi*zeta+41.0*xi)/(den*den*den*den) );
                        }
                      case 20:
                      case 21:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 0, i - 20);
                          return c[0] * RealGradient( -0.5*zeta*(3.0*eta-zeta+1.0)*(eta+zeta-1.0)/(den*den), 0.0, 0.5*(2.0*xi+3.0*zeta)*(3.0*eta-zeta+1.0)*(eta+zeta-1.0)/(den*den) ) +
                                 c[1] * RealGradient( zeta*(3.0*eta-zeta+1.0)*(eta+zeta-1.0)/(den*den), 0.0, 2.0*xi-2.0*eta*xi*(3.0*eta+2.0*zeta-2.0)/(den*den) );
                        }
                      case 22:
                      case 23:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 1, i - 22);
                          return c[0] * RealGradient( 0.0, -zeta*(3.0*xi-zeta+1.0)*(eta+3.0*zeta-3.0)/(den*den), (3.0*xi-zeta+1.0)*(eta*eta+3.0*eta*zeta+zeta-1.0)/(den*den) ) +
                                 c[1] * RealGradient( 0.0, 2.0*eta*zeta*(3.0*xi-zeta+1.0)/(den*den), -2.0*(3.0*xi-zeta+1.0)*(eta*eta+zeta-1.0)/(den*den) );
                        }
                      case 24:
                      case 25:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 2, i - 24);
                          return c[0] * RealGradient( -0.5*zeta*(eta-zeta+1.0)*(3.0*eta+zeta-1.0)/(den*den), 0.0, 0.5*(2.0*xi-3.0*zeta)*(eta-zeta+1.0)*(3.0*eta+zeta-1.0)/(den*den) ) +
                                 c[1] * RealGradient( zeta*(eta-zeta+1.0)*(3.0*eta+zeta-1.0)/(den*den), 0.0, 2.0*xi-2.0*eta*xi*(3.0*eta-2.0*zeta+2.0)/(den*den) );
                        }
                      case 26:
                      case 27:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 3, i - 26);
                          return c[0] * RealGradient( 0.0, -zeta*(3.0*xi+zeta-1.0)*(eta-3.0*zeta+3.0)/(den*den), (3.0*xi+zeta-1.0)*(eta*eta-3.0*eta*zeta+zeta-1.0)/(den*den) ) +
                                 c[1] * RealGradient( 0.0, 2.0*eta*zeta*(3.0*xi+zeta-1.0)/(den*den), -2.0*(eta*eta+zeta-1.0)*(3.0*xi+zeta-1.0)/(den*den) );
                        }
                      case 28:
                        return RealGradient( 0.0, 0.0, 7.5*zeta-7.5*zeta*eta*eta/(den*den) );
                      case 29:
                        return RealGradient( 0.0, -15.0*xi*zeta/den, -15.0*eta*xi*zeta/(den*den) );
                      case 30:
                        return RealGradient( 3.75*zeta*(3.0*zeta-1.0)*(eta-zeta+1.0)*(eta+zeta-1.0)/(den*den*den), 7.5*eta*xi*zeta*(3.0*zeta-1.0)/(den*den*den), -22.5*xi+7.5*xi*(6.0*eta*eta*zeta*zeta-4.0*zeta*eta*eta+eta*eta-3.0*zeta*zeta*zeta+8.0*zeta*zeta-7.0*zeta+2.0)/(den*den*den*den) );
                      case 31:
                        return RealGradient( 7.5*zeta*(6.0*zeta-5.0)*(eta-zeta+1.0)*(eta+zeta-1.0)/(den*den*den), 7.5*eta*xi*zeta*(3.0*zeta-1.0)/(den*den*den), -90.0*xi+3.75*xi*(30.0*eta*eta*zeta*zeta-23.0*zeta*eta*eta-eta*eta-27.0*zeta*zeta*zeta+79.0*zeta*zeta-77.0*zeta+25.0)/(den*den*den*den) );
                      case 32:
                        return RealGradient( 3.75*zeta*(3.0*zeta-1.0)*(eta-zeta+1.0)*(eta+zeta-1.0)/(den*den*den), 15.0*eta*xi*zeta*(6.0*zeta-5.0)/(den*den*den), -22.5*xi+3.75*xi*(30.0*eta*eta*zeta*zeta-23.0*zeta*eta*eta-eta*eta-9.0*zeta*zeta*zeta+25.0*zeta*zeta-23.0*zeta+7.0)/(den*den*den*den) );
                      case 33:
                        return RealGradient( -18.75*zeta*(3.0*zeta-1.0)*(eta-zeta+1.0)*(eta+zeta-1.0)/(den*den*den), -37.5*eta*xi*zeta*(3.0*zeta-1.0)/(den*den*den), 112.5*xi-7.5*xi*(30.0*eta*eta*zeta*zeta-17.0*zeta*eta*eta+2.0*eta*eta-18.0*zeta*zeta*zeta+49.0*zeta*zeta-44.0*zeta+13.0)/(den*den*den*den) );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      } // switch(ii)
                  } // j = 0

                  // d()/deta
                case 1:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient( 0.25*(3.0*eta+zeta-1.0)*(4.0*xi*zeta-3.0*xi+4.0*zeta*zeta-5.0*zeta+1.0)/(den*den), 0.0, -0.25*(3.0*eta+zeta-1.0)*(4.0*xi*xi+4.0*xi*zeta-xi+zeta-1.0)/(den*den) );
                      case 1:
                        return sign * RealGradient( -0.25*(3.0*eta+zeta-1.0)*(4.0*xi*zeta-3.0*xi-4.0*zeta*zeta+5.0*zeta-1.0)/(den*den), 0.0, 0.25*(3.0*eta+zeta-1.0)*(4.0*xi*xi-4.0*xi*zeta+xi+zeta-1.0)/(den*den) );
                      case 2:
                        return sign * RealGradient( 0.0, 0.125*(4.0*zeta-3.0)*(xi-zeta+1.0)*(3.0*xi+zeta-1.0)/(den*den), -0.125*(xi-zeta+1.0)*(3.0*xi+zeta-1.0)*(8.0*eta+4.0*zeta-1.0)/(den*den) );
                      case 3:
                        return sign * RealGradient( 0.0, -0.125*(4.0*zeta-3.0)*(xi-zeta+1.0)*(3.0*xi+zeta-1.0)/(den*den), 0.125*(xi-zeta+1.0)*(8.0*eta-4.0*zeta+1.0)*(3.0*xi+zeta-1.0)/(den*den) );
                      case 4:
                        return sign * RealGradient( 0.25*(3.0*eta-zeta+1.0)*(4.0*xi*zeta-3.0*xi-4.0*zeta*zeta+5.0*zeta-1.0)/(den*den), 0.0, -0.25*(3.0*eta-zeta+1.0)*(4.0*xi*xi-4.0*xi*zeta+xi+zeta-1.0)/(den*den) );
                      case 5:
                        return sign * RealGradient( -0.25*(3.0*eta-zeta+1.0)*(4.0*xi*zeta-3.0*xi+4.0*zeta*zeta-5.0*zeta+1.0)/(den*den), 0.0, 0.25*(3.0*eta-zeta+1.0)*(4.0*xi*xi+4.0*xi*zeta-xi+zeta-1.0)/(den*den) );
                      case 6:
                        return sign * RealGradient( 0.0, 0.125*(4.0*zeta-3.0)*(3.0*xi-zeta+1.0)*(xi+zeta-1.0)/(den*den), -0.125*(3.0*xi-zeta+1.0)*(xi+zeta-1.0)*(8.0*eta+4.0*zeta-1.0)/(den*den) );
                      case 7:
                        return sign * RealGradient( 0.0, -0.125*(4.0*zeta-3.0)*(3.0*xi-zeta+1.0)*(xi+zeta-1.0)/(den*den), 0.125*(3.0*xi-zeta+1.0)*(8.0*eta-4.0*zeta+1.0)*(xi+zeta-1.0)/(den*den) );
                      case 8:
                        return sign * RealGradient( -0.25*zeta*(3.0*eta*xi*zeta+7.0*eta*xi-2.0*eta*zeta*zeta+4.0*eta*zeta-2.0*eta-4.0*xi*zeta*zeta+8.0*xi*zeta-4.0*xi-4.0*zeta*zeta*zeta+10.0*zeta*zeta-8.0*zeta+2.0)/(den*den*den), -0.125*zeta*(xi+zeta-1.0)*(3.0*xi*zeta+7.0*xi-11.0*zeta*zeta+12.0*zeta-1.0)/(den*den*den), 0.25*(11.0*eta+8.0*xi+3.0*zeta)-0.25*(6.0*eta*xi*xi*zeta*zeta+16.0*eta*zeta*xi*xi-7.0*eta*xi*xi-10.0*eta*xi*zeta*zeta*zeta+28.0*eta*xi*zeta*zeta-26.0*eta*xi*zeta+8.0*eta*xi-15.0*eta*zeta*zeta*zeta+40.0*eta*zeta*zeta-35.0*eta*zeta+10.0*eta-5.0*xi*xi*zeta*zeta*zeta+14.0*xi*xi*zeta*zeta-13.0*zeta*xi*xi+4.0*xi*xi-6.0*xi*zeta*zeta*zeta+16.0*xi*zeta*zeta-14.0*xi*zeta+4.0*xi)/(den*den*den*den) );
                      case 9:
                        return sign * RealGradient( 0.5*zeta*(15.0*eta*xi*zeta-5.0*eta*xi-2.0*eta*zeta*zeta+4.0*eta*zeta-2.0*eta-4.0*zeta*zeta*zeta+10.0*zeta*zeta-8.0*zeta+2.0)/(den*den*den), 1.25*zeta*(3.0*zeta-1.0)*(xi-zeta+1.0)*(xi+zeta-1.0)/(den*den*den), -0.5*(15.0*eta+8.0*xi+3.0*zeta-1.0)+0.25*(60.0*eta*xi*xi*zeta*zeta-37.0*eta*zeta*xi*xi+7.0*eta*xi*xi-4.0*eta*xi*zeta*zeta*zeta+8.0*eta*xi*zeta*zeta-4.0*eta*xi*zeta-33.0*eta*zeta*zeta*zeta+89.0*eta*zeta*zeta-79.0*eta*zeta+23.0*eta-2.0*xi*xi*zeta*zeta*zeta+4.0*xi*xi*zeta*zeta-2.0*zeta*xi*xi-18.0*xi*zeta*zeta*zeta+50.0*xi*zeta*zeta-46.0*xi*zeta+14.0*xi)/(den*den*den*den) );
                      case 10:
                        return sign * RealGradient( -0.25*zeta*(3.0*eta*xi*zeta+7.0*eta*xi+2.0*eta*zeta*zeta-4.0*eta*zeta+2.0*eta-4.0*xi*zeta*zeta+8.0*xi*zeta-4.0*xi+4.0*zeta*zeta*zeta-10.0*zeta*zeta+8.0*zeta-2.0)/(den*den*den), -0.125*zeta*(xi-zeta+1.0)*(3.0*xi*zeta+7.0*xi+11.0*zeta*zeta-12.0*zeta+1.0)/(den*den*den), 0.25*(11.0*eta-8.0*xi+3.0*zeta)-0.25*(6.0*eta*xi*xi*zeta*zeta+16.0*eta*zeta*xi*xi-7.0*eta*xi*xi+10.0*eta*xi*zeta*zeta*zeta-28.0*eta*xi*zeta*zeta+26.0*eta*xi*zeta-8.0*eta*xi-15.0*eta*zeta*zeta*zeta+40.0*eta*zeta*zeta-35.0*eta*zeta+10.0*eta-5.0*xi*xi*zeta*zeta*zeta+14.0*xi*xi*zeta*zeta-13.0*zeta*xi*xi+4.0*xi*xi+6.0*xi*zeta*zeta*zeta-16.0*xi*zeta*zeta+14.0*xi*zeta-4.0*xi)/(den*den*den*den) );
                      case 11:
                        return sign * RealGradient( 0.5*zeta*(15.0*eta*xi*zeta-5.0*eta*xi+2.0*eta*zeta*zeta-4.0*eta*zeta+2.0*eta+4.0*zeta*zeta*zeta-10.0*zeta*zeta+8.0*zeta-2.0)/(den*den*den), 1.25*zeta*(3.0*zeta-1.0)*(xi-zeta+1.0)*(xi+zeta-1.0)/(den*den*den), -0.5*(15.0*eta-8.0*xi+3.0*zeta-1.0)+0.25*(60.0*eta*xi*xi*zeta*zeta-37.0*eta*zeta*xi*xi+7.0*eta*xi*xi+4.0*eta*xi*zeta*zeta*zeta-8.0*eta*xi*zeta*zeta+4.0*eta*xi*zeta-33.0*eta*zeta*zeta*zeta+89.0*eta*zeta*zeta-79.0*eta*zeta+23.0*eta-2.0*xi*xi*zeta*zeta*zeta+4.0*xi*xi*zeta*zeta-2.0*zeta*xi*xi+18.0*xi*zeta*zeta*zeta-50.0*xi*zeta*zeta+46.0*xi*zeta-14.0*xi)/(den*den*den*den) );
                      case 12:
                        return sign * RealGradient( -0.25*zeta*(3.0*eta*xi*zeta+7.0*eta*xi+2.0*eta*zeta*zeta-4.0*eta*zeta+2.0*eta+4.0*xi*zeta*zeta-8.0*xi*zeta+4.0*xi-4.0*zeta*zeta*zeta+10.0*zeta*zeta-8.0*zeta+2.0)/(den*den*den), -0.125*zeta*(xi-zeta+1.0)*(3.0*xi*zeta+7.0*xi+11.0*zeta*zeta-12.0*zeta+1.0)/(den*den*den), 0.25*(11.0*eta+8.0*xi-3.0*zeta)-0.25*(6.0*eta*xi*xi*zeta*zeta+16.0*eta*zeta*xi*xi-7.0*eta*xi*xi+10.0*eta*xi*zeta*zeta*zeta-28.0*eta*xi*zeta*zeta+26.0*eta*xi*zeta-8.0*eta*xi-15.0*eta*zeta*zeta*zeta+40.0*eta*zeta*zeta-35.0*eta*zeta+10.0*eta+5.0*xi*xi*zeta*zeta*zeta-14.0*xi*xi*zeta*zeta+13.0*zeta*xi*xi-4.0*xi*xi-6.0*xi*zeta*zeta*zeta+16.0*xi*zeta*zeta-14.0*xi*zeta+4.0*xi)/(den*den*den*den) );
                      case 13:
                        return sign * RealGradient( 0.5*zeta*(15.0*eta*xi*zeta-5.0*eta*xi+2.0*eta*zeta*zeta-4.0*eta*zeta+2.0*eta-4.0*zeta*zeta*zeta+10.0*zeta*zeta-8.0*zeta+2.0)/(den*den*den), 1.25*zeta*(3.0*zeta-1.0)*(xi-zeta+1.0)*(xi+zeta-1.0)/(den*den*den), -0.5*(15.0*eta+8.0*xi-3.0*zeta+1.0)+0.25*(60.0*eta*xi*xi*zeta*zeta-37.0*eta*zeta*xi*xi+7.0*eta*xi*xi+4.0*eta*xi*zeta*zeta*zeta-8.0*eta*xi*zeta*zeta+4.0*eta*xi*zeta-33.0*eta*zeta*zeta*zeta+89.0*eta*zeta*zeta-79.0*eta*zeta+23.0*eta+2.0*xi*xi*zeta*zeta*zeta-4.0*xi*xi*zeta*zeta+2.0*zeta*xi*xi-18.0*xi*zeta*zeta*zeta+50.0*xi*zeta*zeta-46.0*xi*zeta+14.0*xi)/(den*den*den*den) );
                      case 14:
                        return sign * RealGradient( -0.25*zeta*(3.0*eta*xi*zeta+7.0*eta*xi-2.0*eta*zeta*zeta+4.0*eta*zeta-2.0*eta+4.0*xi*zeta*zeta-8.0*xi*zeta+4.0*xi+4.0*zeta*zeta*zeta-10.0*zeta*zeta+8.0*zeta-2.0)/(den*den*den), -0.125*zeta*(xi+zeta-1.0)*(3.0*xi*zeta+7.0*xi-11.0*zeta*zeta+12.0*zeta-1.0)/(den*den*den), 0.25*(11.0*eta-8.0*xi-3.0*zeta)-0.25*(6.0*eta*xi*xi*zeta*zeta+16.0*eta*zeta*xi*xi-7.0*eta*xi*xi-10.0*eta*xi*zeta*zeta*zeta+28.0*eta*xi*zeta*zeta-26.0*eta*xi*zeta+8.0*eta*xi-15.0*eta*zeta*zeta*zeta+40.0*eta*zeta*zeta-35.0*eta*zeta+10.0*eta+5.0*xi*xi*zeta*zeta*zeta-14.0*xi*xi*zeta*zeta+13.0*zeta*xi*xi-4.0*xi*xi+6.0*xi*zeta*zeta*zeta-16.0*xi*zeta*zeta+14.0*xi*zeta-4.0*xi)/(den*den*den*den) );
                      case 15:
                        return sign * RealGradient( 0.5*zeta*(15.0*eta*xi*zeta-5.0*eta*xi-2.0*eta*zeta*zeta+4.0*eta*zeta-2.0*eta+4.0*zeta*zeta*zeta-10.0*zeta*zeta+8.0*zeta-2.0)/(den*den*den), 1.25*zeta*(3.0*zeta-1.0)*(xi-zeta+1.0)*(xi+zeta-1.0)/(den*den*den), -0.5*(15.0*eta-8.0*xi-3.0*zeta+1.0)+0.25*(60.0*eta*xi*xi*zeta*zeta-37.0*eta*zeta*xi*xi+7.0*eta*xi*xi-4.0*eta*xi*zeta*zeta*zeta+8.0*eta*xi*zeta*zeta-4.0*eta*xi*zeta-33.0*eta*zeta*zeta*zeta+89.0*eta*zeta*zeta-79.0*eta*zeta+23.0*eta+2.0*xi*xi*zeta*zeta*zeta-4.0*xi*xi*zeta*zeta+2.0*zeta*xi*xi+18.0*xi*zeta*zeta*zeta-50.0*xi*zeta*zeta+46.0*xi*zeta-14.0*xi)/(den*den*den*den) );
                      case 16:
                      case 17:
                      case 18:
                      case 19:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 4, i - 16);
                          return c[0] * RealGradient( 3.0*eta*xi*zeta*(3.0*zeta-1.0)/(den*den*den), 0.75*(xi-zeta+1.0)*(xi+zeta-1.0)*(24.0*zeta*zeta-23.0*zeta+3.0)/(den*den*den), -0.75*(48.0*eta+5.0*zeta-1.0)+0.75*(60.0*eta*xi*xi*zeta*zeta-55.0*eta*zeta*xi*xi+7.0*eta*xi*xi-45.0*eta*zeta*zeta*zeta+131.0*eta*zeta*zeta-127.0*eta*zeta+41.0*eta+5.0*xi*xi*zeta*zeta*zeta-11.0*xi*xi*zeta*zeta+7.0*zeta*xi*xi-xi*xi)/(den*den*den*den) ) +
                                 c[1] * RealGradient( -3.0*eta*xi*zeta*(3.0*zeta-1.0)/(den*den*den), -0.75*(xi-zeta+1.0)*(xi+zeta-1.0)*(24.0*zeta*zeta-23.0*zeta+3.0)/(den*den*den), 0.75*(48.0*eta-5.0*zeta+1.0)-0.75*(60.0*eta*xi*xi*zeta*zeta-55.0*eta*zeta*xi*xi+7.0*eta*xi*xi-45.0*eta*zeta*zeta*zeta+131.0*eta*zeta*zeta-127.0*eta*zeta+41.0*eta-5.0*xi*xi*zeta*zeta*zeta+11.0*xi*xi*zeta*zeta-7.0*zeta*xi*xi+xi*xi)/(den*den*den*den) ) +
                                 c[2] * RealGradient( -7.5*eta+1.5*eta*(24.0*xi*zeta*zeta-23.0*xi*zeta+3.0*xi+4.0*zeta*zeta-8.0*zeta+4.0)/(den*den*den), 1.5*zeta*(3.0*zeta-1.0)*(xi-zeta+1.0)*(xi+zeta-1.0)/(den*den*den), -9.0*eta+0.75*eta*(60.0*xi*xi*zeta*zeta-55.0*zeta*xi*xi+7.0*xi*xi+10.0*xi*zeta*zeta*zeta-22.0*xi*zeta*zeta+14.0*xi*zeta-2.0*xi-15.0*zeta*zeta*zeta+41.0*zeta*zeta-37.0*zeta+11.0)/(den*den*den*den) ) +
                                 c[3] * RealGradient( -7.5*eta-1.5*eta*(24.0*xi*zeta*zeta-23.0*xi*zeta+3.0*xi-4.0*zeta*zeta+8.0*zeta-4.0)/(den*den*den), -1.5*zeta*(3.0*zeta-1.0)*(xi-zeta+1.0)*(xi+zeta-1.0)/(den*den*den), 9.0*eta-0.75*eta*(60.0*xi*xi*zeta*zeta-55.0*zeta*xi*xi+7.0*xi*xi-10.0*xi*zeta*zeta*zeta+22.0*xi*zeta*zeta-14.0*xi*zeta+2.0*xi-15.0*zeta*zeta*zeta+41.0*zeta*zeta-37.0*zeta+11.0)/(den*den*den*den) );
                        }
                      case 20:
                      case 21:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 0, i - 20);
                          return c[0] * RealGradient( -zeta*(3.0*eta+zeta-1.0)*(xi+3.0*zeta-3.0)/(den*den), 0.0, (3.0*eta+zeta-1.0)*(xi*xi+3.0*xi*zeta+zeta-1.0)/(den*den) ) +
                                 c[1] * RealGradient( 2.0*xi*zeta*(3.0*eta+zeta-1.0)/(den*den), 0.0, -2.0*(xi*xi+zeta-1.0)*(3.0*eta+zeta-1.0)/(den*den) );
                        }
                      case 22:
                      case 23:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 1, i - 22);
                          return c[0] * RealGradient( 0.0, -0.5*zeta*(xi-zeta+1.0)*(3.0*xi+zeta-1.0)/(den*den), 0.5*(2.0*eta+3.0*zeta)*(xi-zeta+1.0)*(3.0*xi+zeta-1.0)/(den*den) ) +
                                 c[1] * RealGradient( 0.0, zeta*(xi-zeta+1.0)*(3.0*xi+zeta-1.0)/(den*den), -2.0*eta*(xi-zeta+1.0)*(3.0*xi+zeta-1.0)/(den*den) );
                        }
                      case 24:
                      case 25:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 2, i - 24);
                          return c[0] * RealGradient( -zeta*(3.0*eta-zeta+1.0)*(xi-3.0*zeta+3.0)/(den*den), 0.0, (3.0*eta-zeta+1.0)*(xi*xi-3.0*xi*zeta+zeta-1.0)/(den*den) ) +
                                 c[1] * RealGradient( 2.0*xi*zeta*(3.0*eta-zeta+1.0)/(den*den), 0.0, -2.0*(3.0*eta-zeta+1.0)*(xi*xi+zeta-1.0)/(den*den) );
                        }
                      case 26:
                      case 27:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 3, i - 26);
                          return c[0] * RealGradient( 0.0, -0.5*zeta*(3.0*xi-zeta+1.0)*(xi+zeta-1.0)/(den*den), 0.5*(2.0*eta-3.0*zeta)*(3.0*xi-zeta+1.0)*(xi+zeta-1.0)/(den*den) ) +
                                 c[1] * RealGradient( 0.0, zeta*(3.0*xi-zeta+1.0)*(xi+zeta-1.0)/(den*den), -2.0*eta*(3.0*xi-zeta+1.0)*(xi+zeta-1.0)/(den*den) );
                        }
                      case 28:
                        return RealGradient( -15.0*eta*zeta/den, 0.0, -15.0*eta*xi*zeta/(den*den) );
                      case 29:
                        return RealGradient( 0.0, 0.0, 7.5*zeta-7.5*zeta*xi*xi/(den*den) );
                      case 30:
                        return RealGradient( 7.5*eta*xi*zeta*(3.0*zeta-1.0)/(den*den*den), 3.75*zeta*(3.0*zeta-1.0)*(xi-zeta+1.0)*(xi+zeta-1.0)/(den*den*den), -22.5*eta+7.5*eta*(6.0*xi*xi*zeta*zeta-4.0*zeta*xi*xi+xi*xi-3.0*zeta*zeta*zeta+8.0*zeta*zeta-7.0*zeta+2.0)/(den*den*den*den) );
                      case 31:
                        return RealGradient( 15.0*eta*xi*zeta*(6.0*zeta-5.0)/(den*den*den), 3.75*zeta*(3.0*zeta-1.0)*(xi-zeta+1.0)*(xi+zeta-1.0)/(den*den*den), -22.5*eta+3.75*eta*(30.0*xi*xi*zeta*zeta-23.0*zeta*xi*xi-xi*xi-9.0*zeta*zeta*zeta+25.0*zeta*zeta-23.0*zeta+7.0)/(den*den*den*den) );
                      case 32:
                        return RealGradient( 7.5*eta*xi*zeta*(3.0*zeta-1.0)/(den*den*den), 7.5*zeta*(6.0*zeta-5.0)*(xi-zeta+1.0)*(xi+zeta-1.0)/(den*den*den), -90.0*eta+3.75*eta*(30.0*xi*xi*zeta*zeta-23.0*zeta*xi*xi-xi*xi-27.0*zeta*zeta*zeta+79.0*zeta*zeta-77.0*zeta+25.0)/(den*den*den*den) );
                      case 33:
                        return RealGradient( -37.5*eta*xi*zeta*(3.0*zeta-1.0)/(den*den*den), -18.75*zeta*(3.0*zeta-1.0)*(xi-zeta+1.0)*(xi+zeta-1.0)/(den*den*den), 112.5*eta-7.5*eta*(30.0*xi*xi*zeta*zeta-17.0*zeta*xi*xi+2.0*xi*xi-18.0*zeta*zeta*zeta+49.0*zeta*zeta-44.0*zeta+13.0)/(den*den*den*den) );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      } // switch(ii)
                  } // j = 1

                  // d()/dzeta
                case 2:
                  {
                    switch(ii)
                      {
                      case 0:
                        return sign * RealGradient( 0.125*(8.0*eta-4.0*xi-8.0*zeta+5.0)+0.125*eta*(12.0*eta*xi*zeta-6.0*eta*xi+9.0*eta*zeta-9.0*eta+2.0*xi*zeta-2.0*xi)/(den*den*den), 0.0, 0.125*(4.0*xi+1.0)-0.125*eta*(24.0*eta*xi*xi+12.0*eta*xi*zeta+6.0*eta*xi+3.0*eta*zeta-3.0*eta+8.0*zeta*xi*xi-8.0*xi*xi+6.0*xi*zeta-6.0*xi)/(den*den*den) );
                      case 1:
                        return sign * RealGradient( 0.125*(8.0*eta+4.0*xi-8.0*zeta+5.0)-0.125*eta*(12.0*eta*xi*zeta-6.0*eta*xi-9.0*eta*zeta+9.0*eta+2.0*xi*zeta-2.0*xi)/(den*den*den), 0.0, 0.125*(4.0*xi-1.0)+0.125*eta*(24.0*eta*xi*xi-12.0*eta*xi*zeta-6.0*eta*xi+3.0*eta*zeta-3.0*eta+8.0*zeta*xi*xi-8.0*xi*xi-6.0*xi*zeta+6.0*xi)/(den*den*den) );
                      case 2:
                        return sign * RealGradient( 0.0, -0.125*(4.0*eta+8.0*xi+8.0*zeta-5.0)+0.125*xi*(12.0*eta*xi*zeta-6.0*eta*xi-2.0*eta*zeta+2.0*eta+9.0*xi*zeta-9.0*xi)/(den*den*den), 0.125*(4.0*eta+1.0)-0.125*xi*(24.0*xi*eta*eta-8.0*zeta*eta*eta+8.0*eta*eta+12.0*eta*xi*zeta+6.0*eta*xi-6.0*eta*zeta+6.0*eta+3.0*xi*zeta-3.0*xi)/(den*den*den) );
                      case 3:
                        return sign * RealGradient( 0.0, 0.125*(4.0*eta-8.0*xi-8.0*zeta+5.0)-0.125*xi*(12.0*eta*xi*zeta-6.0*eta*xi-2.0*eta*zeta+2.0*eta-9.0*xi*zeta+9.0*xi)/(den*den*den), 0.125*(4.0*eta-1.0)+0.125*xi*(24.0*xi*eta*eta-8.0*zeta*eta*eta+8.0*eta*eta-12.0*eta*xi*zeta-6.0*eta*xi+6.0*eta*zeta-6.0*eta+3.0*xi*zeta-3.0*xi)/(den*den*den) );
                      case 4:
                        return sign * RealGradient( 0.125*(8.0*eta-4.0*xi+8.0*zeta-5.0)+0.125*eta*(12.0*eta*xi*zeta-6.0*eta*xi-9.0*eta*zeta+9.0*eta-2.0*xi*zeta+2.0*xi)/(den*den*den), 0.0, -0.125*(4.0*xi-1.0)-0.125*eta*(24.0*eta*xi*xi-12.0*eta*xi*zeta-6.0*eta*xi+3.0*eta*zeta-3.0*eta-8.0*zeta*xi*xi+8.0*xi*xi+6.0*xi*zeta-6.0*xi)/(den*den*den) );
                      case 5:
                        return sign * RealGradient( 0.125*(8.0*eta+4.0*xi+8.0*zeta-5.0)-0.125*eta*(12.0*eta*xi*zeta-6.0*eta*xi+9.0*eta*zeta-9.0*eta-2.0*xi*zeta+2.0*xi)/(den*den*den), 0.0, -0.125*(4.0*xi+1.0)+0.125*eta*(24.0*eta*xi*xi+12.0*eta*xi*zeta+6.0*eta*xi+3.0*eta*zeta-3.0*eta-8.0*zeta*xi*xi+8.0*xi*xi-6.0*xi*zeta+6.0*xi)/(den*den*den) );
                      case 6:
                        return sign * RealGradient( 0.0, -0.125*(4.0*eta-8.0*xi+8.0*zeta-5.0)+0.125*xi*(12.0*eta*xi*zeta-6.0*eta*xi+2.0*eta*zeta-2.0*eta+9.0*xi*zeta-9.0*xi)/(den*den*den), 0.125*(4.0*eta+1.0)-0.125*xi*(24.0*xi*eta*eta+8.0*zeta*eta*eta-8.0*eta*eta+12.0*eta*xi*zeta+6.0*eta*xi+6.0*eta*zeta-6.0*eta+3.0*xi*zeta-3.0*xi)/(den*den*den) );
                      case 7:
                        return sign * RealGradient( 0.0, 0.125*(4.0*eta+8.0*xi-8.0*zeta+5.0)-0.125*xi*(12.0*eta*xi*zeta-6.0*eta*xi+2.0*eta*zeta-2.0*eta-9.0*xi*zeta+9.0*xi)/(den*den*den), 0.125*(4.0*eta-1.0)+0.125*xi*(24.0*xi*eta*eta+8.0*zeta*eta*eta-8.0*eta*eta-12.0*eta*xi*zeta-6.0*eta*xi-6.0*eta*zeta+6.0*eta+3.0*xi*zeta-3.0*xi)/(den*den*den) );
                      case 8:
                        return sign * RealGradient( -0.125*(8.0*eta+11.0*xi+12.0*zeta-2.0)-0.125*(3.0*xi*eta*eta*zeta*zeta+20.0*xi*zeta*eta*eta+7.0*xi*eta*eta-2.0*eta*eta*zeta*zeta+4.0*zeta*eta*eta-2.0*eta*eta-8.0*eta*xi*zeta*zeta+16.0*eta*xi*zeta-8.0*eta*xi-4.0*eta*zeta*zeta+8.0*eta*zeta-4.0*eta-10.0*xi*zeta*zeta+20.0*xi*zeta-10.0*xi)/(den*den*den*den), -0.125*(11.0*eta+8.0*xi+12.0*zeta-2.0)-0.125*(3.0*eta*xi*xi*zeta*zeta+20.0*eta*zeta*xi*xi+7.0*eta*xi*xi-8.0*eta*xi*zeta*zeta+16.0*eta*xi*zeta-8.0*eta*xi-10.0*eta*zeta*zeta+20.0*eta*zeta-10.0*eta-2.0*xi*xi*zeta*zeta+4.0*zeta*xi*xi-2.0*xi*xi-4.0*xi*zeta*zeta+8.0*xi*zeta-4.0*xi)/(den*den*den*den), 0.25*(3.0*eta+3.0*xi+1.0)-0.125*(12.0*eta*eta*xi*xi*zeta*zeta+60.0*zeta*eta*eta*xi*xi-12.0*eta*eta*xi*xi-10.0*xi*eta*eta*zeta*zeta*zeta+26.0*xi*eta*eta*zeta*zeta-22.0*xi*zeta*eta*eta+6.0*xi*eta*eta-15.0*eta*eta*zeta*zeta*zeta+35.0*eta*eta*zeta*zeta-25.0*zeta*eta*eta+5.0*eta*eta-10.0*eta*xi*xi*zeta*zeta*zeta+26.0*eta*xi*xi*zeta*zeta-22.0*eta*zeta*xi*xi+6.0*eta*xi*xi-12.0*eta*xi*zeta*zeta*zeta+28.0*eta*xi*zeta*zeta-20.0*eta*xi*zeta+4.0*eta*xi-15.0*xi*xi*zeta*zeta*zeta+35.0*xi*xi*zeta*zeta-25.0*zeta*xi*xi+5.0*xi*xi)/(den*den*den*den*den) );
                      case 9:
                        return sign * RealGradient( 0.25*(8.0*eta+15.0*xi+12.0*zeta-2.0)+0.25*(15.0*xi*eta*eta*zeta*zeta+20.0*xi*zeta*eta*eta-5.0*xi*eta*eta-2.0*eta*eta*zeta*zeta+4.0*zeta*eta*eta-2.0*eta*eta-4.0*eta*zeta*zeta+8.0*eta*zeta-4.0*eta-10.0*xi*zeta*zeta+20.0*xi*zeta-10.0*xi)/(den*den*den*den), 0.25*(15.0*eta+8.0*xi+12.0*zeta-2.0)+0.25*(15.0*eta*xi*xi*zeta*zeta+20.0*eta*zeta*xi*xi-5.0*eta*xi*xi-10.0*eta*zeta*zeta+20.0*eta*zeta-10.0*eta-2.0*xi*xi*zeta*zeta+4.0*zeta*xi*xi-2.0*xi*xi-4.0*xi*zeta*zeta+8.0*xi*zeta-4.0*xi)/(den*den*den*den), -0.125*(12.0*eta+12.0*xi+5.0)+0.125*(120.0*eta*eta*xi*xi*zeta*zeta+9.0*zeta*eta*eta*xi*xi-9.0*eta*eta*xi*xi-4.0*xi*eta*eta*zeta*zeta*zeta+4.0*xi*eta*eta*zeta*zeta+4.0*xi*zeta*eta*eta-4.0*xi*eta*eta-33.0*eta*eta*zeta*zeta*zeta+79.0*eta*eta*zeta*zeta-59.0*zeta*eta*eta+13.0*eta*eta-4.0*eta*xi*xi*zeta*zeta*zeta+4.0*eta*xi*xi*zeta*zeta+4.0*eta*zeta*xi*xi-4.0*eta*xi*xi-36.0*eta*xi*zeta*zeta*zeta+92.0*eta*xi*zeta*zeta-76.0*eta*xi*zeta+20.0*eta*xi-33.0*xi*xi*zeta*zeta*zeta+79.0*xi*xi*zeta*zeta-59.0*zeta*xi*xi+13.0*xi*xi)/(den*den*den*den*den) );
                      case 10:
                        return sign * RealGradient( 0.125*(8.0*eta-11.0*xi+12.0*zeta-2.0)-0.125*(3.0*xi*eta*eta*zeta*zeta+20.0*xi*zeta*eta*eta+7.0*xi*eta*eta+2.0*eta*eta*zeta*zeta-4.0*zeta*eta*eta+2.0*eta*eta-8.0*eta*xi*zeta*zeta+16.0*eta*xi*zeta-8.0*eta*xi+4.0*eta*zeta*zeta-8.0*eta*zeta+4.0*eta-10.0*xi*zeta*zeta+20.0*xi*zeta-10.0*xi)/(den*den*den*den), -0.125*(11.0*eta-8.0*xi+12.0*zeta-2.0)-0.125*(3.0*eta*xi*xi*zeta*zeta+20.0*eta*zeta*xi*xi+7.0*eta*xi*xi+8.0*eta*xi*zeta*zeta-16.0*eta*xi*zeta+8.0*eta*xi-10.0*eta*zeta*zeta+20.0*eta*zeta-10.0*eta-2.0*xi*xi*zeta*zeta+4.0*zeta*xi*xi-2.0*xi*xi+4.0*xi*zeta*zeta-8.0*xi*zeta+4.0*xi)/(den*den*den*den), 0.25*(3.0*eta-3.0*xi+1.0)-0.125*(12.0*eta*eta*xi*xi*zeta*zeta+60.0*zeta*eta*eta*xi*xi-12.0*eta*eta*xi*xi+10.0*xi*eta*eta*zeta*zeta*zeta-26.0*xi*eta*eta*zeta*zeta+22.0*xi*zeta*eta*eta-6.0*xi*eta*eta-15.0*eta*eta*zeta*zeta*zeta+35.0*eta*eta*zeta*zeta-25.0*zeta*eta*eta+5.0*eta*eta-10.0*eta*xi*xi*zeta*zeta*zeta+26.0*eta*xi*xi*zeta*zeta-22.0*eta*zeta*xi*xi+6.0*eta*xi*xi+12.0*eta*xi*zeta*zeta*zeta-28.0*eta*xi*zeta*zeta+20.0*eta*xi*zeta-4.0*eta*xi-15.0*xi*xi*zeta*zeta*zeta+35.0*xi*xi*zeta*zeta-25.0*zeta*xi*xi+5.0*xi*xi)/(den*den*den*den*den) );
                      case 11:
                        return sign * RealGradient( -0.25*(8.0*eta-15.0*xi+12.0*zeta-2.0)+0.25*(15.0*xi*eta*eta*zeta*zeta+20.0*xi*zeta*eta*eta-5.0*xi*eta*eta+2.0*eta*eta*zeta*zeta-4.0*zeta*eta*eta+2.0*eta*eta+4.0*eta*zeta*zeta-8.0*eta*zeta+4.0*eta-10.0*xi*zeta*zeta+20.0*xi*zeta-10.0*xi)/(den*den*den*den), 0.25*(15.0*eta-8.0*xi+12.0*zeta-2.0)+0.25*(15.0*eta*xi*xi*zeta*zeta+20.0*eta*zeta*xi*xi-5.0*eta*xi*xi-10.0*eta*zeta*zeta+20.0*eta*zeta-10.0*eta-2.0*xi*xi*zeta*zeta+4.0*zeta*xi*xi-2.0*xi*xi+4.0*xi*zeta*zeta-8.0*xi*zeta+4.0*xi)/(den*den*den*den), -0.125*(12.0*eta-12.0*xi+5.0)+0.125*(120.0*eta*eta*xi*xi*zeta*zeta+9.0*zeta*eta*eta*xi*xi-9.0*eta*eta*xi*xi+4.0*xi*eta*eta*zeta*zeta*zeta-4.0*xi*eta*eta*zeta*zeta-4.0*xi*zeta*eta*eta+4.0*xi*eta*eta-33.0*eta*eta*zeta*zeta*zeta+79.0*eta*eta*zeta*zeta-59.0*zeta*eta*eta+13.0*eta*eta-4.0*eta*xi*xi*zeta*zeta*zeta+4.0*eta*xi*xi*zeta*zeta+4.0*eta*zeta*xi*xi-4.0*eta*xi*xi+36.0*eta*xi*zeta*zeta*zeta-92.0*eta*xi*zeta*zeta+76.0*eta*xi*zeta-20.0*eta*xi-33.0*xi*xi*zeta*zeta*zeta+79.0*xi*xi*zeta*zeta-59.0*zeta*xi*xi+13.0*xi*xi)/(den*den*den*den*den) );
                      case 12:
                        return sign * RealGradient( -0.125*(8.0*eta+11.0*xi-12.0*zeta+2.0)-0.125*(3.0*xi*eta*eta*zeta*zeta+20.0*xi*zeta*eta*eta+7.0*xi*eta*eta+2.0*eta*eta*zeta*zeta-4.0*zeta*eta*eta+2.0*eta*eta+8.0*eta*xi*zeta*zeta-16.0*eta*xi*zeta+8.0*eta*xi-4.0*eta*zeta*zeta+8.0*eta*zeta-4.0*eta-10.0*xi*zeta*zeta+20.0*xi*zeta-10.0*xi)/(den*den*den*den), -0.125*(11.0*eta+8.0*xi-12.0*zeta+2.0)-0.125*(3.0*eta*xi*xi*zeta*zeta+20.0*eta*zeta*xi*xi+7.0*eta*xi*xi+8.0*eta*xi*zeta*zeta-16.0*eta*xi*zeta+8.0*eta*xi-10.0*eta*zeta*zeta+20.0*eta*zeta-10.0*eta+2.0*xi*xi*zeta*zeta-4.0*zeta*xi*xi+2.0*xi*xi-4.0*xi*zeta*zeta+8.0*xi*zeta-4.0*xi)/(den*den*den*den), -0.25*(3.0*eta+3.0*xi-1.0)-0.125*(12.0*eta*eta*xi*xi*zeta*zeta+60.0*zeta*eta*eta*xi*xi-12.0*eta*eta*xi*xi+10.0*xi*eta*eta*zeta*zeta*zeta-26.0*xi*eta*eta*zeta*zeta+22.0*xi*zeta*eta*eta-6.0*xi*eta*eta-15.0*eta*eta*zeta*zeta*zeta+35.0*eta*eta*zeta*zeta-25.0*zeta*eta*eta+5.0*eta*eta+10.0*eta*xi*xi*zeta*zeta*zeta-26.0*eta*xi*xi*zeta*zeta+22.0*eta*zeta*xi*xi-6.0*eta*xi*xi-12.0*eta*xi*zeta*zeta*zeta+28.0*eta*xi*zeta*zeta-20.0*eta*xi*zeta+4.0*eta*xi-15.0*xi*xi*zeta*zeta*zeta+35.0*xi*xi*zeta*zeta-25.0*zeta*xi*xi+5.0*xi*xi)/(den*den*den*den*den) );
                      case 13:
                        return sign * RealGradient( 0.25*(8.0*eta+15.0*xi-12.0*zeta+2.0)+0.25*(15.0*xi*eta*eta*zeta*zeta+20.0*xi*zeta*eta*eta-5.0*xi*eta*eta+2.0*eta*eta*zeta*zeta-4.0*zeta*eta*eta+2.0*eta*eta-4.0*eta*zeta*zeta+8.0*eta*zeta-4.0*eta-10.0*xi*zeta*zeta+20.0*xi*zeta-10.0*xi)/(den*den*den*den), 0.25*(15.0*eta+8.0*xi-12.0*zeta+2.0)+0.25*(15.0*eta*xi*xi*zeta*zeta+20.0*eta*zeta*xi*xi-5.0*eta*xi*xi-10.0*eta*zeta*zeta+20.0*eta*zeta-10.0*eta+2.0*xi*xi*zeta*zeta-4.0*zeta*xi*xi+2.0*xi*xi-4.0*xi*zeta*zeta+8.0*xi*zeta-4.0*xi)/(den*den*den*den), 0.125*(12.0*eta+12.0*xi-5.0)+0.125*(120.0*eta*eta*xi*xi*zeta*zeta+9.0*zeta*eta*eta*xi*xi-9.0*eta*eta*xi*xi+4.0*xi*eta*eta*zeta*zeta*zeta-4.0*xi*eta*eta*zeta*zeta-4.0*xi*zeta*eta*eta+4.0*xi*eta*eta-33.0*eta*eta*zeta*zeta*zeta+79.0*eta*eta*zeta*zeta-59.0*zeta*eta*eta+13.0*eta*eta+4.0*eta*xi*xi*zeta*zeta*zeta-4.0*eta*xi*xi*zeta*zeta-4.0*eta*zeta*xi*xi+4.0*eta*xi*xi-36.0*eta*xi*zeta*zeta*zeta+92.0*eta*xi*zeta*zeta-76.0*eta*xi*zeta+20.0*eta*xi-33.0*xi*xi*zeta*zeta*zeta+79.0*xi*xi*zeta*zeta-59.0*zeta*xi*xi+13.0*xi*xi)/(den*den*den*den*den) );
                      case 14:
                        return sign * RealGradient( 0.125*(8.0*eta-11.0*xi-12.0*zeta+2.0)-0.125*(3.0*xi*eta*eta*zeta*zeta+20.0*xi*zeta*eta*eta+7.0*xi*eta*eta-2.0*eta*eta*zeta*zeta+4.0*zeta*eta*eta-2.0*eta*eta+8.0*eta*xi*zeta*zeta-16.0*eta*xi*zeta+8.0*eta*xi+4.0*eta*zeta*zeta-8.0*eta*zeta+4.0*eta-10.0*xi*zeta*zeta+20.0*xi*zeta-10.0*xi)/(den*den*den*den), -0.125*(11.0*eta-8.0*xi-12.0*zeta+2.0)-0.125*(3.0*eta*xi*xi*zeta*zeta+20.0*eta*zeta*xi*xi+7.0*eta*xi*xi-8.0*eta*xi*zeta*zeta+16.0*eta*xi*zeta-8.0*eta*xi-10.0*eta*zeta*zeta+20.0*eta*zeta-10.0*eta+2.0*xi*xi*zeta*zeta-4.0*zeta*xi*xi+2.0*xi*xi+4.0*xi*zeta*zeta-8.0*xi*zeta+4.0*xi)/(den*den*den*den), -0.25*(3.0*eta-3.0*xi-1.0)-0.125*(12.0*eta*eta*xi*xi*zeta*zeta+60.0*zeta*eta*eta*xi*xi-12.0*eta*eta*xi*xi-10.0*xi*eta*eta*zeta*zeta*zeta+26.0*xi*eta*eta*zeta*zeta-22.0*xi*zeta*eta*eta+6.0*xi*eta*eta-15.0*eta*eta*zeta*zeta*zeta+35.0*eta*eta*zeta*zeta-25.0*zeta*eta*eta+5.0*eta*eta+10.0*eta*xi*xi*zeta*zeta*zeta-26.0*eta*xi*xi*zeta*zeta+22.0*eta*zeta*xi*xi-6.0*eta*xi*xi+12.0*eta*xi*zeta*zeta*zeta-28.0*eta*xi*zeta*zeta+20.0*eta*xi*zeta-4.0*eta*xi-15.0*xi*xi*zeta*zeta*zeta+35.0*xi*xi*zeta*zeta-25.0*zeta*xi*xi+5.0*xi*xi)/(den*den*den*den*den) );
                      case 15:
                        return sign * RealGradient( -0.25*(8.0*eta-15.0*xi-12.0*zeta+2.0)+0.25*(15.0*xi*eta*eta*zeta*zeta+20.0*xi*zeta*eta*eta-5.0*xi*eta*eta-2.0*eta*eta*zeta*zeta+4.0*zeta*eta*eta-2.0*eta*eta+4.0*eta*zeta*zeta-8.0*eta*zeta+4.0*eta-10.0*xi*zeta*zeta+20.0*xi*zeta-10.0*xi)/(den*den*den*den), 0.25*(15.0*eta-8.0*xi-12.0*zeta+2.0)+0.25*(15.0*eta*xi*xi*zeta*zeta+20.0*eta*zeta*xi*xi-5.0*eta*xi*xi-10.0*eta*zeta*zeta+20.0*eta*zeta-10.0*eta+2.0*xi*xi*zeta*zeta-4.0*zeta*xi*xi+2.0*xi*xi+4.0*xi*zeta*zeta-8.0*xi*zeta+4.0*xi)/(den*den*den*den), 0.125*(12.0*eta-12.0*xi-5.0)+0.125*(120.0*eta*eta*xi*xi*zeta*zeta+9.0*zeta*eta*eta*xi*xi-9.0*eta*eta*xi*xi-4.0*xi*eta*eta*zeta*zeta*zeta+4.0*xi*eta*eta*zeta*zeta+4.0*xi*zeta*eta*eta-4.0*xi*eta*eta-33.0*eta*eta*zeta*zeta*zeta+79.0*eta*eta*zeta*zeta-59.0*zeta*eta*eta+13.0*eta*eta+4.0*eta*xi*xi*zeta*zeta*zeta-4.0*eta*xi*xi*zeta*zeta-4.0*eta*zeta*xi*xi+4.0*eta*xi*xi+36.0*eta*xi*zeta*zeta*zeta-92.0*eta*xi*zeta*zeta+76.0*eta*xi*zeta-20.0*eta*xi-33.0*xi*xi*zeta*zeta*zeta+79.0*xi*xi*zeta*zeta-59.0*zeta*xi*xi+13.0*xi*xi)/(den*den*den*den*den) );
                      case 16:
                      case 17:
                      case 18:
                      case 19:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 4, i - 16);
                          return c[0] * RealGradient( 4.5*xi+1.5*xi*(3.0*eta*eta*zeta*zeta+4.0*zeta*eta*eta-eta*eta-2.0*zeta*zeta+4.0*zeta-2.0)/(den*den*den*den), 1.5*(12.0*eta+5.0*zeta-3.0)+1.5*(12.0*eta*xi*xi*zeta*zeta+eta*zeta*xi*xi-7.0*eta*xi*xi-2.0*eta*zeta*zeta+4.0*eta*zeta-2.0*eta+2.0*xi*xi*zeta*zeta-4.0*zeta*xi*xi+2.0*xi*xi)/(den*den*den*den), -1.875*(2.0*eta+1.0)+0.375*(120.0*eta*eta*xi*xi*zeta*zeta-45.0*zeta*eta*eta*xi*xi-27.0*eta*eta*xi*xi-45.0*eta*eta*zeta*zeta*zeta+127.0*eta*eta*zeta*zeta-119.0*zeta*eta*eta+37.0*eta*eta+10.0*eta*xi*xi*zeta*zeta*zeta-14.0*eta*xi*xi*zeta*zeta-2.0*eta*zeta*xi*xi+6.0*eta*xi*xi-15.0*xi*xi*zeta*zeta*zeta+37.0*xi*xi*zeta*zeta-29.0*zeta*xi*xi+7.0*xi*xi)/(den*den*den*den*den) ) +
                                 c[1] * RealGradient( -4.5*xi-1.5*xi*(3.0*eta*eta*zeta*zeta+4.0*zeta*eta*eta-eta*eta-2.0*zeta*zeta+4.0*zeta-2.0)/(den*den*den*den), -1.5*(12.0*eta-5.0*zeta+3.0)-1.5*(12.0*eta*xi*xi*zeta*zeta+eta*zeta*xi*xi-7.0*eta*xi*xi-2.0*eta*zeta*zeta+4.0*eta*zeta-2.0*eta-2.0*xi*xi*zeta*zeta+4.0*zeta*xi*xi-2.0*xi*xi)/(den*den*den*den), -1.875*(2.0*eta-1.0)-0.375*(120.0*eta*eta*xi*xi*zeta*zeta-45.0*zeta*eta*eta*xi*xi-27.0*eta*eta*xi*xi-45.0*eta*eta*zeta*zeta*zeta+127.0*eta*eta*zeta*zeta-119.0*zeta*eta*eta+37.0*eta*eta-10.0*eta*xi*xi*zeta*zeta*zeta+14.0*eta*xi*xi*zeta*zeta+2.0*eta*zeta*xi*xi-6.0*eta*xi*xi-15.0*xi*xi*zeta*zeta*zeta+37.0*xi*xi*zeta*zeta-29.0*zeta*xi*xi+7.0*xi*xi)/(den*den*den*den*den) ) +
                                 c[2] * RealGradient( 1.5*(12.0*xi+5.0*zeta-3.0)+1.5*(12.0*xi*eta*eta*zeta*zeta+xi*zeta*eta*eta-7.0*xi*eta*eta+2.0*eta*eta*zeta*zeta-4.0*zeta*eta*eta+2.0*eta*eta-2.0*xi*zeta*zeta+4.0*xi*zeta-2.0*xi)/(den*den*den*den), 4.5*eta+1.5*eta*(3.0*xi*xi*zeta*zeta+4.0*zeta*xi*xi-xi*xi-2.0*zeta*zeta+4.0*zeta-2.0)/(den*den*den*den), -1.875*(2.0*xi+1.0)+0.375*(120.0*eta*eta*xi*xi*zeta*zeta-45.0*zeta*eta*eta*xi*xi-27.0*eta*eta*xi*xi+10.0*xi*eta*eta*zeta*zeta*zeta-14.0*xi*eta*eta*zeta*zeta-2.0*xi*zeta*eta*eta+6.0*xi*eta*eta-15.0*eta*eta*zeta*zeta*zeta+37.0*eta*eta*zeta*zeta-29.0*zeta*eta*eta+7.0*eta*eta-45.0*xi*xi*zeta*zeta*zeta+127.0*xi*xi*zeta*zeta-119.0*zeta*xi*xi+37.0*xi*xi)/(den*den*den*den*den) ) +
                                 c[3] * RealGradient( -1.5*(12.0*xi-5.0*zeta+3.0)-1.5*(12.0*xi*eta*eta*zeta*zeta+xi*zeta*eta*eta-7.0*xi*eta*eta-2.0*eta*eta*zeta*zeta+4.0*zeta*eta*eta-2.0*eta*eta-2.0*xi*zeta*zeta+4.0*xi*zeta-2.0*xi)/(den*den*den*den), -4.5*eta-1.5*eta*(3.0*xi*xi*zeta*zeta+4.0*zeta*xi*xi-xi*xi-2.0*zeta*zeta+4.0*zeta-2.0)/(den*den*den*den), -1.875*(2.0*xi-1.0)-0.375*(120.0*eta*eta*xi*xi*zeta*zeta-45.0*zeta*eta*eta*xi*xi-27.0*eta*eta*xi*xi-10.0*xi*eta*eta*zeta*zeta*zeta+14.0*xi*eta*eta*zeta*zeta+2.0*xi*zeta*eta*eta-6.0*xi*eta*eta-15.0*eta*eta*zeta*zeta*zeta+37.0*eta*eta*zeta*zeta-29.0*zeta*eta*eta+7.0*eta*eta-45.0*xi*xi*zeta*zeta*zeta+127.0*xi*xi*zeta*zeta-119.0*zeta*xi*xi+37.0*xi*xi)/(den*den*den*den*den) );
                        }
                      case 20:
                      case 21:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 0, i - 20);
                          return c[0] * RealGradient( -0.5*(6.0*eta-xi-6.0*zeta+3.0)-0.5*eta*(3.0*eta*xi*zeta+3.0*eta*xi+9.0*eta*zeta-9.0*eta+2.0*xi*zeta-2.0*xi)/(den*den*den), 0.0, -0.5*(3.0*xi+1.0)+0.5*eta*(6.0*eta*xi*xi+9.0*eta*xi*zeta+9.0*eta*xi+3.0*eta*zeta-3.0*eta+2.0*zeta*xi*xi-2.0*xi*xi+6.0*xi*zeta-6.0*xi)/(den*den*den) ) +
                                 c[1] * RealGradient( -xi+eta*xi*(3.0*eta*zeta+3.0*eta+2.0*zeta-2.0)/(den*den*den), 0.0, 1.0-eta*(6.0*eta*xi*xi+3.0*eta*zeta-3.0*eta+2.0*zeta*xi*xi-2.0*xi*xi)/(den*den*den) );
                        }
                      case 22:
                      case 23:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 1, i - 22);
                          return c[0] * RealGradient( 0.0, 0.5*(eta+6.0*xi+6.0*zeta-3.0)-0.5*xi*(3.0*eta*xi*zeta+3.0*eta*xi-2.0*eta*zeta+2.0*eta+9.0*xi*zeta-9.0*xi)/(den*den*den), -0.5*(3.0*eta+1.0)+0.5*xi*(6.0*xi*eta*eta-2.0*zeta*eta*eta+2.0*eta*eta+9.0*eta*xi*zeta+9.0*eta*xi-6.0*eta*zeta+6.0*eta+3.0*xi*zeta-3.0*xi)/(den*den*den) ) +
                                 c[1] * RealGradient( 0.0, -eta+eta*xi*(3.0*xi*zeta+3.0*xi-2.0*zeta+2.0)/(den*den*den), 1.0-xi*(6.0*xi*eta*eta-2.0*zeta*eta*eta+2.0*eta*eta+3.0*xi*zeta-3.0*xi)/(den*den*den) );
                        }
                      case 24:
                      case 25:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 2, i - 24);
                          return c[0] * RealGradient( -0.5*(6.0*eta-xi+6.0*zeta-3.0)-0.5*eta*(3.0*eta*xi*zeta+3.0*eta*xi-9.0*eta*zeta+9.0*eta-2.0*xi*zeta+2.0*xi)/(den*den*den), 0.0, 0.5*(3.0*xi-1.0)+0.5*eta*(6.0*eta*xi*xi-9.0*eta*xi*zeta-9.0*eta*xi+3.0*eta*zeta-3.0*eta-2.0*zeta*xi*xi+2.0*xi*xi+6.0*xi*zeta-6.0*xi)/(den*den*den) ) +
                                 c[1] * RealGradient( -xi+eta*xi*(3.0*eta*zeta+3.0*eta-2.0*zeta+2.0)/(den*den*den), 0.0, 1.0-eta*(6.0*eta*xi*xi+3.0*eta*zeta-3.0*eta-2.0*zeta*xi*xi+2.0*xi*xi)/(den*den*den) );
                        }
                      case 26:
                      case 27:
                        {
                          const std::array<Real, 4> c =
                            nedelec_one_face_coefficients(*elem, 3, i - 26);
                          return c[0] * RealGradient( 0.0, 0.5*(eta+6.0*xi-6.0*zeta+3.0)-0.5*xi*(3.0*eta*xi*zeta+3.0*eta*xi+2.0*eta*zeta-2.0*eta-9.0*xi*zeta+9.0*xi)/(den*den*den), 0.5*(3.0*eta-1.0)+0.5*xi*(6.0*xi*eta*eta+2.0*zeta*eta*eta-2.0*eta*eta-9.0*eta*xi*zeta-9.0*eta*xi-6.0*eta*zeta+6.0*eta+3.0*xi*zeta-3.0*xi)/(den*den*den) ) +
                                 c[1] * RealGradient( 0.0, -eta+eta*xi*(3.0*xi*zeta+3.0*xi+2.0*zeta-2.0)/(den*den*den), 1.0-xi*(6.0*xi*eta*eta+2.0*zeta*eta*eta-2.0*eta*eta+3.0*xi*zeta-3.0*xi)/(den*den*den) );
                        }
                      case 28:
                        return RealGradient( -7.5*(2.0*zeta-1.0)-7.5*eta*eta/(den*den), 0.0, 7.5*xi-7.5*xi*eta*eta*(zeta+1.0)/(den*den*den) );
                      case 29:
                        return RealGradient( 0.0, -7.5*(2.0*zeta-1.0)-7.5*xi*xi/(den*den), 7.5*eta-7.5*eta*xi*xi*(zeta+1.0)/(den*den*den) );
                      case 30:
                        return RealGradient( 11.25*xi+3.75*xi*(3.0*eta*eta*zeta*zeta+4.0*zeta*eta*eta-eta*eta-2.0*zeta*zeta+4.0*zeta-2.0)/(den*den*den*den), 11.25*eta+3.75*eta*(3.0*xi*xi*zeta*zeta+4.0*zeta*xi*xi-xi*xi-2.0*zeta*zeta+4.0*zeta-2.0)/(den*den*den*den), -7.5+3.75*(12.0*eta*eta*xi*xi*zeta*zeta-3.0*eta*eta*zeta*zeta*zeta+7.0*eta*eta*zeta*zeta-5.0*zeta*eta*eta+eta*eta-3.0*xi*xi*zeta*zeta*zeta+7.0*xi*xi*zeta*zeta-5.0*zeta*xi*xi+xi*xi)/(den*den*den*den*den) );
                      case 31:
                        return RealGradient( 45.0*xi+7.5*xi*(6.0*eta*eta*zeta*zeta+2.0*zeta*eta*eta-5.0*eta*eta-zeta*zeta+2.0*zeta-1.0)/(den*den*den*den), 11.25*eta+3.75*eta*(3.0*xi*xi*zeta*zeta+4.0*zeta*xi*xi-xi*xi-2.0*zeta*zeta+4.0*zeta-2.0)/(den*den*den*den), -1.875+1.875*(60.0*eta*eta*xi*xi*zeta*zeta-9.0*zeta*eta*eta*xi*xi-27.0*eta*eta*xi*xi-9.0*eta*eta*zeta*zeta*zeta+23.0*eta*eta*zeta*zeta-19.0*zeta*eta*eta+5.0*eta*eta-27.0*xi*xi*zeta*zeta*zeta+77.0*xi*xi*zeta*zeta-73.0*zeta*xi*xi+23.0*xi*xi)/(den*den*den*den*den) );
                      case 32:
                        return RealGradient( 11.25*xi+3.75*xi*(3.0*eta*eta*zeta*zeta+4.0*zeta*eta*eta-eta*eta-2.0*zeta*zeta+4.0*zeta-2.0)/(den*den*den*den), 45.0*eta+7.5*eta*(6.0*xi*xi*zeta*zeta+2.0*zeta*xi*xi-5.0*xi*xi-zeta*zeta+2.0*zeta-1.0)/(den*den*den*den), -1.875+1.875*(60.0*eta*eta*xi*xi*zeta*zeta-9.0*zeta*eta*eta*xi*xi-27.0*eta*eta*xi*xi-27.0*eta*eta*zeta*zeta*zeta+77.0*eta*eta*zeta*zeta-73.0*zeta*eta*eta+23.0*eta*eta-9.0*xi*xi*zeta*zeta*zeta+23.0*xi*xi*zeta*zeta-19.0*zeta*xi*xi+5.0*xi*xi)/(den*den*den*den*den) );
                      case 33:
                        return RealGradient( -56.25*xi-18.75*xi*(3.0*eta*eta*zeta*zeta+4.0*zeta*eta*eta-eta*eta-2.0*zeta*zeta+4.0*zeta-2.0)/(den*den*den*den), -56.25*eta-18.75*eta*(3.0*xi*xi*zeta*zeta+4.0*zeta*xi*xi-xi*xi-2.0*zeta*zeta+4.0*zeta-2.0)/(den*den*den*den), 26.25-3.75*(60.0*eta*eta*xi*xi*zeta*zeta+9.0*zeta*eta*eta*xi*xi-9.0*eta*eta*xi*xi-18.0*eta*eta*zeta*zeta*zeta+44.0*eta*eta*zeta*zeta-34.0*zeta*eta*eta+8.0*eta*eta-18.0*xi*xi*zeta*zeta*zeta+44.0*xi*xi*zeta*zeta-34.0*zeta*xi*xi+8.0*xi*xi)/(den*den*den*den*den) );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      } // switch(ii)
                  } // j = 2

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 3D element type!: " << Utility::enum_to_string(elem->type()));
          }
      }

