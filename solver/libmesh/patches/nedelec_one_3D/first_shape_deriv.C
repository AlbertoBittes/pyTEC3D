          case PRISM15:
          case PRISM18:
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
                    switch(i)
                      {
                      case 0:
                        return sign * RealGradient( 0.0, 0.5*(zeta-1.0), 0.0 );
                      case 1:
                        return sign * RealGradient( 0.0, 0.5*(zeta-1.0), 0.0 );
                      case 2:
                        return sign * RealGradient( 0.0, -0.5*(zeta-1.0), 0.0 );
                      case 3:
                        return sign * RealGradient( 0.0, 0.0, 0.5 );
                      case 4:
                        return sign * RealGradient( 0.0, 0.0, -0.5 );
                      case 6:
                        return sign * RealGradient( 0.0, -0.5*(zeta+1.0), 0.0 );
                      case 7:
                        return sign * RealGradient( 0.0, -0.5*(zeta+1.0), 0.0 );
                      case 8:
                        return sign * RealGradient( 0.0, 0.5*(zeta+1.0), 0.0 );
                      case 5:
                        return RealGradient();

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      } // switch(i)
                  } // j = 0

                  // d()/deta
                case 1:
                  {
                    switch(i)
                      {
                      case 0:
                        return sign * RealGradient( -0.5*(zeta-1.0), 0.0, 0.0 );
                      case 1:
                        return sign * RealGradient( -0.5*(zeta-1.0), 0.0, 0.0 );
                      case 2:
                        return sign * RealGradient( 0.5*(zeta-1.0), 0.0, 0.0 );
                      case 3:
                        return sign * RealGradient( 0.0, 0.0, 0.5 );
                      case 5:
                        return sign * RealGradient( 0.0, 0.0, -0.5 );
                      case 6:
                        return sign * RealGradient( 0.5*(zeta+1.0), 0.0, 0.0 );
                      case 7:
                        return sign * RealGradient( 0.5*(zeta+1.0), 0.0, 0.0 );
                      case 8:
                        return sign * RealGradient( -0.5*(zeta+1.0), 0.0, 0.0 );
                      case 4:
                        return RealGradient();

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      } // switch(i)
                  } // j = 1

                  // d()/dzeta
                case 2:
                  {
                    switch(i)
                      {
                      case 0:
                        return sign * RealGradient( -0.5*(eta-1.0), 0.5*xi, 0.0 );
                      case 1:
                        return sign * RealGradient( -0.5*eta, 0.5*xi, 0.0 );
                      case 2:
                        return sign * RealGradient( 0.5*eta, -0.5*(xi-1.0), 0.0 );
                      case 6:
                        return sign * RealGradient( 0.5*(eta-1.0), -0.5*xi, 0.0 );
                      case 7:
                        return sign * RealGradient( 0.5*eta, -0.5*xi, 0.0 );
                      case 8:
                        return sign * RealGradient( -0.5*eta, 0.5*(xi-1.0), 0.0 );
                      case 3:
                      case 4:
                      case 5:
                        return RealGradient();

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      } // switch(i)
                  } // j = 2

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          case PYRAMID13:
          case PYRAMID14:
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
                    switch(i)
                      {
                      case 0:
                        return sign * RealGradient( 0.0, 0.0, -0.25+0.25*eta/den );
                      case 1:
                        return sign * RealGradient( 0.0, -0.25, -0.25*eta/den );
                      case 2:
                        return sign * RealGradient( 0.0, 0.0, 0.25+0.25*eta/den );
                      case 3:
                        return sign * RealGradient( 0.0, 0.25, 0.25*eta/den );
                      case 4:
                        return sign * RealGradient( 0.0, 0.25*zeta/den, 0.25+0.25*eta*(2.0*zeta-1.0)/(den*den) );
                      case 5:
                        return sign * RealGradient( 0.0, 0.25-0.25/den, -0.25-0.25*eta*(2.0*zeta-1.0)/(den*den) );
                      case 6:
                        return sign * RealGradient( 0.0, 0.25*zeta/den, -0.25+0.25*eta*(2.0*zeta-1.0)/(den*den) );
                      case 7:
                        return sign * RealGradient( 0.0, 0.25-0.25/den, 0.25-0.25*eta*(2.0*zeta-1.0)/(den*den) );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      } // switch(i)
                  } // j = 0

                  // d()/deta
                case 1:
                  {
                    switch(i)
                      {
                      case 0:
                        return sign * RealGradient( 0.25, 0.0, 0.25*xi/den );
                      case 1:
                        return sign * RealGradient( 0.0, 0.0, -0.25-0.25*xi/den );
                      case 2:
                        return sign * RealGradient( 0.25, 0.0, 0.25*xi/den );
                      case 3:
                        return sign * RealGradient( 0.0, 0.0, -0.25+0.25*xi/den );
                      case 4:
                        return sign * RealGradient( 0.25*zeta/den, 0.0, 0.25+0.25*xi*(2.0*zeta-1.0)/(den*den) );
                      case 5:
                        return sign * RealGradient( 0.25-0.25/den, 0.0, 0.25-0.25*xi*(2.0*zeta-1.0)/(den*den) );
                      case 6:
                        return sign * RealGradient( 0.25*zeta/den, 0.0, -0.25+0.25*xi*(2.0*zeta-1.0)/(den*den) );
                      case 7:
                        return sign * RealGradient( 0.25-0.25/den, 0.0, -0.25-0.25*xi*(2.0*zeta-1.0)/(den*den) );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      } // switch(i)
                  } // j = 1

                  // d()/dzeta
                case 2:
                  {
                    switch(i)
                      {
                      case 0:
                        return sign * RealGradient( 0.25, 0.0, 0.25*eta*xi/(den*den) );
                      case 1:
                        return sign * RealGradient( 0.0, 0.25, -0.25*eta*xi/(den*den) );
                      case 2:
                        return sign * RealGradient( -0.25, 0.0, 0.25*eta*xi/(den*den) );
                      case 3:
                        return sign * RealGradient( 0.0, 0.25, 0.25*eta*xi/(den*den) );
                      case 4:
                        return sign * RealGradient( -0.25+0.25*eta/(den*den), -0.25+0.25*xi/(den*den), 0.5*eta*xi*zeta/(den*den*den) );
                      case 5:
                        return sign * RealGradient( 0.25-0.25*eta/(den*den), -0.25-0.25*xi/(den*den), -0.5*eta*xi*zeta/(den*den*den) );
                      case 6:
                        return sign * RealGradient( 0.25+0.25*eta/(den*den), 0.25+0.25*xi/(den*den), 0.5*eta*xi*zeta/(den*den*den) );
                      case 7:
                        return sign * RealGradient( -0.25-0.25*eta/(den*den), 0.25-0.25*xi/(den*den), -0.5*eta*xi*zeta/(den*den*den) );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      } // switch(i)
                  } // j = 2

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

