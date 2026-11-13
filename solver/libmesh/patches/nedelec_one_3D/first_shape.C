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
              switch(i)
                {
                case 0:
                  return sign * RealGradient( -0.5*(eta-1.0)*(zeta-1.0), 0.5*xi*(zeta-1.0), 0.0 );
                case 1:
                  return sign * RealGradient( -0.5*eta*(zeta-1.0), 0.5*xi*(zeta-1.0), 0.0 );
                case 2:
                  return sign * RealGradient( 0.5*eta*(zeta-1.0), -0.5*(xi-1.0)*(zeta-1.0), 0.0 );
                case 3:
                  return sign * RealGradient( 0.0, 0.0, 0.5*(eta+xi-1.0) );
                case 4:
                  return sign * RealGradient( 0.0, 0.0, -0.5*xi );
                case 5:
                  return sign * RealGradient( 0.0, 0.0, -0.5*eta );
                case 6:
                  return sign * RealGradient( 0.5*(zeta+1.0)*(eta-1.0), -0.5*xi*(zeta+1.0), 0.0 );
                case 7:
                  return sign * RealGradient( 0.5*eta*(zeta+1.0), -0.5*xi*(zeta+1.0), 0.0 );
                case 8:
                  return sign * RealGradient( -0.5*eta*(zeta+1.0), 0.5*(zeta+1.0)*(xi-1.0), 0.0 );

                default:
                  libmesh_error_msg("Invalid i = " << i);
                } // switch(i)
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

              // The shape functions are rational in den = 1-zeta. At the apex,
              // where den vanishes, they only have limits along each direction:
              // take those along the pyramid axis.
              if (1.0 - zeta < 1.e-6)
                {
                  switch(i)
                    {
                    case 4:
                      return sign * RealGradient( -0.25, -0.25, -0.25 );
                    case 5:
                      return sign * RealGradient( 0.25, -0.25, -0.25 );
                    case 6:
                      return sign * RealGradient( 0.25, 0.25, -0.25 );
                    case 7:
                      return sign * RealGradient( -0.25, 0.25, -0.25 );
                    case 0:
                    case 1:
                    case 2:
                    case 3:
                      return RealGradient();

                    default:
                      libmesh_error_msg("Invalid i = " << i);
                    } // switch(i)
                }

              const Real den = 1.0 - zeta;

              switch(i)
                {
                case 0:
                  return sign * RealGradient( 0.25*(eta+zeta-1.0), 0.0, -0.25*xi+0.25*eta*xi/den );
                case 1:
                  return sign * RealGradient( 0.0, -0.25*(xi-zeta+1.0), -0.25*eta-0.25*eta*xi/den );
                case 2:
                  return sign * RealGradient( 0.25*(eta-zeta+1.0), 0.0, 0.25*xi+0.25*eta*xi/den );
                case 3:
                  return sign * RealGradient( 0.0, 0.25*(xi+zeta-1.0), -0.25*eta+0.25*eta*xi/den );
                case 4:
                  return sign * RealGradient( 0.25*zeta*(eta+zeta-1.0)/den, 0.25*zeta*(xi+zeta-1.0)/den, 0.25*(eta+xi-1.0)+0.25*eta*xi*(2.0*zeta-1.0)/(den*den) );
                case 5:
                  return sign * RealGradient( 0.25*(eta+zeta)-0.25*eta/den, 0.25*(xi-zeta)-0.25*xi/den, 0.25*(eta-xi-1.0)-0.25*eta*xi*(2.0*zeta-1.0)/(den*den) );
                case 6:
                  return sign * RealGradient( 0.25*zeta*(eta-zeta+1.0)/den, 0.25*zeta*(xi-zeta+1.0)/den, -0.25*(eta+xi+1.0)+0.25*eta*xi*(2.0*zeta-1.0)/(den*den) );
                case 7:
                  return sign * RealGradient( 0.25*(eta-zeta)-0.25*eta/den, 0.25*(xi+zeta)-0.25*xi/den, -0.25*(eta-xi+1.0)-0.25*eta*xi*(2.0*zeta-1.0)/(den*den) );

                default:
                  libmesh_error_msg("Invalid i = " << i);
                } // switch(i)
            }

