#include <iostream>
#include <algorithm>
#include "BiQuadraticInterpolation.hpp"

BiQuadraticInterpolation::BiQuadraticInterpolation(std::vector<double> xx,
                                                   std::vector<double> yy,
                                                   std::function<double(double, double)> f) : n(xx.size()-1), m(yy.size()-1), z(xx), t(yy), coeffs(n*m), Interpolation()
                                                   {
                                                       x = std::vector<double>(n);
                                                       y = std::vector<double>(m);
                                                       dx = std::vector<double>(n+1);
                                                       dy = std::vector<double>(m+1);

                                                       for (int i = 0; i < x.size(); i++)
                                                       {
                                                           x[i] = 0.5*(z[i+1] + z[i]);
                                                       }

                                                       for (int j = 0; j < y.size(); j++)
                                                       {
                                                           y[j] = 0.5*(t[j+1] + t[j]);
                                                       }

                                                       for (int i = 1; i < dx.size()-1; i++)
                                                       {
                                                           dx[i] = x[i] - x[i-1];
                                                       }
                                                       dx[0] = dx[1];
                                                       dx[n] = dx[n-1];

                                                       for (int j = 1; j < dy.size()-1; j++)
                                                       {
                                                           dy[j] = y[j] - y[j-1];
                                                       }
                                                       dy[0] = dy[1];
                                                       dy[m] = dy[m-1];

                                                       calcCoeffs(f);
                                                   }


                                                   BiQuadraticInterpolation::BiQuadraticInterpolation(std::vector<double> xx,
                                                                                                      std::vector<double> yy,
                                                                                                      std::function<double(double, double)> f,
                                                                                                      std::function<double(double, double)> fx,
                                                                                                      std::function<double(double, double)> fy,
                                                                                                      std::function<double(double, double)> fxy) : n(xx.size()), m(yy.size()), x(xx), y(yy), coeffs(n*m), Interpolation()
                                                                                                      {
                                                                                                          dx = std::vector<double>(n+1);
                                                                                                          dy = std::vector<double>(m+1);
                                                                                                          z = std::vector<double>(n+1);
                                                                                                          t = std::vector<double>(m+1);

                                                                                                          for (int i = 1; i < dx.size()-1; i++)
                                                                                                          {
                                                                                                              dx[i] = x[i] - x[i-1];
                                                                                                          }
                                                                                                          dx[0] = dx[1];
                                                                                                          dx[n] = dx[n-1];

                                                                                                          for (int j = 1; j < dy.size()-1; j++)
                                                                                                          {
                                                                                                              dy[j] = y[j] - y[j-1];
                                                                                                          }
                                                                                                          dy[0] = dy[1];
                                                                                                          dy[m] = dy[m-1];

                                                                                                          z[0] = x[0] - 0.5*dx[0];
                                                                                                          for (int i = 0; i < z.size()-1; i++)
                                                                                                          {
                                                                                                              z[i+1] = x[i] + 0.5*dx[i+1];
                                                                                                          }

                                                                                                          t[0] = y[0] - 0.5*dy[0];
                                                                                                          for (int j = 0; j < t.size()-1; j++)
                                                                                                          {
                                                                                                              t[j+1] = y[j] + 0.5*dy[j+1];
                                                                                                          }

                                                                                                          calcCoeffs(f, fx, fy, fxy);
                                                                                                      }


                                                                                                      BiQuadraticInterpolation::BiQuadraticInterpolation(std::vector<int> gridSizeX_, std::vector<int> gridSizeY_,
                                                                                                                                                         std::vector<double> boundaryX_, std::vector<double> boundaryY_,
                                                                                                                                                         Transformation transformationX_, Transformation transformationY_,
                                                                                                                                                         std::function<double(double, double)> f) : n(0), m(0), x(0), y(0), coeffs(0),
                                                                                                                                                         Interpolation(gridSizeX_, gridSizeY_, boundaryX_, boundaryY_, transformationX_, transformationY_)
                                                                                                                                                         {
                                                                                                                                                             std::vector<double> boundaryXAux = boundaryX;
                                                                                                                                                             std::vector<double> boundaryYAux = boundaryY;
                                                                                                                                                             for (int i = 0; i < boundaryXAux.size(); i++)
                                                                                                                                                             {
                                                                                                                                                                 boundaryXAux[i] = transform(boundaryXAux[i], transformationX);
                                                                                                                                                             }

                                                                                                                                                             for (int j = 0; j < boundaryYAux.size(); j++)
                                                                                                                                                             {
                                                                                                                                                                 boundaryYAux[j] = transform(boundaryYAux[j], transformationY);
                                                                                                                                                             }

                                                                                                                                                             if(boundaryXAux[0] > boundaryXAux[1])
                                                                                                                                                             {
                                                                                                                                                                 std::reverse(boundaryXAux.begin(), boundaryXAux.end());
                                                                                                                                                             }

                                                                                                                                                             if(boundaryYAux[0] > boundaryYAux[1])
                                                                                                                                                             {
                                                                                                                                                                 std::reverse(boundaryYAux.begin(), boundaryYAux.end());
                                                                                                                                                             }

                                                                                                                                                             double sizeX = 0.0;
                                                                                                                                                             for (int i = 0; i < gridSizeX.size(); i++)
                                                                                                                                                             {
                                                                                                                                                                 sizeX += gridSizeX[i];
                                                                                                                                                             }

                                                                                                                                                             double sizeY = 0.0;
                                                                                                                                                             for (int j = 0; j < gridSizeY.size(); j++)
                                                                                                                                                             {
                                                                                                                                                                 sizeY += gridSizeY[j];
                                                                                                                                                             }

                                                                                                                                                             x = std::vector<double>(sizeX);
                                                                                                                                                             z = std::vector<double>(sizeX+1);

                                                                                                                                                             y = std::vector<double>(sizeY);
                                                                                                                                                             t = std::vector<double>(sizeY+1);

                                                                                                                                                             dz = std::vector<double>(gridSizeX.size());
                                                                                                                                                             dt = std::vector<double>(gridSizeY.size());

                                                                                                                                                             for (int i = 0; i < dz.size(); i++)
                                                                                                                                                             {
                                                                                                                                                                 dz[i] = (boundaryXAux[i+1] - boundaryXAux[i])/gridSizeX[i];
                                                                                                                                                             }

                                                                                                                                                             for (int j = 0; j < dt.size(); j++)
                                                                                                                                                             {
                                                                                                                                                                 dt[j] = (boundaryYAux[j+1] - boundaryYAux[j])/gridSizeY[j];
                                                                                                                                                             }

                                                                                                                                                             int idx = 0;
                                                                                                                                                             z[idx] = backTransform(boundaryXAux[0], transformationX);
                                                                                                                                                             x[idx] = backTransform(boundaryXAux[0]+0.5*dz[0], transformationX);
                                                                                                                                                             idx++;

                                                                                                                                                             for (int ii = 0; ii < gridSizeX.size(); ii++)
                                                                                                                                                             {
                                                                                                                                                                 for (int i = 0; i < gridSizeX[ii]; i++)
                                                                                                                                                                 {
                                                                                                                                                                     z[idx] = backTransform(boundaryXAux[ii] + dz[ii]*(i+1), transformationX);
                                                                                                                                                                     x[idx] = backTransform(boundaryXAux[ii] + dz[ii]*(i+1) + 0.5*dz[ii], transformationX);
                                                                                                                                                                     idx++;
                                                                                                                                                                 }
                                                                                                                                                             }

                                                                                                                                                             idx = 0;
                                                                                                                                                             t[idx] = backTransform(boundaryYAux[0], transformationY);
                                                                                                                                                             y[idx] = backTransform(boundaryYAux[0]+0.5*dt[0], transformationY);
                                                                                                                                                             idx++;

                                                                                                                                                             for (int jj = 0; jj < gridSizeY.size(); jj++)
                                                                                                                                                             {
                                                                                                                                                                 for (int j = 0; j < gridSizeY[jj]; j++)
                                                                                                                                                                 {
                                                                                                                                                                     t[idx] = backTransform(boundaryYAux[jj] + dt[jj]*(j+1), transformationY);
                                                                                                                                                                     y[idx] = backTransform(boundaryYAux[jj] + dt[jj]*(j+1) + 0.5*dt[jj], transformationY);
                                                                                                                                                                     idx++;
                                                                                                                                                                 }
                                                                                                                                                             }

                                                                                                                                                             if(z[0] > z[1])
                                                                                                                                                             {
                                                                                                                                                                 std::reverse(x.begin(), x.end());
                                                                                                                                                                 std::reverse(z.begin(), z.end());
                                                                                                                                                             }
                                                                                                                                                             if(t[0] > t[1])
                                                                                                                                                             {
                                                                                                                                                                 std::reverse(y.begin(), y.end());
                                                                                                                                                                 std::reverse(t.begin(), t.end());
                                                                                                                                                             }


                                                                                                                                                             ////////
                                                                                                                                                             n = x.size();
                                                                                                                                                             m = y.size();
                                                                                                                                                             dx = std::vector<double>(n+1);
                                                                                                                                                             dy = std::vector<double>(m+1);
                                                                                                                                                             coeffs = std::vector<Mat3x3>(m*n);

                                                                                                                                                             for (int i = 1; i < dx.size()-1; i++)
                                                                                                                                                             {
                                                                                                                                                                 dx[i] = x[i] - x[i-1];
                                                                                                                                                             }
                                                                                                                                                             dx[0] = dx[1];
                                                                                                                                                             dx[n] = dx[n-1];

                                                                                                                                                             for (int j = 1; j < dy.size()-1; j++)
                                                                                                                                                             {
                                                                                                                                                                 dy[j] = y[j] - y[j-1];
                                                                                                                                                             }
                                                                                                                                                             dy[0] = dy[1];
                                                                                                                                                             dy[m] = dy[m-1];

                                                                                                                                                             calcCoeffs(f);
                                                                                                                                                         }



                                                                                                                                                         void BiQuadraticInterpolation::calcCoeffs(std::function<double(double, double)> f,
                                                                                                                                                                                                   std::function<double(double, double)> fx,
                                                                                                                                                                                                   std::function<double(double, double)> fy,
                                                                                                                                                                                                   std::function<double(double, double)> fxy)
                                                                                                                                                         {
                                                                                                                                                             Matrix u(n, m);
                                                                                                                                                             Matrix p(n+1, m);
                                                                                                                                                             Matrix q(n, m+1);
                                                                                                                                                             Matrix r(n+1, m+1);

                                                                                                                                                             for (int i = 0; i < n; i++)
                                                                                                                                                             {
                                                                                                                                                                 for (int j = 0; j < m; j++)
                                                                                                                                                                 {
                                                                                                                                                                     u[i][j] = f(x[i], y[j]);
                                                                                                                                                                 }
                                                                                                                                                             }

                                                                                                                                                             for (int j = 0; j < m; j++)
                                                                                                                                                             {
                                                                                                                                                                 p[0][j] = fx(z[0], y[j]);
                                                                                                                                                                 p[n][j] = fx(z[n], y[j]);
                                                                                                                                                             }

                                                                                                                                                             for (int i = 0; i < n; i++)
                                                                                                                                                             {
                                                                                                                                                                 q[i][0] = fy(x[i], t[0]);
                                                                                                                                                                 q[i][m] = fy(x[i], t[m]);
                                                                                                                                                             }

                                                                                                                                                             r[0][0] = fxy(z[0], t[0]);
                                                                                                                                                             r[n][0] = fxy(z[n], t[0]);
                                                                                                                                                             r[0][m] = fxy(z[0], t[m]);
                                                                                                                                                             r[n][m] = fxy(z[n], t[m]);


                                                                                                                                                             /////// solution for p
                                                                                                                                                             for (int j = 0; j < m; j++)
                                                                                                                                                             {
                                                                                                                                                                 std::vector<double> U(n-2);
                                                                                                                                                                 std::vector<double> D(n-1);
                                                                                                                                                                 std::vector<double> L(n-2);
                                                                                                                                                                 std::vector<double> b(n-1);
                                                                                                                                                                 int i = 0;

                                                                                                                                                                 U[i] = 1.0/(dx[i+1] + dx[i+2]);
                                                                                                                                                                 D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));

                                                                                                                                                                 b[i] = (-1.0/(dx[i] + dx[i+1]))*p[0][j] + (4.0/dx[i+1])*((u[i+1][j] - u[i][j])/dx[i+1]);

                                                                                                                                                                 i++;
                                                                                                                                                                 for ( ; i < n-2; i++)
                                                                                                                                                                 {
                                                                                                                                                                     U[i] = 1.0/(dx[i+1] + dx[i+2]);
                                                                                                                                                                     D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));
                                                                                                                                                                     L[i-1] = 1.0/(dx[i] + dx[i+1]);
                                                                                                                                                                     b[i] = (4.0/dx[i+1])*((u[i+1][j] - u[i][j])/dx[i+1]);
                                                                                                                                                                 }

                                                                                                                                                                 D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));
                                                                                                                                                                 L[i-1] = 1.0/(dx[i]+dx[i+1]);
                                                                                                                                                                 b[i] = (-1.0/(dx[i+1] + dx[i+2]))*p[n][j] + (4.0/dx[i+1])*((u[i+1][j] - u[i][j])/dx[i+1]);

                                                                                                                                                                 std::vector<double> pj = solveTridiagonal(L, D, U, b);
                                                                                                                                                                 for (int i = 1; i < pj.size()+1; i++)
                                                                                                                                                                 {
                                                                                                                                                                     p[i][j] = pj[i-1];
                                                                                                                                                                 }
                                                                                                                                                             }

                                                                                                                                                             ///////// solution for q
                                                                                                                                                             for (int i = 0; i < n; i++)
                                                                                                                                                             {
                                                                                                                                                                 std::vector<double> U(m-2);
                                                                                                                                                                 std::vector<double> D(m-1);
                                                                                                                                                                 std::vector<double> L(m-2);
                                                                                                                                                                 std::vector<double> b(m-1);
                                                                                                                                                                 int j = 0;

                                                                                                                                                                 U[j] = 1.0/(dy[j+1] + dy[j+2]);
                                                                                                                                                                 D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));

                                                                                                                                                                 b[j] = (-1.0/(dy[j] + dy[j+1]))*q[i][0] + (4.0/dy[j+1])*((u[i][j+1] - u[i][j])/dy[j+1]);

                                                                                                                                                                 j++;
                                                                                                                                                                 for ( ; j < m-2; j++)
                                                                                                                                                                 {
                                                                                                                                                                     U[j] = 1.0/(dy[j+1] + dy[j+2]);
                                                                                                                                                                     D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));
                                                                                                                                                                     L[j-1] = 1.0/(dy[j] + dy[j+1]);
                                                                                                                                                                     b[j] = (4.0/dy[j+1])*((u[i][j+1] - u[i][j])/dy[j+1]);
                                                                                                                                                                 }

                                                                                                                                                                 D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));
                                                                                                                                                                 L[j-1] = 1.0/(dy[j] + dy[j+1]);
                                                                                                                                                                 b[j] = (-1.0/(dy[j+1] + dy[j+2]))*q[i][m] + (4.0/dy[j+1])*((u[i][j+1] - u[i][j])/dy[j+1]);

                                                                                                                                                                 std::vector<double> qi = solveTridiagonal(L, D, U, b);
                                                                                                                                                                 for (int j = 1; j < qi.size()+1; j++)
                                                                                                                                                                 {
                                                                                                                                                                     q[i][j] = qi[j-1];
                                                                                                                                                                 }
                                                                                                                                                             }


                                                                                                                                                             ///////// solution for ri0 and rim
                                                                                                                                                             for (int j = 0; j < m+1; j += m)
                                                                                                                                                             {
                                                                                                                                                                 std::vector<double> U(n-2);
                                                                                                                                                                 std::vector<double> D(n-1);
                                                                                                                                                                 std::vector<double> L(n-2);
                                                                                                                                                                 std::vector<double> b(n-1);
                                                                                                                                                                 int i = 0;

                                                                                                                                                                 U[i] = 1.0/(dx[i+1] + dx[i+2]);
                                                                                                                                                                 D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));

                                                                                                                                                                 b[i] = (-1.0/(dx[i] + dx[i+1]))*r[i][j] + (4.0/dx[i+1])*((q[i+1][j] - q[i][j])/dx[i+1]);

                                                                                                                                                                 i++;
                                                                                                                                                                 for ( ; i < n-2; i++)
                                                                                                                                                                 {
                                                                                                                                                                     U[i] = 1.0/(dx[i+1] + dx[i+2]);
                                                                                                                                                                     D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));
                                                                                                                                                                     L[i-1] = 1.0/(dx[i] + dx[i+1]);
                                                                                                                                                                     b[i] = (4.0/dx[i+1])*((q[i+1][j] - q[i][j])/dx[i+1]);
                                                                                                                                                                 }

                                                                                                                                                                 D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));
                                                                                                                                                                 L[i-1] = 1.0/(dx[i]+dx[i+1]);
                                                                                                                                                                 b[i] = (-1.0/(dx[i+1] + dx[i+2]))*r[i+2][j] + (4.0/dx[i+1])*((q[i+1][j] - q[i][j])/dx[i+1]);

                                                                                                                                                                 std::vector<double> rj = solveTridiagonal(L, D, U, b);
                                                                                                                                                                 for (int i = 1; i < rj.size()+1; i++)
                                                                                                                                                                 {
                                                                                                                                                                     r[i][j] = rj[i-1];
                                                                                                                                                                 }
                                                                                                                                                             }


                                                                                                                                                             //////solution for r
                                                                                                                                                             for (int i = 0; i < n+1; i++)
                                                                                                                                                             {
                                                                                                                                                                 std::vector<double> U(m-2);
                                                                                                                                                                 std::vector<double> D(m-1);
                                                                                                                                                                 std::vector<double> L(m-2);
                                                                                                                                                                 std::vector<double> b(m-1);
                                                                                                                                                                 int j = 0;

                                                                                                                                                                 U[j] = 1.0/(dy[j+1] + dy[j+2]);
                                                                                                                                                                 D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));

                                                                                                                                                                 b[j] = (-1.0/(dy[j] + dy[j+1]))*r[i][0] + (4.0/dy[j+1])*((p[i][j+1] - p[i][j])/dy[j+1]);

                                                                                                                                                                 j++;
                                                                                                                                                                 for ( ; j < m-2; j++)
                                                                                                                                                                 {
                                                                                                                                                                     U[j] = 1.0/(dy[j+1] + dy[j+2]);
                                                                                                                                                                     D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));
                                                                                                                                                                     L[j-1] = 1.0/(dy[j] + dy[j+1]);
                                                                                                                                                                     b[j] = (4.0/dy[j+1])*((p[i][j+1] - p[i][j])/dy[j+1]);
                                                                                                                                                                 }

                                                                                                                                                                 D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));
                                                                                                                                                                 L[j-1] = 1.0/(dy[j] + dy[j+1]);
                                                                                                                                                                 b[j] = (-1.0/(dy[j+1] + dy[j+2]))*r[i][m] + (4.0/dy[j+1])*((p[i][j+1] - p[i][j])/dy[j+1]);

                                                                                                                                                                 std::vector<double> ri = solveTridiagonal(L, D, U, b);
                                                                                                                                                                 for (int j = 1; j < ri.size()+1; j++)
                                                                                                                                                                 {
                                                                                                                                                                     r[i][j] = ri[j-1];
                                                                                                                                                                 }
                                                                                                                                                             }

                                                                                                                                                             int idx = 0;
                                                                                                                                                             for(int j = 0; j < m; j++)
                                                                                                                                                             {
                                                                                                                                                                 for(int i = 0; i < n; i++)
                                                                                                                                                                 {
                                                                                                                                                                     Mat3x3 VxInv({0 , 1.0, 0,
                                                                                                                                                                         dx[i+1]/(dx[i] + dx[i+1]), 0, dx[i]/(dx[i] + dx[i+1]),
                                                                                                                                                                                  -1.0/(dx[i] + dx[i+1]), 0, 1.0/(dx[i] + dx[i+1])});

                                                                                                                                                                     Mat3x3 VyInv({0 , 1.0, 0,
                                                                                                                                                                         dy[j+1]/(dy[j] + dy[j+1]), 0, dy[j]/(dy[j] + dy[j+1]),
                                                                                                                                                                                  -1.0/(dy[j] + dy[j+1]), 0, 1.0/(dy[j] + dy[j+1])});

                                                                                                                                                                     Mat3x3 C({r[i][j], p[i][j], r[i][j+1],
                                                                                                                                                                         q[i][j], u[i][j], q[i][j+1],
                                                                                                                                                                         r[i+1][j], p[i+1][j], r[i+1][j+1]});

                                                                                                                                                                     coeffs[idx] = VyInv*C*(VyInv.transpose());
                                                                                                                                                                     idx++;
                                                                                                                                                                 }
                                                                                                                                                             }
                                                                                                                                                         }


                                                                                                                                                         ////////////////////

                                                                                                                                                         void BiQuadraticInterpolation::calcCoeffs(std::function<double(double, double)> f)
                                                                                                                                                         {
                                                                                                                                                             Matrix u(n, m);
                                                                                                                                                             Matrix p(n+1, m);
                                                                                                                                                             Matrix q(n, m+1);
                                                                                                                                                             Matrix r(n+1, m+1);

                                                                                                                                                             for (int i = 0; i < n; i++)
                                                                                                                                                             {
                                                                                                                                                                 for (int j = 0; j < m; j++)
                                                                                                                                                                 {
                                                                                                                                                                     u[i][j] = f(x[i], y[j]);
                                                                                                                                                                 }
                                                                                                                                                             }

                                                                                                                                                             const double divider = 100.0;

                                                                                                                                                             for (int j = 0; j < m; j++)
                                                                                                                                                             {
                                                                                                                                                                 p[0][j] = calcDerivativeX(f, z[0], y[j], dx[0]/divider);
                                                                                                                                                                 p[n][j] = calcDerivativeX(f, z[n], y[j], dx[n]/divider);
                                                                                                                                                             }

                                                                                                                                                             for (int i = 0; i < n; i++)
                                                                                                                                                             {
                                                                                                                                                                 q[i][0] = calcDerivativeY(f, x[i], t[0], dy[0]/divider);
                                                                                                                                                                 q[i][m] = calcDerivativeY(f, x[i], t[m], dy[m]/divider);
                                                                                                                                                             }

                                                                                                                                                             r[0][0] = calcDerivativeXY(f, z[0], t[0], dx[0]/divider, dy[0]/divider);
                                                                                                                                                             r[n][0] = calcDerivativeXY(f, z[n], t[0], dx[n]/divider, dy[0]/divider);
                                                                                                                                                             r[0][m] = calcDerivativeXY(f, z[0], t[m], dx[0]/divider, dy[m]/divider);
                                                                                                                                                             r[n][m] = calcDerivativeXY(f, z[n], t[m], dx[n]/divider, dy[m]/divider);


                                                                                                                                                             /////// solution for p
                                                                                                                                                             for (int j = 0; j < m; j++)
                                                                                                                                                             {
                                                                                                                                                                 std::vector<double> U(n-2);
                                                                                                                                                                 std::vector<double> D(n-1);
                                                                                                                                                                 std::vector<double> L(n-2);
                                                                                                                                                                 std::vector<double> b(n-1);
                                                                                                                                                                 int i = 0;

                                                                                                                                                                 U[i] = 1.0/(dx[i+1] + dx[i+2]);
                                                                                                                                                                 D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));

                                                                                                                                                                 b[i] = (-1.0/(dx[i] + dx[i+1]))*p[0][j] + (4.0/dx[i+1])*((u[i+1][j] - u[i][j])/dx[i+1]);

                                                                                                                                                                 i++;
                                                                                                                                                                 for ( ; i < n-2; i++)
                                                                                                                                                                 {
                                                                                                                                                                     U[i] = 1.0/(dx[i+1] + dx[i+2]);
                                                                                                                                                                     D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));
                                                                                                                                                                     L[i-1] = 1.0/(dx[i] + dx[i+1]);
                                                                                                                                                                     b[i] = (4.0/dx[i+1])*((u[i+1][j] - u[i][j])/dx[i+1]);
                                                                                                                                                                 }

                                                                                                                                                                 D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));
                                                                                                                                                                 L[i-1] = 1.0/(dx[i]+dx[i+1]);
                                                                                                                                                                 b[i] = (-1.0/(dx[i+1] + dx[i+2]))*p[n][j] + (4.0/dx[i+1])*((u[i+1][j] - u[i][j])/dx[i+1]);

                                                                                                                                                                 std::vector<double> pj = solveTridiagonal(L, D, U, b);
                                                                                                                                                                 for (int i = 1; i < pj.size()+1; i++)
                                                                                                                                                                 {
                                                                                                                                                                     p[i][j] = pj[i-1];
                                                                                                                                                                 }
                                                                                                                                                             }

                                                                                                                                                             ///////// solution for q
                                                                                                                                                             for (int i = 0; i < n; i++)
                                                                                                                                                             {
                                                                                                                                                                 std::vector<double> U(m-2);
                                                                                                                                                                 std::vector<double> D(m-1);
                                                                                                                                                                 std::vector<double> L(m-2);
                                                                                                                                                                 std::vector<double> b(m-1);
                                                                                                                                                                 int j = 0;

                                                                                                                                                                 U[j] = 1.0/(dy[j+1] + dy[j+2]);
                                                                                                                                                                 D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));

                                                                                                                                                                 b[j] = (-1.0/(dy[j] + dy[j+1]))*q[i][0] + (4.0/dy[j+1])*((u[i][j+1] - u[i][j])/dy[j+1]);

                                                                                                                                                                 j++;
                                                                                                                                                                 for ( ; j < m-2; j++)
                                                                                                                                                                 {
                                                                                                                                                                     U[j] = 1.0/(dy[j+1] + dy[j+2]);
                                                                                                                                                                     D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));
                                                                                                                                                                     L[j-1] = 1.0/(dy[j] + dy[j+1]);
                                                                                                                                                                     b[j] = (4.0/dy[j+1])*((u[i][j+1] - u[i][j])/dy[j+1]);
                                                                                                                                                                 }

                                                                                                                                                                 D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));
                                                                                                                                                                 L[j-1] = 1.0/(dy[j] + dy[j+1]);
                                                                                                                                                                 b[j] = (-1.0/(dy[j+1] + dy[j+2]))*q[i][m] + (4.0/dy[j+1])*((u[i][j+1] - u[i][j])/dy[j+1]);

                                                                                                                                                                 std::vector<double> qi = solveTridiagonal(L, D, U, b);
                                                                                                                                                                 for (int j = 1; j < qi.size()+1; j++)
                                                                                                                                                                 {
                                                                                                                                                                     q[i][j] = qi[j-1];
                                                                                                                                                                 }
                                                                                                                                                             }


                                                                                                                                                             ///////// solution for ri0 and rim
                                                                                                                                                             for (int j = 0; j < m+1; j += m)
                                                                                                                                                             {
                                                                                                                                                                 std::vector<double> U(n-2);
                                                                                                                                                                 std::vector<double> D(n-1);
                                                                                                                                                                 std::vector<double> L(n-2);
                                                                                                                                                                 std::vector<double> b(n-1);
                                                                                                                                                                 int i = 0;

                                                                                                                                                                 U[i] = 1.0/(dx[i+1] + dx[i+2]);
                                                                                                                                                                 D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));

                                                                                                                                                                 b[i] = (-1.0/(dx[i] + dx[i+1]))*r[i][j] + (4.0/dx[i+1])*((q[i+1][j] - q[i][j])/dx[i+1]);

                                                                                                                                                                 i++;
                                                                                                                                                                 for ( ; i < n-2; i++)
                                                                                                                                                                 {
                                                                                                                                                                     U[i] = 1.0/(dx[i+1] + dx[i+2]);
                                                                                                                                                                     D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));
                                                                                                                                                                     L[i-1] = 1.0/(dx[i] + dx[i+1]);
                                                                                                                                                                     b[i] = (4.0/dx[i+1])*((q[i+1][j] - q[i][j])/dx[i+1]);
                                                                                                                                                                 }

                                                                                                                                                                 D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));
                                                                                                                                                                 L[i-1] = 1.0/(dx[i]+dx[i+1]);
                                                                                                                                                                 b[i] = (-1.0/(dx[i+1] + dx[i+2]))*r[i+2][j] + (4.0/dx[i+1])*((q[i+1][j] - q[i][j])/dx[i+1]);

                                                                                                                                                                 std::vector<double> rj = solveTridiagonal(L, D, U, b);
                                                                                                                                                                 for (int i = 1; i < rj.size()+1; i++)
                                                                                                                                                                 {
                                                                                                                                                                     r[i][j] = rj[i-1];
                                                                                                                                                                 }
                                                                                                                                                             }


                                                                                                                                                             //////solution for r
                                                                                                                                                             for (int i = 0; i < n+1; i++)
                                                                                                                                                             {
                                                                                                                                                                 std::vector<double> U(m-2);
                                                                                                                                                                 std::vector<double> D(m-1);
                                                                                                                                                                 std::vector<double> L(m-2);
                                                                                                                                                                 std::vector<double> b(m-1);
                                                                                                                                                                 int j = 0;

                                                                                                                                                                 U[j] = 1.0/(dy[j+1] + dy[j+2]);
                                                                                                                                                                 D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));

                                                                                                                                                                 b[j] = (-1.0/(dy[j] + dy[j+1]))*r[i][0] + (4.0/dy[j+1])*((p[i][j+1] - p[i][j])/dy[j+1]);

                                                                                                                                                                 j++;
                                                                                                                                                                 for ( ; j < m-2; j++)
                                                                                                                                                                 {
                                                                                                                                                                     U[j] = 1.0/(dy[j+1] + dy[j+2]);
                                                                                                                                                                     D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));
                                                                                                                                                                     L[j-1] = 1.0/(dy[j] + dy[j+1]);
                                                                                                                                                                     b[j] = (4.0/dy[j+1])*((p[i][j+1] - p[i][j])/dy[j+1]);
                                                                                                                                                                 }

                                                                                                                                                                 D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));
                                                                                                                                                                 L[j-1] = 1.0/(dy[j] + dy[j+1]);
                                                                                                                                                                 b[j] = (-1.0/(dy[j+1] + dy[j+2]))*r[i][m] + (4.0/dy[j+1])*((p[i][j+1] - p[i][j])/dy[j+1]);

                                                                                                                                                                 std::vector<double> ri = solveTridiagonal(L, D, U, b);
                                                                                                                                                                 for (int j = 1; j < ri.size()+1; j++)
                                                                                                                                                                 {
                                                                                                                                                                     r[i][j] = ri[j-1];
                                                                                                                                                                 }
                                                                                                                                                             }

                                                                                                                                                             int idx = 0;
                                                                                                                                                             for(int j = 0; j < m; j++)
                                                                                                                                                             {
                                                                                                                                                                 for(int i = 0; i < n; i++)
                                                                                                                                                                 {
                                                                                                                                                                     Mat3x3 VxInv({0 , 1.0, 0,
                                                                                                                                                                         dx[i+1]/(dx[i] + dx[i+1]), 0, dx[i]/(dx[i] + dx[i+1]),
                                                                                                                                                                                  -1.0/(dx[i] + dx[i+1]), 0, 1.0/(dx[i] + dx[i+1])});

                                                                                                                                                                     Mat3x3 VyInv({0 , 1.0, 0,
                                                                                                                                                                         dy[j+1]/(dy[j] + dy[j+1]), 0, dy[j]/(dy[j] + dy[j+1]),
                                                                                                                                                                                  -1.0/(dy[j] + dy[j+1]), 0, 1.0/(dy[j] + dy[j+1])});

                                                                                                                                                                     Mat3x3 C({r[i][j], p[i][j], r[i][j+1],
                                                                                                                                                                         q[i][j], u[i][j], q[i][j+1],
                                                                                                                                                                         r[i+1][j], p[i+1][j], r[i+1][j+1]});

                                                                                                                                                                     coeffs[idx] = VyInv*C*(VyInv.transpose());
                                                                                                                                                                     idx++;
                                                                                                                                                                 }
                                                                                                                                                             }
                                                                                                                                                         }


                                                                                                                                                         //////////////////

                                                                                                                                                         std::pair<int, int> BiQuadraticInterpolation::findPosition(double xx, double yy) const
                                                                                                                                                         {
                                                                                                                                                             int i = 0;
                                                                                                                                                             int found = 2;
                                                                                                                                                             while(i < z.size()-1)
                                                                                                                                                             {
                                                                                                                                                                 if(xx >= z[i] && xx <= z[i+1])
                                                                                                                                                                 {
                                                                                                                                                                     found--;
                                                                                                                                                                     break;
                                                                                                                                                                 }
                                                                                                                                                                 i = i + 1;
                                                                                                                                                             }

                                                                                                                                                             int j = 0;
                                                                                                                                                             while(t.size()-1)
                                                                                                                                                             {
                                                                                                                                                                 if(yy >= t[j] && yy <= t[j+1])
                                                                                                                                                                 {
                                                                                                                                                                     found--;
                                                                                                                                                                     break;
                                                                                                                                                                 }
                                                                                                                                                                 j = j + 1;
                                                                                                                                                             }

                                                                                                                                                             if(found == 0)
                                                                                                                                                             {
                                                                                                                                                                 return std::pair<int, int>(i, j);
                                                                                                                                                             }

                                                                                                                                                             std::cout << "ERROR: Interpolation out of range" << std::endl;
                                                                                                                                                             return std::pair<int, int>(0, 0);
                                                                                                                                                         }

                                                                                                                                                         double BiQuadraticInterpolation::calc(double xx, double yy) const
                                                                                                                                                         {
                                                                                                                                                             std::pair<int, int> position = findPosition(xx, yy);
                                                                                                                                                             //std::pair<int, int> position = fastFindPosition(xx, yy);

                                                                                                                                                             double v = xx - x[position.first];
                                                                                                                                                             double w = yy - y[position.second];

                                                                                                                                                             const Mat3x3& coeff = coeffs[position.second*n + position.first];

                                                                                                                                                             return coeff[0][0] + w*(coeff[0][1] + w*coeff[0][2]) + v*(coeff[1][0] + w*(coeff[1][1] + w*coeff[1][2]) + v*(coeff[2][0] + w*(coeff[2][1] + w*coeff[2][2])));
                                                                                                                                                         }

                                                                                                                                                         std::pair<int, int> BiQuadraticInterpolation::fastFindPosition(double xx, double yy) const
                                                                                                                                                         {
                                                                                                                                                             int shiftIdxX = 0;
                                                                                                                                                             int shiftIdxY = 0;

                                                                                                                                                             int ii;
                                                                                                                                                             for (ii = 0; ii < gridSizeX.size(); ii++)
                                                                                                                                                             {
                                                                                                                                                                 if (xx < boundaryX[ii+1]) { break; }
                                                                                                                                                                 shiftIdxX += gridSizeX[ii];
                                                                                                                                                             }
                                                                                                                                                             int jj;
                                                                                                                                                             for (jj = 0; jj < gridSizeY.size(); jj++)
                                                                                                                                                             {
                                                                                                                                                                 if (yy < boundaryY[jj+1]) { break; }
                                                                                                                                                                 shiftIdxY += gridSizeY[jj];
                                                                                                                                                             }

                                                                                                                                                             return std::pair<int, int>(std::floor((abs(transform(xx, transformationX) - transform(boundaryX[ii], transformationX)))/dz[ii]) + shiftIdxX,
                                                                                                                                                                                        std::floor((abs(transform(yy, transformationY) - transform(boundaryY[jj], transformationY)))/dt[jj]) + shiftIdxY);
                                                                                                                                                         }

                                                                                                                                                         double BiQuadraticInterpolation::calcFastFind(double xx, double yy) const
                                                                                                                                                         {
                                                                                                                                                             std::pair<int, int> position = fastFindPosition(xx, yy);
                                                                                                                                                             /*std::pair<int, int> position2 = findPosition(xx, yy);
                                                                                                                                                              *
                                                                                                                                                              *   if (position != position2)
                                                                                                                                                              *   {
                                                                                                                                                              *       std::cout << position.first << ", " << position.second << " ; " <<position2.first << " , " <<position2.second << " x: " << xx << " y: " << yy <<  std::endl;
                                                                                                                                                         }*/

                                                                                                                                                             double v = xx - x[position.first];
                                                                                                                                                             double w = yy - y[position.second];

                                                                                                                                                             const Mat3x3& coeff = coeffs[position.second*n + position.first];

                                                                                                                                                             return coeff[0][0] + w*(coeff[0][1] + w*coeff[0][2]) + v*(coeff[1][0] + w*(coeff[1][1] + w*coeff[1][2]) + v*(coeff[2][0] + w*(coeff[2][1] + w*coeff[2][2])));
                                                                                                                                                         }

                                                                                                                                                         double BiQuadraticInterpolation::calcFastFindDiffX(double xx, double yy) const
                                                                                                                                                         {
                                                                                                                                                             std::pair<int, int> position = fastFindPosition(xx, yy);
                                                                                                                                                             double v = xx - x[position.first];
                                                                                                                                                             double w = yy - y[position.second];

                                                                                                                                                             const Mat3x3& coeff = coeffs[position.second*n + position.first];

                                                                                                                                                             return coeff[1][0] + w*(coeff[1][1] + w*coeff[1][2]) + v*2.0*(coeff[2][0] + w*(coeff[2][1] + w*coeff[2][2]));
                                                                                                                                                         }

                                                                                                                                                         double BiQuadraticInterpolation::calcFastFindDiffY(double xx, double yy) const
                                                                                                                                                         {
                                                                                                                                                             std::pair<int, int> position = fastFindPosition(xx, yy);
                                                                                                                                                             double v = xx - x[position.first];
                                                                                                                                                             double w = yy - y[position.second];

                                                                                                                                                             const Mat3x3& coeff = coeffs[position.second*n + position.first];

                                                                                                                                                             return coeff[0][1] + v*(coeff[1][1] + v*coeff[2][1]) + w*2.0*(coeff[0][2] + v*(coeff[1][2] + v*coeff[2][2]));
                                                                                                                                                         }

                                                                                                                                                         double BiQuadraticInterpolation::calcInverseX(double zz, double yy, double guessXX) const
                                                                                                                                                         {
                                                                                                                                                             std::pair<int, int> position = findPosition(guessXX, yy);

                                                                                                                                                             Mat3x3 coeff = coeffs[position.second*n + position.first];

                                                                                                                                                             double dy = yy - y[position.second];

                                                                                                                                                             double A = coeff[2][0] + dy*(coeff[2][1] + coeff[2][2]*dy);
                                                                                                                                                             double B = coeff[1][0] + dy*(coeff[1][1] + coeff[1][2]*dy);
                                                                                                                                                             double C = coeff[0][0] + dy*(coeff[0][1] + coeff[0][2]*dy) - zz;

                                                                                                                                                             return (-B + sqrt(B*B - 4.0*A*C))/(2.0*A) + x[position.first];
                                                                                                                                                         }

                                                                                                                                                         double BiQuadraticInterpolation::calcInverseY(double xx, double zz, double guessYY) const
                                                                                                                                                         {
                                                                                                                                                             std::pair<int, int> position = findPosition(xx, guessYY);

                                                                                                                                                             Mat3x3 coeff = coeffs[position.second*n + position.first];

                                                                                                                                                             double dx = xx - x[position.first];

                                                                                                                                                             double A = coeff[2][0] + dx*(coeff[2][1] + coeff[2][2]*dx);
                                                                                                                                                             double B = coeff[1][0] + dx*(coeff[1][1] + coeff[1][2]*dx);
                                                                                                                                                             double C = coeff[0][0] + dx*(coeff[0][1] + coeff[0][2]*dx) - zz;

                                                                                                                                                             return (-B + sqrt(B*B - 4.0*A*C))/(2.0*A) + y[position.second];
                                                                                                                                                         }


                                                                                                                                                         std::vector<double> BiQuadraticInterpolation::solveTridiagonal(const std::vector<double>& L,
                                                                                                                                                                                                                        const std::vector<double>& D,
                                                                                                                                                                                                                        const std::vector<double>& U,
                                                                                                                                                                                                                        const std::vector<double>& b) const
                                                                                                                                                                                                                        {
                                                                                                                                                                                                                            int n = b.size();
                                                                                                                                                                                                                            std::vector<double> out(n);

                                                                                                                                                                                                                            std::vector<double> UStar(n-1, 0.0);
                                                                                                                                                                                                                            std::vector<double> bStar(n, 0.0);

                                                                                                                                                                                                                            UStar[0] = U[0] / D[0];
                                                                                                                                                                                                                            bStar[0] = b[0] / D[0];

                                                                                                                                                                                                                            for (int i = 1; i < n-1; i++)
                                                                                                                                                                                                                            {
                                                                                                                                                                                                                                double m = 1.0/(D[i] - L[i-1]*UStar[i-1]);
                                                                                                                                                                                                                                UStar[i] = U[i] * m;
                                                                                                                                                                                                                                bStar[i] = (b[i] - L[i-1]*bStar[i-1])*m;
                                                                                                                                                                                                                            }
                                                                                                                                                                                                                            bStar[n-1] = (b[n-1] - L[n-2]*bStar[n-2])*(1.0/(D[n-1] - L[n-2] * UStar[n-2]));

                                                                                                                                                                                                                            out[n-1] = bStar[n-1];

                                                                                                                                                                                                                            for (int i = n-2; i >= 0; i--)
                                                                                                                                                                                                                            {
                                                                                                                                                                                                                                out[i] = bStar[i] - UStar[i]*out[i+1];
                                                                                                                                                                                                                            }

                                                                                                                                                                                                                            return out;
                                                                                                                                                                                                                        }

                                                                                                                                                                                                                        double BiQuadraticInterpolation::calcDerivativeX(std::function<double(double, double)> f, double x, double y, double step) const
                                                                                                                                                                                                                        {
                                                                                                                                                                                                                            return (f(x+step, y) - f(x-step, y))/(2.0*step);
                                                                                                                                                                                                                        }

                                                                                                                                                                                                                        double BiQuadraticInterpolation::calcDerivativeY(std::function<double(double, double)> f, double x, double y, double step) const
                                                                                                                                                                                                                        {
                                                                                                                                                                                                                            return (f(x, y+step) - f(x, y-step))/(2.0*step);
                                                                                                                                                                                                                        }

                                                                                                                                                                                                                        double BiQuadraticInterpolation::calcDerivativeXY(std::function<double(double, double)> f, double x, double y, double stepX, double stepY) const
                                                                                                                                                                                                                        {
                                                                                                                                                                                                                            return (f(x+stepX, y+stepY) - f(x+stepX, y-stepY) - f(x-stepX, y+stepY) + f(x-stepX, y-stepY))/(4.0*stepX*stepY);
                                                                                                                                                                                                                        }
