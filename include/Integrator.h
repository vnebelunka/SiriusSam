#ifndef SIRIUS_INTEGRATOR_H
#define SIRIUS_INTEGRATOR_H

#include <array>
#include "Geometry.h"
#include "Grid.h"

//TODO: разделить на интегратор и сами интегралы по ядрам

//TODO: попробовать привязать к Eigen.

//TODO: переписать квадратуру на функцию с произвольным числом доп параметров и обьединить все квадратуры.

using vec3 = std::array<double, 3>;

std::array<Point, 4> calcBarCoords(Triangle t);

complex<double> integrateGaus(const Triangle &t, complex<double> f(Point x));

complex<double> surfIntegral1D(const std::vector<Triangle>& grid, complex<double> f(Point));

complex<double> integrateGaus(const MarkedTriangle &tx, const MarkedTriangle &ty, double k,
                              complex<double> (*f)(const Point, const Point, double, const MarkedTriangle&, const MarkedTriangle&));

double integral1Divr(const Triangle& t, const Point& a);

complex<double> intFar(const MarkedTriangle &tx, const MarkedTriangle &ty, double k);

complex<double> intNear(const MarkedTriangle &tx, const MarkedTriangle &ty, double k);

complex<double> intF(const MarkedTriangle &t, double k, vec3 Eplr, vec3 v0);

array<complex<double>, 3> intSigma(const MarkedTriangle &t, vec3 tau);

double calcSigma(const Grid &g, vector<complex<double>>& j, double k, vec3 tau);

#endif //SIRIUS_INTEGRATOR_H