#ifndef SIRIUS_INTEGRATOR_H
#define SIRIUS_INTEGRATOR_H

#include <array>
#include "Geometry.h"
#include "Grid.h"
#include "armadillo"

//TODO: разделить на интегратор и сами интегралы по ядрам



complex<double> surfIntegral1D(const std::vector<Triangle>& grid, complex<double> f(vec3 const&));


template<typename ... Args>
complex<double> integrateGaus(Triangle const& t, complex<double> (*f)(vec3 const& x, Args... args), Args... args);

template<typename ... Args>
complex<double> integrateGaus(MarkedTriangle const& tx, MarkedTriangle const& ty,
                              complex<double> (*f)(vec3 const&x, vec3 const&y, MarkedTriangle const&, MarkedTriangle const&, Args...),
                                      Args... args);

double integral1Divr(const Triangle& t, const vec3& a);

complex<double> intFar(const MarkedTriangle &tx, const MarkedTriangle &ty, double k);

complex<double> intNear(MarkedTriangle const& tx, MarkedTriangle const& ty, double k);

complex<double> intF(const MarkedTriangle &t, double k, vec3 const& Eplr, vec3 const& v0);

double calcSigma(const Grid &g, arma::cx_vec const& j, double k, vec3 const& tau);

#include "../src/Integrator.tpp"

#endif //SIRIUS_INTEGRATOR_H