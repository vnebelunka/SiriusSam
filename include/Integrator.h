/**
 * \file
 * Header file for integrator functions
 */
#ifndef SIRIUS_INTEGRATOR_H
#define SIRIUS_INTEGRATOR_H

#include <array>
#include "Geometry.h"
#include "Grid.h"
#include "armadillo"

//TODO: разделить на интегратор и сами интегралы по ядрам


/**
 * Integrate function f on grid g based on Gauss quad formulas
 * @param grid surface of integral
 * @param f kernel function
 * @return \f$ \int\limits_{grid} f dS \f$
 */
complex<double> surfIntegral1D(const std::vector<Triangle>& grid, complex<double> f(vec3 const&));

/**
 * Integrate function with template arguments on triangle surface
 * @tparam Args meta argument types for function
 * @param t triangle surface
 * @param f kernel function f(x, args) [with meta arguments]
 * @param args meta arguments for function
 * @return \f$ \int\limits_{\Delta} f(x, args) dx \f$
 */
template<typename ... Args>
complex<double> integrateGaus(Triangle const& t, complex<double> (*f)(vec3 const& x, Args... args), Args... args);


/**
 * double integral of f function with meta args on pair of triangles
 * @tparam Args meta argument types for function
 * @param tx first triangle surface
 * @param ty second triangle surface
 * @param f kernel function f(x, y, args) [with meta args]
 * @param args meta argument types for function
 * @return \f$ \iint\limits_{\Delta_x, \Delta_y} f(x, y, args) dx dy \f$
 */
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