
#ifndef SIRIUSSAM_MAGNETIC_H
#define SIRIUSSAM_MAGNETIC_H

#include "Integrator.h"
#include "Grid.h"
#include "armadillo"

using arma::cx_mat;
using arma::cx_vec;

/*
complex<double> intFarM(const MarkedTriangle &tx, const MarkedTriangle &ty, double k);

complex<double> intNearM(MarkedTriangle const& tx, MarkedTriangle const& ty, double k);

complex<double> intFM(const MarkedTriangle &t, double k, vec3 const& Eplr, vec3 const& v0);

complex<double> intEdgeM(const Grid &g, const pair<int, int>& e1, const pair<int, int>& e2, const pair<int, int>& v1,
                         const pair<int, int>& v2, double k);
*/
void calcMatrixM(const Grid &g, double k, cx_mat& M);

cx_vec calcFM(const Grid &g, double k, vec3 Eplr, vec3 v0);

#endif //SIRIUSSAM_MAGNETIC_H
