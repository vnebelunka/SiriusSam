
#ifndef SIRIUSSAM_MAGNETIC_H
#define SIRIUSSAM_MAGNETIC_H

#include "Integrator.h"
#include "Grid.h"
#include "armadillo"

using arma::cx_mat;
using arma::cx_vec;

void calcMatrixM(const Grid &g, double k, cx_mat& M);

void calcFM(const Grid &g, double k, vec3 Eplr, vec3 v0, cx_vec &f);

#endif //SIRIUSSAM_MAGNETIC_H
