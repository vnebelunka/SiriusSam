//
// Created by nebil on 10.09.22.
//

#ifndef SIRIUS_ELECTRIC_H
#define SIRIUS_ELECTRIC_H

#include "Geometry.h"
#include "Integrator.h"

using arma::cx_mat;
using arma::cx_vec;


void calcMatrixE(const Grid &g, double k, cx_mat &M);

cx_vec calcFE(const Grid &g, double k, vec3 Eplr, vec3 v0);

#endif //SIRIUS_ELECTRIC_H
