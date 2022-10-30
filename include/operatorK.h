
#ifndef SIRIUSSAM_MAGNETIC_H
#define SIRIUSSAM_MAGNETIC_H

#include "Integrator.h"
#include "Grid.h"
#include "armadillo"

using arma::cx_mat;
using arma::cx_vec;

complex<double>
intEdge_e_Ke(const Grid &g, const pair<int, int> &e1, const pair<int, int> &e2, const pair<int, int> &v1,
             const pair<int, int> &v2, double k);

complex<double>
intEdge_e_Einc(const Grid &g, std::pair<int, int> e, std::pair<int, int> v,
               double k, vec3 const& Eplr, vec3 const& v0);

#endif //SIRIUSSAM_MAGNETIC_H
