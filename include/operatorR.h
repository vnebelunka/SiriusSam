//
// Created by nebil on 10.09.22.
//

#ifndef SIRIUS_OPERATORR_H
#define SIRIUS_OPERATORR_H

#include "Geometry.h"
#include "Integrator.h"

using arma::cx_mat;
using arma::cx_vec;


complex<double>
intEdge_en_Re(const Grid &g, const pair<int, int> &e1, const pair<int, int> &e2, const pair<int, int> &v1,
                              const pair<int, int> &v2, double k);
complex<double>
intEdge_e_Ren(const Grid &g, const pair<int, int> &e1, const pair<int, int> &e2, const pair<int, int> &v1,
              const pair<int, int> &v2, double k);

complex<double>
intEdge_en_Einc(const Grid &g, std::pair<int, int> e, std::pair<int, int> v,
                double k, vec3 const& Eplr, vec3 const& v0);

complex<double> intEdge_e1_e2(const Grid &g, const pair<int, int> &e1, const pair<int, int> &e2, const pair<int, int> &v1,
                              const pair<int, int> &v2);

#endif //SIRIUS_OPERATORR_H
