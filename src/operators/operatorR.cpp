#include "operatorR.h"
#include "armadillo"
#include <cassert>
#include "progressbar.h"


static inline
vec3c gradF(vec3 const& x, vec3 const& y, const double k){
    static complex<double> i = complex<double>(0, 1);
    double r = dist(x, y);
    complex<double> mult = exp(i * k * r) * (i * k * r - 1.) / (r * r * r);
    return (x - y) * (mult / (4 * M_PI));
}

/*
 * (en, R e)
 */
static inline
complex<double> ker_en_Re(vec3 const&x, vec3 const&y, MarkedTriangle const& tx, MarkedTriangle const& ty, const double k){
    vec3c gF = gradF(x, y, k);
    vec3c ker_inner = cross(gF, e(ty, y));
    return cdot(ker_inner, en(tx, x));
}

static inline
complex<double> ker_e_Ren(vec3 const&x, vec3 const&y, MarkedTriangle const& tx, MarkedTriangle const& ty, const double k){
    vec3c gf = gradF(x, y, k);
    vec3c ker_inner = cross(gf, en(ty, y));
    return cdot(ker_inner, en(tx, x));
}

static inline // (en, Re)
complex<double> int_en_Re(MarkedTriangle const& tx, MarkedTriangle const& ty, double k){
    if(tx == ty){
        return 0;
    }
    return integrateGauss(tx, ty, &ker_en_Re, k);
}

static inline
complex<double> int_e_Ren(MarkedTriangle const& tx, MarkedTriangle const& ty, double k){
    return integrateGauss(tx, ty, &ker_e_Ren, k);
}

// (en, R e)
complex<double> intEdge_e_Ren(const Grid &g, const pair<int, int> &e1, const pair<int, int> &e2, const pair<int, int> &v1,
                              const pair<int, int> &v2, double k){
    MarkedTriangle txPlus(g.triangles.find({e1.first, e1.second, v1.first})->second);
    MarkedTriangle txMinus(g.triangles.find({e1.first, e1.second, v1.second})->second);
    MarkedTriangle tyPlus(g.triangles.find({e2.first, e2.second, v2.first})->second);
    MarkedTriangle tyMinus(g.triangles.find({e2.first, e2.second, v2.second})->second);
    complex<double> ans = int_en_Re(txPlus, tyPlus, k) + int_en_Re(txMinus, tyMinus, k);
    ans -= int_e_Ren(txMinus, tyPlus, k) + int_e_Ren(txPlus, tyMinus, k);
    return ans;
}


static inline // (e1, e2)
complex<double> int_e1_e2(MarkedTriangle const& tx, MarkedTriangle const& ty){
    if(tx != ty){
        return 0;
    }
    const array<vec3, 4> &x = tx.barCoords;
    complex<double> ans = 0;
    for(int i = 0; i < 4; ++i){
        ans += dot(e(tx, x[i]), e(ty, x[i])) * w[i];
    }
    ans *= tx.S;
    return ans;
}

// (en, R e)
complex<double> intEdge_en_Re(const Grid &g, const pair<int, int> &e1, const pair<int, int> &e2, const pair<int, int> &v1,
                              const pair<int, int> &v2, double k){
    MarkedTriangle txPlus(g.triangles.find({e1.first, e1.second, v1.first})->second);
    MarkedTriangle txMinus(g.triangles.find({e1.first, e1.second, v1.second})->second);
    MarkedTriangle tyPlus(g.triangles.find({e2.first, e2.second, v2.first})->second);
    MarkedTriangle tyMinus(g.triangles.find({e2.first, e2.second, v2.second})->second);
    complex<double> ans = int_en_Re(txPlus, tyPlus, k) + int_en_Re(txMinus, tyMinus, k);
    ans -= int_en_Re(txMinus, tyPlus, k) + int_en_Re(txPlus, tyMinus, k);
    return ans;
}

complex<double> intEdge_e1_e2(const Grid &g, const pair<int, int> &e1, const pair<int, int> &e2, const pair<int, int> &v1,
                              const pair<int, int> &v2){
    MarkedTriangle txPlus(g.triangles.find({e1.first, e1.second, v1.first})->second);
    MarkedTriangle txMinus(g.triangles.find({e1.first, e1.second, v1.second})->second);
    MarkedTriangle tyPlus(g.triangles.find({e2.first, e2.second, v2.first})->second);
    MarkedTriangle tyMinus(g.triangles.find({e2.first, e2.second, v2.second})->second);
    return int_e1_e2(txPlus, tyPlus) + int_e1_e2(txMinus, tyMinus) - int_e1_e2(txMinus, tyPlus) - int_e1_e2(txPlus, tyMinus);
}


// (en_x(x), E_plr) e^{i k (v0, x)}
static inline // (en, E_inc)
complex<double> ker_en_Einc(vec3 const& x, MarkedTriangle const& t, vec3 const& Eplr, vec3 const& v0, double k){
    static complex<double> i(0., 1.);
    return dot(en(t, x), Eplr) * exp(i * k * dot(v0, vec3(x)));
}

static inline // (en, E_inc)
complex<double> int_en_Einc(MarkedTriangle const& t, double k, vec3 const& Eplr, vec3 const& v0){
    return integrateGauss<MarkedTriangle const &, vec3 const &, vec3 const &, double>(t, &ker_en_Einc, t, Eplr, v0, k);
}

complex<double>
intEdge_en_Einc(const Grid &g, std::pair<int, int> e, std::pair<int, int> v,
                     double k, vec3 const& Eplr, vec3 const& v0){
    MarkedTriangle tPlus(g.triangles.find({e.first, e.second, v.first})->second);
    MarkedTriangle tMinus(g.triangles.find({e.first, e.second, v.second})->second);
    return int_en_Einc(tPlus, k, Eplr, v0) - int_en_Einc(tMinus, k, Eplr, v0);
}