#include "Electric.h"
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
complex<double> ker_operator(vec3 const&x, vec3 const&y, MarkedTriangle const& tx, MarkedTriangle const& ty, const double k){
    vec3c gF = gradF(x, y, k);
    vec3c ker_inner = cross(gF, e(ty, y));
    return cdot(ker_inner, en(tx, x));
}


/*
 * (en, R e)
 */
static inline
complex<double> int_en_Re(MarkedTriangle const& tx, MarkedTriangle const& ty, double k){
    if(tx == ty){
        return 0;
    }
    return integrateGauss(tx, ty, &ker_operator, k);
}


static inline
complex<double> int2(MarkedTriangle const& tx, MarkedTriangle const& ty){
    const array<vec3, 4> &x = tx.barCoords;
    complex<double> ans = 0;
    for(int i = 0; i < 4; ++i){
        ans += 0.5 * dot(e(tx, x[i]), e(ty, x[i])) * w[i];
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
    complex<double> ans2;
    if(txPlus == tyPlus){
        ans2 += int2(txPlus, tyPlus);
    }
    if(txPlus == tyMinus){
        ans2 -= int2(txPlus, tyMinus);
    }
    if(txMinus == tyPlus){
        ans2 -= int2(txMinus, tyPlus);
    }
    if(txMinus == tyMinus){
        ans2 += int2(txMinus, tyMinus);
    }

    return ans2 + ans;
}


// (en_x(x), E_plr) e^{i k (v0, x)}
static inline // (en, E_inc)
complex<double> kerF(vec3 const& x, MarkedTriangle const& t, vec3 const& Eplr, vec3 const& v0, double k){
    static complex<double> i(0., 1.);
    return -dot(en(t, x), Eplr) * exp(i * k * dot(v0, vec3(x)));
}

static inline // (en, E_inc)
complex<double> intF(MarkedTriangle const& t, double k, vec3 const& Eplr, vec3 const& v0){
    return integrateGauss<MarkedTriangle const &, vec3 const &, vec3 const &, double>(t, &kerF, t, Eplr, v0, k);
}

complex<double>
intEdge_en_Einc(const Grid &g, std::pair<int, int> e, std::pair<int, int> v,
                     double k, vec3 const& Eplr, vec3 const& v0){
    MarkedTriangle tPlus(g.triangles.find({e.first, e.second, v.first})->second);
    MarkedTriangle tMinus(g.triangles.find({e.first, e.second, v.second})->second);
    return intF(tPlus, k, Eplr, v0) - intF(tMinus, k, Eplr, v0);
}