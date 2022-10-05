#include "Electric.h"
#include "armadillo"
#include <cassert>
#include "progressbar.h"


static
vec3c gradF(vec3 const& x, vec3 const& y, const double k){
    static complex<double> i = complex<double>(0, 1);
    double r = dist(x, y);
    complex<double> mult = exp(i * k * r) * (i * k * r - 1.) / (r * r);
    return (x - y) * (mult / (4 * M_PI * r));
}

static
complex<double> ker_operator(vec3 const&x, vec3 const&y, MarkedTriangle const& tx, MarkedTriangle const& ty, const double k){
    vec3c gF = gradF(x, y, k);
    vec3c ker_inner = cross(gF, e(ty, y));
    return cdot(ker_inner, en(tx, x));
}

static inline
complex<double> intOperator(MarkedTriangle const& tx, MarkedTriangle const& ty, double k){
    if(tx == ty){
        return 0;
    }
    return integrateGauss(tx, ty, &ker_operator, k);
}

/*static
complex<double> ker2(vec3 const &x, MarkedTriangle const &tx){
    return dot(e(tx, x), e(ty, y));
}*/

complex<double> int2(MarkedTriangle const& tx, MarkedTriangle const& ty){
    const array<vec3, 4> &x = tx.barCoords;
    complex<double> ans = 0;
    for(int i = 0; i < 4; ++i){
        ans += dot(e(tx, x[i]), e(ty, x[i])) * w[i];
    }
    ans *= tx.S;
    return ans;
}

static
complex<double> intEdge(const Grid &g, const pair<int, int> &e1, const pair<int, int> &e2, const pair<int, int> &v1,
                         const pair<int, int> &v2, double k){
    MarkedTriangle txPlus(g.triangles.find({e1.first, e1.second, v1.first})->second);
    MarkedTriangle txMinus(g.triangles.find({e1.first, e1.second, v1.second})->second);
    MarkedTriangle tyPlus(g.triangles.find({e2.first, e2.second, v2.first})->second);
    MarkedTriangle tyMinus(g.triangles.find({e2.first, e2.second, v2.second})->second);
    complex<double> ans = intOperator(txPlus, tyPlus, k) + intOperator(txMinus, tyMinus, k);
    ans -= intOperator(txMinus, tyPlus, k) + intOperator(txPlus, tyMinus, k);
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

void calcMatrixE(const Grid &g, double k, cx_mat &M) {
    int i = 0, j;
    progressbar p(int(g.edges_inner_enum.size()));
    for(auto [e1, v1]: g.edges){
        j = 0;
        if(v1.second == -1){
            continue;
        }
        for(auto [e2, v2]: g.edges){
            if(v2.second == -1){
                continue;
            }
            complex<double> temp = intEdge(g, e1, e2, v1, v2, k);
            M(i,j) = temp;
            ++j;
        }
        ++i;
        p.update();
    }
    std::cerr << std::endl;
}

// (en_x(x), E_plr) e^{i k (v0, x)}
static
complex<double> kerF(vec3 const& x, MarkedTriangle const& t, vec3 const& Eplr, vec3 const& v0, double k){
    static complex<double> i(0., 1.);
    return dot(e(t, x), Eplr) * exp(i * k * dot(v0, vec3(x)));
}

static inline
complex<double> intF(MarkedTriangle const& t, double k, vec3 const& Eplr, vec3 const& v0){
    return integrateGauss<MarkedTriangle const &, vec3 const &, vec3 const &, double>(t, &kerF, t, Eplr, v0, k);
}

void calcFE(const Grid &g, double k, vec3 Eplr, vec3 v0, cx_vec &f) {
    int i = 0;
    progressbar p(g.edges_inner_enum.size());
    for(auto [e, v]: g.edges){
        if(v.second == -1){
            continue;
        }
        MarkedTriangle tPlus(g.points[e.first], g.points[e.second], g.points[v.first]);
        MarkedTriangle tMinus(g.points[e.first], g.points[e.second], g.points[v.second]);
        auto temp = intF(tPlus, k, Eplr, v0) - intF(tMinus, k, Eplr, v0);
        f[i] = temp;
        ++i;
        p.update();
    }
    std::cerr<<std::endl;
}