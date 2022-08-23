#include "Integrator.h"

static const std::array<double, 4> w ({-9./16, 25./48, 25./48, 25./48});

template<typename ... Args>
complex<double> integrateGaus(MarkedTriangle const& tx, MarkedTriangle const& ty,
                              complex<double> (*f)(Point const&x, Point const&y, MarkedTriangle const& tx, MarkedTriangle const& ty, Args...), Args... args){
    const array<Point, 4> x = tx.t.barCoords, y = ty.t.barCoords;
    complex<double> ans = 0;
    for(int i = 0; i < 4; ++i){
        for(int j = 0; j < 4; ++j){
            ans += f(x[i], y[j], tx, ty,args...) * w[i] * w[j];
        }
    }
    ans *= ty.t.S;
    ans *= tx.t.S;
    return ans;
}