#include "Integrator.h"

static const std::array<double, 4> w ({-9./16, 25./48, 25./48, 25./48});

template<typename ... Args>
complex<double> integrateGaus(MarkedTriangle const& tx, MarkedTriangle const& ty,
                              complex<double> (*f)(Point const&x, Point const&y, MarkedTriangle const& tx, MarkedTriangle const& ty, Args...), Args... args){
    const array<Point, 4> x = tx.barCoords, y = ty.barCoords;
    complex<double> ans = 0;
    for(int i = 0; i < 4; ++i){
        for(int j = 0; j < 4; ++j){
            ans += f(x[i], y[j], tx, ty,args...) * w[i] * w[j];
        }
    }
    ans *= ty.S;
    ans *= tx.S;
    return ans;
}

template<typename ... Args>
complex<double> integrateGaus(Triangle const& t, complex<double> (*f)(Point const& x, Args... args),  Args... args) {
    const array<Point, 4> &x = t.barCoords;
    complex<double> ans = 0;
    for(int i = 0; i < 4; ++i){
        ans += f(x[i], args...) * w[i];
    }
    ans *= t.S;
    return ans;
}