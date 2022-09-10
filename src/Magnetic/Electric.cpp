#include "Electric.h"

static
vec3c gradF(vec3 const& x, vec3 const& y, const double k){
    static complex<double> i = complex<double>(0, 1);
    double r = dist(x, y);
    complex<double> mult = exp(i * k * r) * (i * k * r - 1.) / (r * r *r);
    return {x.x * mult, x.y * mult, x.z * mult};
}

static
double ker(vec3 const&x, vec3 const&y, MarkedTriangle const& tx, MarkedTriangle const& ty, const double k){
    vec3c gF = gradF(x, y, k);
    return 0.;
}