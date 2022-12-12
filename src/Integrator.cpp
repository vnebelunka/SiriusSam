#include "../include/Geometry.h"
#include "../include/Integrator.h"
#include <complex>


complex<double> surfIntegral1D(const vector<Triangle> &grid, complex<double> f(vec3 const&)){
    complex<double> integral = 0;
    for(auto const& t: grid) {
        auto temp = integrateGauss(t, f);
        integral += temp;
    }
    return integral;
}

//TODO: нижние 3 функции делают +- одно и то же.

complex<double> calcJ(const Grid &g, arma::cx_vec const& j, const pair<int, int> e1, int v){
    pair<int, int> p = {min(e1.first, e1.second), max(e1.first, e1.second)};
    auto it = g.edges.find(p);
    if(it == g.edges.end()){
        exit(2);
    }
    auto vertices = it->second;
    if(vertices.second == -1){
        return 0;
    }

    auto it2 = g.edges_inner_enum.find(p);
    if(it2 == g.edges_inner_enum.end()){
        exit(3);
    }
    size_t i = it2->second;

    if(vertices.first == v){
        return j[i];
    }
    if(vertices.second == v){
        return -j[i];
    }
    exit(1); // Сюда попадать не стоит
}

vec3c calcFlow(const Grid &g, arma::cx_vec const& j, iTriangle const& t, const char mod){
    vec3 A(g.points[t.iv1]), B(g.points[t.iv2]), C(g.points[t.iv3]), M((A + B + C) / 3);
    pair<int, int> edges[3];
    edges[0] = {t.iv1, t.iv2}, edges[1] = {t.iv2, t.iv3}, edges[2] = {t.iv1, t.iv3};
    complex<double> coeffs[3];
    coeffs[0] = calcJ(g, j, edges[0], t.iv3);
    coeffs[1] = calcJ(g, j, edges[1], t.iv1);
    coeffs[2] = calcJ(g, j, edges[2], t.iv2);
    vec3c e1, e2, e3;
    if(mod == 'n') {
        e1 = en(MarkedTriangle(A, B, C), M) * coeffs[0];
        e2 = en(MarkedTriangle(B, C, A), M) * coeffs[1];
        e3 = en(MarkedTriangle(A, C, B), M) * coeffs[2];
    } else if(mod == 'e') {
        e1 = e(MarkedTriangle(A, B, C), M) * coeffs[0];
        e2 = e(MarkedTriangle(B, C, A), M) * coeffs[1];
        e3 = e(MarkedTriangle(A, C, B), M) * coeffs[2];
    }
    return e1 + e2 + e3;
}

double calcSigmaE(const Grid &g, arma::cx_vec const& je, double k, vec3 const& tau){
    vec3c ans;
    for(auto &t: g.itriangles){
        pair<int, int> edges[3];
        edges[0] = {t.iv1, t.iv2}, edges[1] = {t.iv2, t.iv3}, edges[2] = {t.iv1, t.iv3};
        complex<double> coeffs[3];
        coeffs[0] = calcJ(g, je, edges[0], t.iv3);
        coeffs[1] = calcJ(g, je, edges[1], t.iv1);
        coeffs[2] = calcJ(g, je, edges[2], t.iv2);
        vec3 A = g.points[t.iv1], B = g.points[t.iv2], C = g.points[t.iv3];
        Triangle sigmai (A, B, C);
        vec3c intSigma;

        const array<vec3, 4> &x = sigmai.barCoords;
        for(int i = 0; i < 4; ++i){
            vec3c g = (C - x[i]) * coeffs[0] + (A - x[i]) * coeffs[1] + (B - x[i]) * coeffs[2];
            g = g * (1. / sigmai.S);
            vec3c kervec = g - tau * dot(tau, g);
            vec3c ker = kervec * k * k * exp(-complex(0., 1.) * k * dot(tau, vec3(x[i])));
            intSigma += ker * barweights_4[i];
        }
        intSigma = intSigma * sigmai.S;
        ans += intSigma;
    }
    return normsqr(ans) / (4 * M_PI);
}

double calcSigmaM(const Grid &g, arma::cx_vec const& jm, double k, vec3 const& tau){
    vec3c ans;
    for(auto &t: g.itriangles){
        pair<int, int> edges[3];
        edges[0] = {t.iv1, t.iv2}, edges[1] = {t.iv2, t.iv3}, edges[2] = {t.iv1, t.iv3};
        complex<double> coeffs[3];
        coeffs[0] = calcJ(g, jm, edges[0], t.iv3); // c - marked
        coeffs[1] = calcJ(g, jm, edges[1], t.iv1); // a - marked
        coeffs[2] = calcJ(g, jm, edges[2], t.iv2); // b - marked
        vec3 A = g.points[t.iv1], B = g.points[t.iv2], C = g.points[t.iv3];
        Triangle sigmai (A, B, C);
        vec3c intSigma;

        const array<vec3, 4> &x = sigmai.barCoords;
        for(int i = 0; i < 4; ++i){
            vec3c g = (C - x[i]) * coeffs[0] + (A - x[i]) * coeffs[1] + (B - x[i]) * coeffs[2];
            g = g * (1. / sigmai.S);
            vec3c kervec = cross(g, tau);
            vec3c ker = kervec * k * exp(-complex(0., 1.) * k * dot(tau, vec3(x[i])));
            intSigma += ker * barweights_4[i];
        }
        intSigma = intSigma * sigmai.S;
        ans += intSigma;
    }
    return normsqr(ans) / (4 * M_PI);
}

double calcSigmaEM(const Grid &g, arma::cx_vec const& je, arma::cx_vec const& jm, double k, double w, double eps, vec3 const& tau){
    vec3c ans;
    complex<double> i(0, 1);
    for(auto &t: g.itriangles){
        pair<int, int> edges[3];
        edges[0] = {t.iv1, t.iv2}, edges[1] = {t.iv2, t.iv3}, edges[2] = {t.iv1, t.iv3};
        complex<double> coeffse[3], coeffsm[3];

        vec3 A = g.points[t.iv1], B = g.points[t.iv2], C = g.points[t.iv3];
        Triangle sigmai (A, B, C);
        vec3c intSigma;

        if(!je.empty()) {
            coeffse[0] = calcJ(g, je, edges[0], t.iv3); // c - marked
            coeffse[1] = calcJ(g, je, edges[1], t.iv1); // a - marked
            coeffse[2] = calcJ(g, je, edges[2], t.iv2); // b - marked
        }
        if(!jm.empty()) {
            coeffsm[0] = calcJ(g, jm, edges[0], t.iv3); // c - marked
            coeffsm[1] = calcJ(g, jm, edges[1], t.iv1); // a - marked
            coeffsm[2] = calcJ(g, jm, edges[2], t.iv2); // b - marked
        }
        const array<vec3, 4> &x = sigmai.barCoords;
        for(int m = 0; m < 4; ++m){
            vec3c kervecm, kervece;
            vec3c ker;
            //calcSigmaM part
            if(!jm.empty()) {
                vec3c gm = (C - x[m]) * coeffsm[0] + (A - x[m]) * coeffsm[1] + (B - x[m]) * coeffsm[2];
                gm = gm * (1. / sigmai.S);
                kervecm = cross(gm, tau) * k * i;
                //ker = kervecm * k * exp(-i * k * dot(tau, vec3(x[m])));
                //intSigma = ker * barweights_4[m] * (i / (4 * M_PI));
            }
            //E part
            if(!je.empty()) {
                vec3c ge = (C - x[m]) * coeffse[0] + (A - x[m]) * coeffse[1] + (B - x[m]) * coeffse[2];
                ge = ge * (1. / sigmai.S);
                kervece = ge - tau * dot(tau, ge) * (k * k * i / (w * eps));
                //ker = kervec * k * k * exp(-i * k * dot(tau, vec3(x[m])));
                intSigma += ker * barweights_4[m] * (i / (4 * M_PI * w * eps));
            }
            vec3c kervec = kervece + kervecm;
            intSigma += kervec * exp(-i * k * dot(tau, vec3(x[m]))) * (1 / (4 * M_PI));
        }
        intSigma = intSigma * sigmai.S;
        ans += intSigma;
    }
    return normsqr(ans) * (4 * M_PI);
}

void calcTotalFlow(const Grid &g, arma::cx_vec const& j, const char* gridFname, const char* fieldRFname,
                   const char* fieldIFname, const char* normalFname, char mod='e'){
    FILE *fg = fopen(gridFname, "w");
    FILE *fr = fopen(fieldRFname, "w");
    FILE *fi = fopen(fieldIFname, "w");
    FILE *fn = fopen(normalFname, "w");
    for(auto &t: g.itriangles){
        std::array<int, 3> e({t.iv1, t.iv2, t.iv3});
        std::sort(e.begin(), e.end());

        auto it = g.triangles.find({e[0], e[1], e[2]});
        if(it == g.triangles.end()){
            exit(1);
        }
        MarkedTriangle sigmai(it->second);
        vec3 M = (sigmai.a + sigmai.b + sigmai.c) / 3;
        vec3c f = calcFlow(g, j, t, mod);
        vec3 n = sigmai.norm;
        vec3 fre(f[0].real(), f[1].real(), f[2].real());
        vec3 fim(f[0].imag(), f[1].imag(), f[2].imag());
        fprintf(fg, "%lf %lf %lf\n", M.x, M.y, M.z);
        fprintf(fr, "%lf %lf %lf\n", fre.x, fre.y, fre.z);
        fprintf(fi, "%lf %lf %lf\n", fim.x, fim.y, fim.z);
        fprintf(fn, "%lf %lf %lf\n", n.x, n.y, n.z);
    }
    fclose(fg), fclose(fr), fclose(fi), fclose(fn);
}