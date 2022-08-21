#include "../include/Geometry.h"
#include "../include/Integrator.h"
#include <complex>

using vec3 = std::array<double, 3>;

static const std::array<double, 4> w ({-9./16, 25./48, 25./48, 25./48});

std::array<Point, 4> calcBarCoords(Triangle t) {
    std::array<Point, 4> b_coords;
    b_coords[0] = t.a / 3. + t.b / 3. + t.c / 3.;
    b_coords[1] = t.a * (3. / 5) + t.b / 5. + t.c / 5.;
    b_coords[2] = t.a / 5. + t.b * (3. / 5) + t.c / 5.;
    b_coords[3] = t.a / 5. + t.b / 5 + t.c * (3. / 5.);
    return b_coords;
}

complex<double> integrateGaus(const Triangle &t, complex<double> (*f)(Point)) {
    std::array<Point, 4> x = calcBarCoords(t);
    complex<double> ans = 0;
    for(int i = 0; i < 4; ++i){
        ans += f(x[i]) * w[i];
    }
    ans *= area(t.a, t.b, t.c);
    return ans;
}


complex<double> integrateGaus(const MarkedTriangle &tx, const MarkedTriangle &ty, double k,
                              complex<double> (*f)(const Point, const Point, double,
                                                    const MarkedTriangle&, const MarkedTriangle&)) {
    std::array<Point, 4> x = calcBarCoords(tx.t), y = calcBarCoords(ty.t);
    complex<double> ans = 0;
    for(int i = 0; i < 4; ++i){
        for(int j = 0; j < 4; ++j){
            ans += f(x[i], y[j], k, tx, ty) * w[i] * w[j];
        }
    }
    ans *= ty.S;
    ans *= tx.S;
    return ans;
}

complex<double> surfIntegral1D(const vector<Triangle> &grid, complex<double> f(Point)){
    complex<double> integral = 0;
    for(auto t: grid) {
        auto temp = integrateGaus(t, f);
        integral += temp;
    }
    return integral;
}

// compute rotation
static void lartgp(double x, double y, double &cs, double &sn, double &r) {
    double tmp = sqrt(x * x + y * y);
    cs = x / tmp;
    sn = y / tmp;
    r = tmp;
}

// apply rotation
static void rot(double &x, double &y, double cs, double sn) {
    double tmp = x * cs + y * sn;
    y = y * cs - x * sn;
    x = tmp;
}

static double Integral1Divr_inv(std::array<vec3, 3> ABC, vec3 D){
    double p0, l, l_pl, ratio, r, ad, proj;
    double n[3], edges[3][3], Q[3][2], DABC[7][3], c[3], s[3];
    int ed,  emn;

    // calculate vectors collinear with to edges of the triangle.

    for (int i = 0; i < 3; ++i) {
        edges[i][0] = ABC[i][1] - ABC[i][0]; // AB
        edges[i][1] = ABC[i][2] - ABC[i][1]; // BC
        edges[i][2] = ABC[i][0] - ABC[i][2]; // CA
    }

    // Compute edge length to select largest one and normalize
    // Also compute direction from D to vertices of triangle
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            DABC[j][i] = ABC[j][i] - D[j];           // DA, DB, DC
        }
        // |DA|, |DB|, |DC|
        DABC[3][i] = sqrt(DABC[0][i] * DABC[0][i] + DABC[1][i] * DABC[1][i] + DABC[2][i] * DABC[2][i]);
        Q[i][0] = edges[i][0];
        Q[i][1] = edges[i][1];
        n[i] = 0;
    }

    DABC[4][0] = DABC[3][1]; // |DB|
    DABC[4][1] = DABC[3][2]; // |DC|
    DABC[4][2] = DABC[3][0]; // |DA|

    // Compute QR with pivoting for triangle edges AB, BC
    lartgp(Q[0][0], Q[1][0], c[0], s[0], Q[0][0]);



    lartgp(Q[0][0], Q[2][0], c[1], s[1], Q[0][0]);
    rot(Q[0][1], Q[1][1], c[0], s[0]);
    rot(Q[0][1], Q[2][1], c[1], s[1]);
    lartgp(Q[1][1], Q[2][1], c[2], s[2], Q[1][1]);

    // S_ABC = Q(1, 1) * Q(2, 2) / 2

    // As edges are linearly dependent Q e_3 is orthogonal to the plane containing triangle
    // This vector is unit and can be used as normal
    n[2] = 1.0;

    // Multiplication of Q by e_3 itself
    rot(n[1], n[2], c[2], -s[2]);
    rot(n[0], n[2], c[1], -s[1]);
    rot(n[0], n[1], c[0], -s[0]);

    // Now we will construct QR decomposition using rotations
    // (in fact, reflections can be used as well, but rotations guarantees correct orientation of columns in matrix Q)
    // for matrices [n edges(:, 1)] = Q_1 R_1, [n edges(:, 2)] = Q_2 R_2 and [n edges(:, 3)] = Q_3 R_3
    // It will give us orhonormal basises Q_1, Q_2, Q_3 as it needs
    // Q_1^T DABC(:, 1), Q_2^T DABC(:, 2), Q_3^T DABD(:, 3) and [Q_1^T DABC(:, 2)]_2, [Q_2^T DABC(:, 3)]_2 and [Q_3^T DABC(:, 1)]_2
    // It is easy to note, that each decomposition requires 3 rotations and first two of them is the same for all the three matrices

    // Compute rotation to zero 2nd element of normal
    lartgp(n[0], n[1], c[0], s[0], n[0]); //First rotation to save

    for (int i = 0; i < 3; ++i) {
        // Apply this rotation to edges
        // call drot(3, edges(1, 1), 3, edges(2, 1), 3, c(1), s(1))
        rot(edges[0][i], edges[1][i], c[0], s[0]);    // AB, BC, CA

        // Also we can apply it to target vectors
        // call drot(3, DABC(1, 1), 7, DABC(2, 1), 7, c(1), s(1))
        rot(DABC[0][i], DABC[1][i], c[0], s[0]);      // XA, XB, XC
    }
    // The same here to zero 3rd element of normal
    lartgp(n[0], n[2], c[0], s[0], n[0]); //Second rotation to save
    for (int i = 0; i < 3; ++i) {
        rot(edges[0][i], edges[2][i],  c[0], s[0]);
        rot(DABC[0][i], DABC[2][i], c[0], s[0]);
    }

    // Now compute last rotation for each edges
    lartgp(edges[1][0], edges[2][0], c[0], s[0], edges[1][0]); //Third rotation for AB
    lartgp(edges[1][1], edges[2][1], c[1], s[1], edges[1][1]); //Third rotation for BC
    lartgp(edges[1][2], edges[2][2], c[2], s[2], edges[1][2]); //Third rotation for CA
    // Compute [Q_1^T DABC(:, 2)]_2, [Q_2^T DABC(:, 3)]_2 and [Q_3^T DABC(:, 1)]_2
    DABC[5][0] = c[0] * DABC[1][1] + s[0] * DABC[2][1]; // (XB, AB)
    DABC[5][1] = c[1] * DABC[1][2] + s[1] * DABC[2][2]; // (XC, BC)
    DABC[5][2] = c[2] * DABC[1][0] + s[2] * DABC[2][0]; // (XA, CA)
    DABC[6][0] = -s[0] * DABC[1][1] + c[0] * DABC[2][1]; // (XB, AB_ort)
    DABC[6][1] = -s[1] * DABC[1][2] + c[1] * DABC[2][2]; // (XC, BC_ort)
    DABC[6][2] = -s[2] * DABC[1][0] + c[2] * DABC[2][0]; // (XA, CA_ort)
    // And now apply last rotations
    rot(DABC[1][0], DABC[2][0], c[0], s[0]);
    rot(DABC[1][1], DABC[2][1], c[1], s[1]);
    rot(DABC[1][2], DABC[2][2], c[2], s[2]);

    // Now DABC contains all the values we want
    //
    // loop by edges of triangle ABC:
    //
    double int_1r = 0;
    proj = 0;
    for (ed = 0; ed < 3; ++ed) {
        ad = 0;
        /*
         * calculate the part of integral for function 1/r that connected
         * with edge(ver)
         */

        r = DABC[3][ed];

        if (std::abs(r) > 0 && std::abs(DABC[4][ed]) > 0) {
            //if p0 is very small all parts are small, so putting them zero
            // Projection to normal to edge in plane.
            // Minus because of left orientation (originally normal was in the middle and now
            // it is swapped with the edge)
            p0 = -DABC[2][ed] / r;

            // Projection to edge
            l = DABC[1][ed] / r;

            // ratio of r_pl and r
            ratio = DABC[4][ed] / r;

            // Projection of next vector to the edge
            l_pl = DABC[5][ed] / DABC[4][ed];
            // Note, that r_pl = 0 means that point D is the same as the next vertice, so direction from the current vertice to D is
            // collinear to the edge. Thus othogonal part is 0, so p0 = 0.
            // r = 0 also means p0 = 0.
            // l = -1 or l_pl = -1 means that vector from point D to the current vertice or to the next vertice is is collinear
            // to the edge connecting those two vertices, so orthogonal part of vector from D to current point is zero, so p_0 = 0
            // So, problems here are possible only if p0 = 0, and we have already checked this

            if (l > -1 && l_pl > -1) {
                ad = p0 * std::log(ratio * (1 + l_pl) /  (1 + l));
            }

            // proj = (XA, XB) * |XC| + (XB, XC) * |XA| + (XC, XA) * |XB|
            emn = (ed + 2) % 3; // emn == ed -1 mod 3
            proj = proj + (DABC[0][ed] * DABC[0][ed] + DABC[1][ed] * DABC[5][ed] + DABC[2][ed] * DABC[6][ed]) * DABC[3][emn];

            int_1r = int_1r + r * ad;
        }
    }

    c[0] = std::abs(DABC[0][0]) * Q[0][0] * Q[1][1]; // (XA, XB, XC)
    c[1] = DABC[3][0] * DABC[3][1] * DABC[3][2] + proj; // |XA| * |XB| * |XC| + ...
    if (c[0] > 0 || std::abs(c[1]) > 0) {
        int_1r = int_1r - std::abs(DABC[0][0]) * 2 * std::atan2(c[0], c[1]);
    }
    return int_1r;
}

double integral1Divr(const Triangle& t, const Point& a) {
    std::array<vec3, 3> abc;
    vec3 d({a.x, a.y, a.z});
    abc[0][0] = t.a.x, abc[0][1] = t.b.x, abc[0][2] = t.c.x;
    abc[1][0] = t.a.y, abc[1][1] = t.b.y, abc[1][2] = t.c.y;
    abc[2][0] = t.a.z, abc[2][1] = t.b.z, abc[2][2] = t.c.z;
    return Integral1Divr_inv(abc, d);
}

vec3 e(MarkedTriangle t, Point x){
    return vec3((t.C - x) / t.S);
}

static complex<double> kerFar(const Point x, const Point y, double k, const MarkedTriangle& tx, const MarkedTriangle& ty){
    double r = dist(x, y);
    const complex<double> i(0, 1);
    complex<double>F = exp(i * k * r) / r;
    double Dx = -2 / tx.S, Dy = -2 / ty.S;
    vec3 ex = e(tx, x), ey = e(ty, y);
    return F * (k * k * dot(ex, ey) - Dx * Dy);
}

complex<double> intFar(const MarkedTriangle &tx, const MarkedTriangle &ty, double k){
    return integrateGaus(tx, ty, k, &kerFar);
}


/*
 * [[k^2 (Cx - x, Cy - y) - 4] e^{ikr - 1} / r ] - k^2 / |x - y|dy dx
 */
static complex<double> ker_near_1(const Point x, const Point y, double k, const MarkedTriangle& tx, const MarkedTriangle& ty){
    double r = dist(x, y);
    const complex<double> i(0, 1);
    complex<double> F;
    if(r < std::numeric_limits<double>::epsilon() * norm(y)){
        F = i * k;
    } else {
        F = (exp(i * k * r) - complex(1., 0.)) / r;
    }
    return F * (k * k * dot(vec3(tx.C - x), vec3(ty.C - y)) - 4.) - k * k / 2. * dist(x, y);
}

static complex<double> int_near_1(const MarkedTriangle &tx, const MarkedTriangle &ty, double k){
    return integrateGaus(tx, ty, k, &ker_near_1) / (tx.S * ty.S);
}

static complex<double> int_near_2(const MarkedTriangle& tx, const MarkedTriangle& ty, double k){
    std::array<Point, 4> x = calcBarCoords(tx.t);
    complex<double> ans = 0;
    for(int i = 0; i < 4; ++i){
            double integral_inner = integral1Divr(ty.t, x[i]);
            double ker = (k * k / 2 * (dot(vec3(x[i]), vec3(x[i] - ty.C)) + dot(vec3(ty.C), vec3(tx.C - x[i]))));
            ker -= 2;
            ans += ker * integral_inner * w[i];
    }
    //ans *= area(tx.t.a, tx.t.b, tx.t.c);
    ans /= ty.S;
    return ans;
}

complex<double> intNear(const MarkedTriangle &tx, const MarkedTriangle &ty, double k) {
    auto i1 = int_near_1(tx, ty, k);
    auto i2 = int_near_2(tx, ty, k);
    auto i3 = int_near_2(ty, tx, k);
    return i1 + i2 + i3;
}


// (e_x(x), E_plr) e^{i k (v0, x)}
complex<double> intF(const MarkedTriangle &t, double k, vec3 Eplr, vec3 v0) {
    auto x = calcBarCoords(t.t);
    complex<double> ans = 0;
    for (int i = 0; i < 4; ++i) {
        ans += dot(e(t, x[i]), Eplr) * exp(complex(0., 1.) * k * dot(v0, vec3(x[i]))) * w[i];
    }
    ans *= t.S;
    return ans;
}

using vec3c = array<complex<double>, 3>;

array<complex<double>, 3> intSigma(const MarkedTriangle &t, const vec3 tau){
    auto x = calcBarCoords(t.t);
    vec3c ans;
    complex<double> i(0, 1);
    for (int k = 0; k < 4; ++k) {
        vec3 ei = e(t, x[k]);
        vec3 vecKer = ei - tau * dot(ei, tau);
        vec3c vecKerC;
        vecKerC[0] = {vecKer[0], 0}, vecKerC[1] = {vecKer[1], 0}, vecKerC[2] = {vecKer[2], 0};
        ans += vecKerC * exp(-i * dot(tau, vec3(x[k]))) * w[k];
    }
    ans = ans * t.S;
    return ans;
}

complex<double> calcJ(const Grid &g, vector<complex<double>>& j, const pair<int, int> e1, int v){
    pair<int, int> p = {min(e1.first, e1.second), max(e1.first, e1.second)};
    auto it = g.edges.find(p);
    size_t i = std::distance( g.edges.begin(), it);
    auto vertices = it->second;
    if(vertices.first == v){
        return j[i];
    }
    if(vertices.second == v){
        return -j[i];
    }
    exit(1); // Сюда попадать не стоит
}


double calcSigma(const Grid &g, vector<complex<double>>& j, double k, vec3 tau){
    vec3c ans;
    for(auto &t: g.triangles){
        pair<int, int> edges[3];
        edges[0] = {t.iv1, t.iv2}, edges[1] = {t.iv2, t.iv3}, edges[2] = {t.iv1, t.iv3};
        complex<double> coeffs[3];
        coeffs[0] = calcJ(g, j, edges[0], t.iv3);
        coeffs[1] = calcJ(g, j, edges[1], t.iv1);
        coeffs[2] = calcJ(g, j, edges[2], t.iv2);
        Triangle sigmai (g.points[t.iv1], g.points[t.iv2], g.points[t.iv3]);
        vec3c intSigma;
        auto x = calcBarCoords(sigmai);
        for(int i = 0; i < 4; ++i){
            vec3c g = (sigmai.c - x[i]) * coeffs[0] + (sigmai.a - x[i]) * coeffs[1] + (sigmai.b - x[i]) * coeffs[2];
            g = g * (1. / area(sigmai));
            vec3c kervec = g - tau * dot(tau, g);
            vec3c ker = kervec * k * k * exp(-complex(0., 1.) * k * dot(tau, vec3(x[i])));
            intSigma += ker * w[i];
        }
        intSigma = intSigma * area((sigmai));
        ans += intSigma;
    }
    return normsqr(ans) * 4 * M_PI;
}