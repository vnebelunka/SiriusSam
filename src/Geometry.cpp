/**
 * \file
 * Source file for geometry operations
 */

#include "../include/Geometry.h"
#include <iostream>


double area(const vec3 &a, const vec3 &b, const vec3 &c) {
    vec3 AB = b - a;
    vec3 AC = c - b;
    vec3 n = cross(AB, AC);
    return 0.5 * norm(n);
}

double area(const Triangle &t) {
    return area(t.a, t.b, t.c);
}

vec3 normal(const vec3 &a, const vec3 &b, const vec3 &c) {
    vec3 AB = b - a;
    vec3 AC = c - b;
    vec3 n = cross(AB, AC);
    return n / norm(n);
}

vec3 normal(const Triangle &t){
    return normal(t.a, t.b, t.c);
}

double norm(const vec3 &v) {
    return sqrt(dot(v, v));
}


double normsqr(const array<complex<double>, 3> &p) {
    double ans = 0;
    for (int i = 0; i < 3; ++i) {
        ans += (p[i] * conj(p[i])).real();
    }
    return ans;
}

double dist(const vec3 &x, const vec3 &y) {
    return norm(x - y);
}

vec3 cross(const vec3 &a, const vec3 &b) {
    return {a.y * b.z - a.z * b.y, b.x * a.z - a.x * b.z, a.x * b.y - a.y * b.x};
}

double dot(const vec3 &a, const vec3 &b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

using vec3c = array<complex<double>, 3>;

vec3c operator+(const vec3c &a, const vec3c &b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

vec3c operator*(const vec3c &a, complex<double> mult) {
    return {a[0] * mult, a[1] * mult, a[2] * mult};
}

vec3c operator-(const vec3c &a, const vec3c &b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

vec3c operator+=(vec3c &a, const vec3c &other) {
    a[0] += other[0], a[1] += other[1], a[2] += other[2];
    return a;
}

vec3c operator*(vec3 const &a, complex<double> mult) {
    return {a.x * mult, a.y * mult, a.z * mult};
}

ostream &operator<<(ostream &os, const vec3 &p) {
    os << p.x << " " << p.y << " " << p.z;
    return os;
}

vec3 vec3::operator+=(const vec3 &other) {
    x += other.x, y += other.y, z += other.z;
    return *this;
}

vec3 vec3::operator-=(const vec3 &other) {
    x -= other.x, y -= other.y, z -= other.z;
    return *this;
}


complex<double> dot(vec3 const &a, const array<complex<double>, 3> &b) {
    return a.x * b[0] + a.y * b[1] + a.z * b[2];
}

complex<double> cdot(vec3c const &a, vec3 const &b){
    return a[0] * b.x + a[1] * b.y + a[2] * b.z;
}

vec3c cross(vec3c const &a, vec3 const &b){
    return {a[1] * b.z - a[2] * b.y, a[2] * b.x - a[0] * b.z, a[0] * b.y - a[1] * b.x};
}

vec3c cross(vec3c const& a, vec3c const &b){
    return {a[1] * b[2] -a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}


std::array<vec3, 4> calcBarCoords(const vec3 &a, const vec3 &b, const vec3 &c) {
    std::array<vec3, 4> b_coords;
    b_coords[0] = a / 3. + b / 3. + c / 3.;
    b_coords[1] = a * (3. / 5) + b / 5. + c / 5.;
    b_coords[2] = a / 5. + b * (3. / 5) + c / 5.;
    b_coords[3] = a / 5. + b / 5 + c * (3. / 5.);
    return b_coords;
}