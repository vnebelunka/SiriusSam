#include "../include/Geometry.h"
#include <iostream>

using vec3 = array<double, 3>;


double area(const Point &a, const Point &b, const Point &c) {
    Point AB = b - a;
    Point AC = c - b;
    auto n = cross(vec3(AB), vec3(AC));
    return 0.5 * norm(n);
}

double area(const Triangle &t) {
    return area(t.a, t.b, t.c);
}

vec3 normal(const Point &a, const Point &b, const Point &c) {
    Point AB = b - a;
    Point AC = c - b;
    auto n = cross(vec3(AB), vec3(AC));
    auto modn = norm(n);
    return {n[0] / modn, n[1] / modn, n[2] / modn};
}

double norm(vec3 v) {
    return sqrt(dot(v, v));
}

double norm(Point p) {
    return norm(vec3(p));
}

double normsqr(vec3c p){
    double ans = 0;
    for(int i = 0; i < 3; ++i){
        ans += (p[i] * conj(p[i])).real();
    }
    return ans;
}

double dist(const Point &x, const Point &y) {
    return norm(vec3(x - y));
}

vec3 cross(vec3 a, vec3 b) {
    return {a[1] * b[2] - a[2] * b[1], b[0] * a[2] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

double dot(const vec3 &a, const vec3 &b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

vector<double> operator*(const vector<double> &ac, const double mult) {
    auto temp = ac;
    for (auto &x: temp) {
        x *= mult;
    }
    return temp;
}

vec3 operator*(const vec3 &ac, const double mult){
    auto temp = ac;
    for (auto &x: temp) {
        x *= mult;
    }
    return temp;
}

vector<double> operator/(const vector<double> &lhs, const double div) {
    auto temp = lhs;
    for (auto &x: temp) {
        x /= div;
    }
    return temp;
}

vec3 operator+(const vec3 &lhs, const vec3 &rhs) {
    vec3 res{lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2]};
    return res;
}

vec3 operator-(const vec3 &lhs, const vec3 &rhs) {
    vec3 res{lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2]};
    return res;
}

using vec3c = array<complex<double>, 3>;

vec3c operator+(const vec3c& a,const vec3c& b){
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

vec3c operator*(const vec3c& a, complex<double> mult){
    return {a[0] * mult, a[1] * mult, a[2] * mult};
}

vec3c operator-(const vec3c& a,const vec3c& b){
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

vec3c operator+=(vec3c&a, const vec3c&other){
    a[0] += other[0];
    a[1] += other[1];
    a[2] += other[2];
    return a;
}

vec3c operator*(Point a, complex<double> mult){
    return {a.x * mult, a.y * mult, a.z * mult};
}

ostream &operator<<(ostream &os, const Point &p) {
    os << p.x << " " << p.y << " " << p.z;
    return os;
}

vec3c operator*(const vec3 a, complex<double> mult) {
    return{a[0] * mult, a[1] * mult, a[2] * mult};
}

complex<double> dot(const array<double, 3> &a, const array<complex<double>, 3> &b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
