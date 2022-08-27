#ifndef SIRIUS_GEOMETRY_H
#define SIRIUS_GEOMETRY_H

#include <vector>
#include <cmath>
#include <array>
#include <complex>
#include <armadillo>

using namespace std;

using vec3 = std::array<double, 3>;

//TODO: разделить на файл с примитивами и работой с векторами + перейти к шаблонам.

//TODO: Triangle и Marked triangle можно обьединить скорее всего

//TODO: мб стоит переписать на Eigen

struct Point {
    double x, y, z;

    Point(double x1 = 0, double y1 = 0, double z1 = 0) {
        x = x1; y = y1, z = z1;
    }

    Point operator+=(const Point& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    Point& operator=(const Point& other) = default;

    friend ostream& operator<<(ostream& os, const Point& p);

    Point operator-=(const Point& other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }

    friend Point operator+(const Point& a, const Point& b) {
        Point tmp = a;
        tmp += b;
        return tmp;
    }

    friend vec3 operator-(const Point& a, const Point& b) {
        return {a.x - b.x, a.y - b.y, a.z - b.z};
    }

    Point operator/=(const double other) {
        x /= other;
        y /= other;
        z /= other;
        return *this;
    }

    Point operator*=(const double other){
        x *= other;
        y *= other;
        z *= other;
        return *this;
    }

    friend Point operator/(const Point& a, const double b) {
        Point tmp = a;
        tmp /= b;
        return tmp;
    }

    friend Point operator*(const Point& a, const double b){
        Point tmp = a;
        tmp *= b;
        return tmp;
    }

    explicit operator vec3() const{
        return {x, y, z};
    }

    explicit operator const array <double, 3>(){
        return array{x, y, z};
    }

    explicit operator vector<double>(){
        return vector{x, y, z};
    }
};


double area(const Point &a, const Point &b, const Point &c);


std::array<Point, 4> calcBarCoords(Point const& a, Point const& b, Point const& c);

struct Triangle{
    Point a, b, c;
    double S;
    array<Point, 4> barCoords;
    Triangle(const Point& a,const Point& b, const Point& c) : a(a), b(b), c(c){
        S = area(a,b,c);
        barCoords = calcBarCoords(a, b, c);
    }
    Triangle& operator=(const Triangle& other) = default;
};


struct iTriangle{
    int iv1, iv2, iv3;
    iTriangle(int iv1, int iv2, int iv3) : iv1(iv1), iv2(iv2), iv3(iv3){}
};

struct MarkedTriangle : Triangle{
    Point C;
    MarkedTriangle(Point const& v1, Point const& v2, Point const& v3) : Triangle(v1, v2, v3), C(v3){}
    MarkedTriangle(const Triangle &otherT): Triangle(otherT){C = otherT.c;}
};

vec3 operator+(const vec3 &lhs , const vec3 &rhs);


vec3 operator*(const vec3 &ac, const double mult);

vec3 operator-(const vec3 &lhs, const vec3 &rhs);

vec3 operator/(const vec3 &lhs , const double div);

double dot(const vec3 &a , const vec3 &b);
complex<double> dot(const vec3 &a, const array<complex<double>, 3> &b);

vec3 cross(vec3 a, vec3 b);

double norm(vec3 v);
double norm(Point p);
double normsqr(array<complex<double>, 3> p);

vec3 normal(const Point &a, const Point &b, const Point &c);

double area(const Point &a, const Point &b, const Point &c);
double area(const Triangle &t);

double dist(const Point& x, const Point& y);

using vec3c = array<complex<double>, 3>;

vec3c operator+(const vec3c& a,const vec3c& b);

vec3c operator-(const vec3c& a,const vec3c& b);

vec3c operator*(const vec3c& a, complex<double> mult);


vec3c operator*(const vec3 a, complex<double> mult);

vec3c operator+=(vec3c&a, const vec3c&other);

vec3c operator*(Point a, complex<double> mult);
#endif //SIRIUS_GEOMETRY_H
