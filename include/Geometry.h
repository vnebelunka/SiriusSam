#ifndef SIRIUS_GEOMETRY_H
#define SIRIUS_GEOMETRY_H

#include <vector>
#include <cmath>
#include <array>
#include <complex>
#include <armadillo>

using namespace std;


//TODO: разделить на файл с примитивами и работой с векторами + перейти к шаблонам.

//TODO: Triangle и Marked triangle можно обьединить скорее всего


struct vec3 {
    double x, y, z;

    vec3(double x = 0, double y = 0, double z = 0): x(x), y(y), z(z){}

    vec3& operator=(const vec3& other) = default;

    friend ostream& operator<<(ostream& os, const vec3& p);

    vec3 operator+=(const vec3& other);
    vec3 operator+(const vec3& other) const{
        return {x + other.x, y + other.y, z + other.z};
    }

    vec3 operator-=(const vec3& other);
    vec3 operator-(const vec3& other) const{
        return {x - other.x, y - other.y, z - other.z};
    }

    vec3 operator/=(const double other) {
        x /= other, y /= other, z /= other;
        return *this;
    }
    vec3 operator/(const double div) const{
        return {x / div, y / div, z / div};
    }

    vec3 operator*=(const double other){
        x *= other, y *= other, z *= other;
        return *this;
    }
    vec3 operator*(const double mul) const{
        return {x * mul, y * mul, z * mul};
    }

    explicit operator array <double, 3>() const{
        return array{x, y, z};
    }

    explicit operator vector<double>() const{
        return vector{x, y, z};
    }
};


double area(const vec3 &a, const vec3 &b, const vec3 &c);
std::array<vec3, 4> calcBarCoords(vec3 const& a, vec3 const& b, vec3 const& c);

struct Triangle{
    vec3 a, b, c;
    double S;
    array<vec3, 4> barCoords;
    Triangle(const vec3& a, const vec3& b, const vec3& c) : a(a), b(b), c(c){
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
    vec3 C;
    MarkedTriangle(vec3 const& v1, vec3 const& v2, vec3 const& v3) : Triangle(v1, v2, v3), C(v3){}
    MarkedTriangle(const Triangle &otherT): Triangle(otherT){C = otherT.c;}
};


double dot(const vec3 &a , const vec3 &b);
complex<double> dot(vec3 const &a, const array<complex<double>, 3> &b);

vec3 cross(const vec3 &a, const vec3 &b);

double norm(const vec3 &v);
double normsqr(const array<complex<double>, 3> &p);

vec3 normal(const vec3 &a, const vec3 &b, const vec3 &c);

double area(const vec3 &a, const vec3 &b, const vec3 &c);
double area(const Triangle &t);

double dist(const vec3& x, const vec3& y);

using vec3c = array<complex<double>, 3>;

vec3c operator+(const vec3c& a,const vec3c& b);

vec3c operator-(const vec3c& a,const vec3c& b);

vec3c operator*(const vec3c& a, complex<double> mult);


vec3c operator*(vec3 const& a, complex<double> mult);

vec3c operator+=(vec3c&a, const vec3c&other);

#endif //SIRIUS_GEOMETRY_H
