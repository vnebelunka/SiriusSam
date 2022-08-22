#ifndef SIRIUS_GEOMETRY_H
#define SIRIUS_GEOMETRY_H

#include <vector>
#include <cmath>
#include <array>
#include <complex>

using namespace std;

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

    friend Point operator-(const Point& a, const Point& b) {
        Point tmp = a;
        tmp -= b;
        return tmp;
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

    explicit operator array<double, 3>() const{
        return array{x, y, z};
    }

    explicit operator const array <double, 3>(){
        return array{x, y, z};
    }

    explicit operator vector<double>(){
        return vector{x, y, z};
    }
};

struct Triangle{
    Point a, b, c;
    Triangle(Point a, Point b, Point c) : a(a), b(b), c(c){}
    Triangle(std::array<Point, 3> a) : a(a[0]), b(a[1]), c(a[2]){}
};

double area(const Triangle &t);

struct iTriangle{
    int iv1, iv2, iv3;
    iTriangle(int iv1, int iv2, int iv3) : iv1(iv1), iv2(iv2), iv3(iv3){}
};


//TODO: отнаследоваться
struct MarkedTriangle{
    const Triangle t;
    Point C;
    double S;
    MarkedTriangle(Point v1, Point v2, Point v3) : t(v1, v2, v3), C(v3){
        S = area(t);
    }
};

array<double, 3> operator+(const array<double, 3> &lhs , const array<double, 3> &rhs);

vector<double> operator/(const vector<double> &lhs , const double div);

vector<double>& operator*(vector<double> &ac , const double mult);

array<double, 3> operator*(const array<double, 3> &ac, const double mult);

array<double, 3> operator-(const array<double, 3> &lhs, const array<double, 3> &rhs);

double dot(const array<double, 3> &a , const array<double, 3> &b);
complex<double> dot(const array<double, 3> &a, const array<complex<double>, 3> &b);

array<double, 3> cross(std::array<double, 3> a, std::array<double, 3> b);

double norm(array<double, 3> v);
double norm(Point p);
double normsqr(array<complex<double>, 3> p);

array<double, 3> normal(const Point &a, const Point &b, const Point &c);

double area(const Point &a, const Point &b, const Point &c);

double dist(const Point& x, const Point& y);

using vec3c = array<complex<double>, 3>;

vec3c operator+(const vec3c& a,const vec3c& b);

vec3c operator-(const vec3c& a,const vec3c& b);

vec3c operator*(const vec3c& a, complex<double> mult);

using vec3 = array<double, 3>;

vec3c operator*(const vec3 a, complex<double> mult);

vec3c operator+=(vec3c&a, const vec3c&other);

vec3c operator*(Point a, complex<double> mult);
#endif //SIRIUS_GEOMETRY_H
