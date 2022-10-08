/**
 * \file
 * Header file for all Point and vector operations
 */

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


/**
 * \brief Struct for 3d point
 */
struct vec3 {


    double x, y, z;
    vec3(double x = 0, double y = 0, double z = 0): x(x), y(y), z(z){}

    vec3& operator=(const vec3& other) = default;

    friend ostream& operator<<(ostream& os, const vec3& p);

    bool operator ==(const vec3& other) const{return x == other.x && y == other.y && z == other.z;}
    bool operator <(const vec3& other) const{
        if(x < other.x) {return true;}
        if(x > other.x) {return false;}
        if(y < other.y) {return true;}
        if(y > other.y) {return false;}
        if(z < other.z) {return true;}
        if(z > other.z) {return false;}
        return false;
    }
    /*
     * Arithmetic operations
     */
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

    /*
     * Convert constructors
     */
    explicit operator array <double, 3>() const{
        return array{x, y, z};
    }

    explicit operator vector<double>() const{
        return vector{x, y, z};
    }
};


/**
 * Area of triangle with (a, b, c) vertices.
 * @param a 1st point of the triangle
 * @param b 2nd point of the triangle
 * @param c 3rd point of the triangle
 * @return area of triangle
 */
double area(const vec3 &a, const vec3 &b, const vec3 &c);

/**
 * Calculates the barycentric coordinates of triangle (a, b, c)
 * @param a 1st point of the triangle
 * @param b 2nd point of the triangle
 * @param c 3rd point of the
 * @return array of 4 points (barycentric coordinates).
 * \TODO: make as template for num of points
 */
std::array<vec3, 4> calcBarCoords(vec3 const& a, vec3 const& b, vec3 const& c);



/**
 * Dot product of two vectors \f$(a,b) = a_1 b_1 + a_2 b_2 + a_3 b_3\f$
 * @param[in] a first vector
 * @param[in] b second vector
 * @return dot product of two vectors
 */
double dot(const vec3 &a , const vec3 &b);


/**
 * \brief Triangle struct with precalced area and Barycentric coordinates.
 */
struct Triangle{
    vec3 a, b, c; ///< Points of triangle
    array<vec3, 3> sorted;
    double S; ///< area of Triangle
    array<vec3, 4> barCoords; ///< Barycentric coordinates
    Triangle(const vec3& a, const vec3& b, const vec3& c): a(a), b(b), c(c){
        S = area(a,b,c);
        sorted[0] = a, sorted[1] = b, sorted[2] = c;
        sort(sorted.begin(), sorted.end());
        barCoords = calcBarCoords(a, b, c);
    }
    Triangle& operator=(const Triangle& other) = default;
    bool operator==(Triangle const &other) const{
        return sorted[0] == other.sorted[0] && sorted[1] == other.sorted[1] && sorted[2] == other.sorted[2];
    }
};

/**
 * \brief Triangle in grid as 3 indexes in point array
 */
struct iTriangle{
    int iv1, iv2, iv3;
    iTriangle(int iv1, int iv2, int iv3) : iv1(iv1), iv2(iv2), iv3(iv3){}
};

/**
 * \brief Triangle with 1 marked point for RWG function
 */
struct MarkedTriangle : Triangle{
    vec3 marked; ///< Marked point of triangle
    MarkedTriangle(vec3 const& v1, vec3 const& v2, vec3 const& v3) : marked(v3), Triangle(v1, v2, v3){}
    MarkedTriangle(const Triangle &otherT): Triangle(otherT){ marked = otherT.c;}
    bool operator == (const MarkedTriangle &otherT) const{
        return sorted[0] == otherT.sorted[0] && sorted[1] == otherT.sorted[1] && sorted[2] == otherT.sorted[2];
    } //TODO: можно лучше.
};


/**
 * complex dot of 2 vectors (2nd is vector with complex coords) \f$ (a,b) = a_1 \overline b_1 + a_2 \overline b_2 + a_3 \overline b_3 \f$
 * @param a real vector
 * @param b complex vector
 * @return complex dot \f$ (a,b) = a_1 \overline b_1 + a_2 \overline b_2 + a_3 \overline b_3 \f$
 */
complex<double> dot(vec3 const &a, const array<complex<double>, 3> &b);

/**
 *  cross product of two vectors.
 * @param a first vector
 * @param b second vector
 * @return cross product of two vectors
 */
vec3 cross(const vec3 &a, const vec3 &b);

/**
 * return 2nd norm of vector \f$|v| = \sqrt{(v, v)} \f$
 * @param[in] v vector
 * @return 2nd norm of vector
 */
double norm(const vec3 &v);

/**
 * squqred norm of complex vector
 * @param p vector
 * @return complex norm of p
 */
double normsqr(const array<complex<double>, 3> &p);

/**
 * Normal vector of plane on points a, b, c
 * @param a 1st point
 * @param b 2nd point
 * @param c 3rd point
 * @return normal vector \f$ n \perp abc \f$
 */
vec3 normal(const vec3 &a, const vec3 &b, const vec3 &c);

vec3 normal(const Triangle &t);


/**
 * area of triangle
 * @param t triangle
 * @return area of t
 */
double area(const Triangle &t);

/**
 * euclidean distance between two points
 * @param x 1st point
 * @param y 2nd point
 * @return euclidean distance \f$ \rho(x, y) \f$
 */
double dist(const vec3& x, const vec3& y);


/*
 * Operators for complex 3d vectors
 */
using vec3c = array<complex<double>, 3>;

vec3c operator+(const vec3c& a,const vec3c& b);

vec3c operator-(const vec3c& a,const vec3c& b);

vec3c operator*(const vec3c& a, complex<double> mult);


vec3c operator*(vec3 const& a, complex<double> mult);

vec3c operator+=(vec3c&a, const vec3c&other);


complex<double> cdot(vec3c const &a, vec3 const &b);

vec3c cross(vec3c const&a, vec3 const &b);

#endif //SIRIUS_GEOMETRY_H
