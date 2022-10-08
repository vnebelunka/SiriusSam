#ifndef SIRIUS_GRID_H
#define SIRIUS_GRID_H

#include <string>
#include <vector>
#include <set>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <algorithm>
#include <map>
#include "Geometry.h"
#include <complex>


struct Grid{
    std::ifstream data_file;
    size_t num_points = 0;
    std::vector<vec3> points = {};
    size_t num_triangles = 0;
    std::vector<iTriangle> itriangles = {};

    std::map<std::pair<int, int>, std::pair<int, int>> edges;
    std::map<std::tuple<int, int, int>, MarkedTriangle> triangles = {};
    std::map<std::pair<int,int>, int> edges_inner_enum;

    Grid(const char *file_name);

    void read();

    void insert_edge(int e1, int e2, int v);

    void get_unique_edges();

    std::array<double, 3> get_point_coord(int i);

    std::array<vec3, 3> get_points(iTriangle t);

    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> get_coordinates();

    void print_edges();

    //void plot_grid();

    std::vector<Triangle> get_unique_tringles();

    void getMarkedTriangles();

    ~Grid(){data_file.close();}

    double max_edge(const iTriangle &t) const{
        double d1 = dist(points[t.iv1], points[t.iv2]);
        double d2 = dist(points[t.iv2], points[t.iv3]);
        double d3 = dist(points[t.iv1], points[t.iv3]);
        return max(d1, max(d2, d3));
    }

    void enum_inner_edges(){
        int i = 0;
        for(auto [e, v]: edges){
            if(v.second != -1){
                edges_inner_enum[e] = i++;
            }
        }
    }

    //TODO: переписать на функцию от 4 треугольников: они уже найдены в calcMatrixM
    bool check_dist(pair<int, int> p1, pair<int, int> p2) const;

    double diametr_grid();
};

//TODO: сделать как header
extern void test_unique_edges();
extern void test_integral_1d();
extern void test_equiv_integral();
#endif //SIRIUS_GRID_H