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

    //TODO: мб лучше структурой.
    std::map<std::pair<int, int>, std::pair<int, int>> edges;
    std::map<std::tuple<int, int, int>, Triangle> triangles = {};

    Grid(const std::string &file_name);

    void read();

    void insert_edge(int e1, int e2, int v);

    void get_unique_edges();

    std::array<double, 3> get_point_coord(int i);

    std::array<vec3, 3> get_points(iTriangle t);

    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> get_coordinates();

    void print_edges();

    //void plot_grid();

    std::vector<Triangle> get_unique_tringles();

    void getMarkedTriangles(){
        for(auto [e, v] : edges){
            triangles.insert(std::make_pair<std::tuple<int, int, int>, Triangle>
                ({e.first, e.second, v.first}, {points[e.first], points[e.second], points[v.first]})
                );
            triangles.insert(std::make_pair<std::tuple<int, int, int>, Triangle>
                                     ({e.first, e.second, v.second}, {points[e.first], points[e.second], points[v.second]})
            );
        }
    }

    ~Grid(){data_file.close();}

    double max_edge(const iTriangle &t) const{
        double d1 = dist(points[t.iv1], points[t.iv2]);
        double d2 = dist(points[t.iv2], points[t.iv3]);
        double d3 = dist(points[t.iv1], points[t.iv3]);
        return max(d1, max(d2, d3));
    }


    //TODO: лучше проверять это через диаметры.
    //TODO: переписать на функцию от 4 треугольников: они уже найдены в calcMatrix
    bool check_dist(pair<int, int> p1, pair<int, int> p2) const{
        vec3 m1 = (points[p1.first] + points[p1.second]) / 2;
        vec3 m2 = (points[p2.first] + points[p2.second]) / 2;
        auto v1 = edges.find(p1)->second;
        auto v2 =  edges.find(p2)->second;
        double d1 = max_edge({p1.first, p1.second, v1.first});
        double d2 = max_edge({p1.first, p1.second, v1.second});
        double d3 = max_edge({p2.first, p2.second, v2.first});
        double d4 = max_edge({p2.first, p2.second, v2.second});
        double d = max(max(d1, d2), max(d3, d4));
        return dist(m1, m2) > 2 * d;
    }

    double diametr_grid(){
        double ans = -1;
        for(auto &t : itriangles){
            double temp = max_edge(t);
            if(temp > ans) {
                ans = temp;
                cout << t.iv1<<" " << t.iv2<<" " << t.iv3 << '\n';
            }
        }
        return ans;
    }
};

//TODO: сделать как header
extern void test_unique_edges();
extern void test_integral_1d();
extern void test_equiv_integral();
#endif //SIRIUS_GRID_H