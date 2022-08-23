#include "../include/Grid.h"
#include "../include/Integrator.h"
#include <iostream>

static bool find_triangle_by_edge(Grid &r, int e1, int e2, int v){
    if(r.edges.find({e1, e2}) == r.edges.end()){
        std::cout << "no edge: " << e1 << " " << e2<< "!\n";
        return false;
    } else {
        auto v_edge = r.edges[{e1, e2}];
        if(v_edge.first != v && v_edge.second != v){
            std::cout << "no triangle" << e1 << " " << e2 << " " << v << '\n';
            return false;
        }
    }
    return true;
}

void test_unique_edges(){
    Grid r("./data/sphere_5cm_reg.obj");
    r.read();
    r.get_unique_edges();
    std::set<std::array<int, 3>> unique_triangles;
    for(auto &t: r.itriangles){
        std::array<int, 3> temp({t.iv1, t.iv2, t.iv3});
        std::sort(temp.begin(), temp.end());
        unique_triangles.insert(temp);
    }
    /*
     * Все треугольники "ребро - вершина существуют"
    */
    for(auto [key, value]: r.edges){
        std::array<int, 3> t({key.first, key.second, value.first});
        std::sort(t.begin(), t.end());
        if(unique_triangles.find(t) == unique_triangles.end()){
            std::cout << "tringle: " << t[0] << " " << t[1] << " " << t[2] << "not in grid!\n";
            return;
        }
        if(value.second != -1){
            t = {key.first, key.second, value.second};
            std::sort(t.begin(), t.end());
            if(unique_triangles.find(t) == unique_triangles.end()){
                std::cout << "tringle: " << t[0] << " " << t[1] << " " << t[2] << "not in grid!\n";
                return;
            }
        }
    }
    /*
     * Ребро каждого треугольника учтено.
     */
    for(auto &t : unique_triangles){
        find_triangle_by_edge(r, t[0], t[1], t[2]);
        find_triangle_by_edge(r, t[0], t[2], t[1]);
        find_triangle_by_edge(r, t[1], t[2], t[0]);
    }
}

void test_integral_1d(){
    cout << "Testing func for Surface 1d Integral with Gaus-4 method\n";
    cout << "-------------------------------------------------------\n";
    Grid r("../data/sphere_5cm_reg.obj");
    r.read();
    r.get_unique_edges();
    auto triangles = r.get_unique_tringles();
    auto int1 = surfIntegral1D(triangles, [](Point p) -> complex<double> { return {1., 0}; });
    cout << "Test1: f == 1: " << int1 << " as res; " << 4 * M_PI * 0.5 * 0.5 << " as true value.\n";
    auto int2 = surfIntegral1D(triangles, [](Point p) -> complex<double> { return norm(p); });
    cout << "Test2: f == r: " << int2 << " as res; " << 4 * M_PI * 0.5 * 0.5 * 0.5 << " as true value.\n";
    auto int3 = surfIntegral1D(triangles, [](Point p) -> complex<double> { return p.x + p.y + abs(p.z); });
    cout << "Test3: f = x + y + |z|: " << int3 << " as res; " << 2 * M_PI * 0.5 * 0.5 * 0.5 << " as true value.\n";
    auto int4 = surfIntegral1D(triangles,
                               [](Point p) -> complex<double> {
                                   Point x(1, 0, 0);
                                   return 1. / dist(p, x);
                               });
    cout << "Test4: f = 1/|a-x|, a = (100, 0, 0): " << int4 << " as res; " << 4 * M_PI * 0.5 * 0.5  << " as true value.\n";
    auto int5 = surfIntegral1D(triangles,
                               [](Point p) -> complex<double> {
                                   Point x(0.1, 0.1, 0.1);
                                   return 1. / dist(p, x);
                               });
    cout << "Test5: f = 1/|a-x|, a = (10, 10, 10): " << int5 << " as res; " << 4 * M_PI * 50 << " as true value.\n";
    cout << "Test ended\n";
    cout << "----------" << endl;
}


void test_equiv_integral(){
    Grid r("./data/sphere_5cm_reg.obj");
    r.read();
    r.get_unique_edges();
    auto triangles = r.get_unique_tringles();
    double maxdif = 0;
    int k = 0;
    for(auto t : triangles){
        auto res1 = integrateGaus(t, [](Point p) -> complex<double> {
            Point x(100, 0, 0);
            return 1. / dist(p, x);
        });
        Point x(100, 0, 0);
        auto res2 = integral1Divr(t, x);
        auto temp = abs(res1 - res2) / max(res1.real(), res2);
        if(temp > maxdif){
            maxdif = temp;
        }
        ++k;
    }
    cout << maxdif;
}