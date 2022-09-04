#include <iostream>
#include <complex>
#include <chrono>
#include <armadillo>
#include <sys/stat.h>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>

#include "Grid.h"
#include "Integrator.h"
//#include "matplot/matplot.h"

//TODO: сделать namespace
//TODO: документация

using namespace std::chrono;
using arma::cx_mat;
using arma::cx_vec;

complex<double> intEdge(const Grid &g, const pair<int, int>& e1, const pair<int, int>& e2, const pair<int, int>& v1,
                        const pair<int, int>& v2, double k);

void calcMatrix(const Grid &g, double k, cx_mat& M){
    size_t n = g.edges.size();
    size_t partition = n / 10;
    int i = 0, j;
    for(auto [e1, v1]: g.edges){
        j = 0;
        for(auto [e2, v2]: g.edges){
            M(i,j) = intEdge(g, e1, e2, v1, v2, k);
            ++j;
        }
        if(i % partition == 0){
            spdlog::info("{:3.0}%", double(100 * i) / n);
        }
        ++i;
    }
}

complex<double> intEdge(const Grid &g, const pair<int, int>& e1, const pair<int, int>& e2, const pair<int, int>& v1,
                       const pair<int, int>& v2, double k){
    MarkedTriangle txPlus(g.triangles.find({e1.first, e1.second, v1.first})->second);
    MarkedTriangle txMinus(g.triangles.find({e1.first, e1.second, v1.second})->second);
    MarkedTriangle tyPlus(g.triangles.find({e2.first, e2.second, v2.first})->second);
    MarkedTriangle tyMinus(g.triangles.find({e2.first, e2.second, v2.second})->second);
    complex<double> ans;
    if(g.check_dist(e1, e2)){
        ans = intFar(txPlus, tyPlus, k) + intFar(txMinus, tyMinus, k);
        ans -= intFar(txPlus, tyMinus, k) + intFar(txMinus, tyPlus, k);
    } else {
        ans = intNear(txPlus, tyPlus, k) + intNear(txMinus, tyMinus, k);
        ans -= intNear(txPlus, tyMinus, k) + intNear(txMinus, tyPlus, k);
    }
    return ans;
}

cx_vec calcF(const Grid &g, double k, vec3 Eplr, vec3 v0){
    cx_vec f(g.edges.size());
    int i = 0;
    for(auto [e, v]: g.edges){
        MarkedTriangle tPlus(g.points[e.first], g.points[e.second], g.points[v.first]);
        MarkedTriangle tMinus(g.points[e.first], g.points[e.second], g.points[v.second]);
        auto temp = intF(tPlus, k, Eplr, v0) - intF(tMinus, k, Eplr, v0);
        f[i] = temp;
        ++i;
    }
    return f;
}

int main(int argc, char* argv[]){
    spdlog::info("Welcome to spdlog!");
    auto logger = spdlog::basic_logger_mt("basic_logger", "logs/log.txt");
    char *dname = argv[1];
    double k = atof(argv[2]);
    logger->info("Input data: Filename = {}\nk = {}", dname, k);
    spdlog::info("Parsing {} file", dname);
    Grid g(dname);
    g.read();
    g.get_unique_edges();
    g.getMarkedTriangles();
    size_t n = g.edges.size();
    double d = g.diametr_grid();
    logger->info("Diameter of grid (max Edge length) = d = {}, lambda/d = {}\n", d, d / (2. * M_PI / k));
    logger->info("Num of Triangles: {}, Num of Edges: {}", g.triangles.size(), g.edges.size());
    spdlog::info("Starting calculation of Matrix coefficients");
    cx_mat A(n, n);
    calcMatrix(g, k, A);
    spdlog::info("Saving matrix");
    A.save("./logs/matrix.txt", arma::arma_ascii);
    spdlog::info("Starting calculation of right side");
    cx_vec f = calcF(g, k, {0, 1, 0}, {-1, 0, 0});
    spdlog::info("Saving right side");
    f.save("./logs/f.txt", arma::arma_ascii);
    spdlog::info("Solving linear system");
    cx_vec j;
    arma::solve(j, A, f);
    spdlog::info("Saving system solution");
    j.save("./logs/j.txt", arma::arma_ascii);
    spdlog::info("Calculating Radar cross-section");
    int parts = 360;
    vector<double> sigma(parts), x(parts);
    for(int i = 0; i < parts; ++i){
        double alpha = M_PI * i / parts;
        vec3 tau({cos(alpha), sin(alpha), 0});
        sigma[i] = calcSigma(g, j, k, tau);
        x[i] = i / 2.;
        if(i % 72 == 0) {
            spdlog::info("{:2}%", double(100 * i) / 360);
        }
    }
    spdlog::info("Saving Radar cross-section at ./logs/sigma.txt");
    ofstream out("./logs/sigma.txt");
    for(int i = 0; i < parts; ++i){
        out << x[i] << " ";
    }
    out << '\n';
    for(int i = 0; i < parts; ++i){
        out << sigma[i] << " ";
    }
}