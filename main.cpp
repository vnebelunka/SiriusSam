#include <iostream>
#include <complex>
#include <armadillo>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>

#include "Grid.h"
#include "Integrator.h"
#include "Magnetic.h"
#include "Electric.h"
//#include "matplot/matplot.h"

//TODO: сделать namespace
//TODO: документация


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
    calcMatrixE(g, k, A);
    spdlog::info("Saving matrix");
    A.save("./logs/matrix.txt", arma::arma_ascii);
    spdlog::info("Starting calculation of right side");
    cx_vec f = calcFE(g, k, {0, 1, 0}, {-1, 0, 0});
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
    calcTotalFlow(g,j,"./logs/points.txt", "./logs/greal.txt", "./logs/gimag.txt");
}