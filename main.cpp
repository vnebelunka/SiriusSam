#include <iostream>
#include <complex>
#include <armadillo>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>

#include "Grid.h"
#include "Integrator.h"
#include "Magnetic.h"
#include "Electric.h"
#include "progressbar.h"


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
    g.enum_inner_edges();
    size_t n = g.edges_inner_enum.size();
    double d = g.diametr_grid();
    logger->info("Diameter of grid (max Edge length) = d = {}, lambda/d = {}\n", d, d / (2. * M_PI / k));
    logger->info("Num of Triangles: {}, Num of Edges: {}", g.triangles.size(), g.edges.size());
    spdlog::info("Starting calculation of Matrix coefficients");
    /*
     * Matrix calculation
     */
    cx_mat A(n, n);
    calcMatrixE(g, k, A);
    spdlog::info("Saving matrix");
    A.save("./logs/matrix.txt", arma::arma_ascii);
    /*
     * Calcing free vector f
     */
    spdlog::info("Starting calculation of right side");
    cx_vec f(n);
    calcFE(g, k, {0, 1, 0}, {-1, 0, 0}, f);
    spdlog::info("Saving right side");
    f.save("./logs/f.txt", arma::arma_ascii);
    /*
     * Solving System
     */
    spdlog::info("Solving linear system");
    cx_vec j;
    arma::solve(j, A, f);
    spdlog::info("Saving system solution");
    j.save("./logs/j.txt", arma::arma_ascii);
    /*
     * Calculating radar cross-section
     */
    spdlog::info("Calculating Radar cross-section");
    progressbar p(360);
    int parts = 360;
    vector<double> sigma(parts), x(parts);
    for(int i = 0; i < parts; ++i){
        double alpha = M_PI * i / parts;
        vec3 tau({cos(alpha), sin(alpha), 0});
        sigma[i] = calcSigmaE(g, j, k, tau);
        x[i] = i / 2.;
        p.update();
    }
    std::cerr << std::endl;
    spdlog::info("Saving Radar cross-section at ./logs/sigma.txt");
    ofstream out("./logs/sigma.txt");
    for(int i = 0; i < parts; ++i){
        out << x[i] << " ";
    }
    out << '\n';
    for(int i = 0; i < parts; ++i){
        out << sigma[i] << " ";
    }
    /*
     * Calculating vector flow
     */
    calcTotalFlow(g,j,"./logs/points.txt", "./logs/greal.txt", "./logs/gimag.txt", "./logs/norms.txt");
}