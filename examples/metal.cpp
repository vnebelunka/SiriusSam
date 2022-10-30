#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>

#include "Grid.h"
#include "Integrator.h"
#include "Magnetic.h"
#include "Electric.h"
#include "progressbar.h"


void metal_proceed(shared_ptr<spdlog::logger> &logger, Grid &g, double k, char mod='M'){
    /*
     * Matrix calculation
     */
    size_t n = g.edges_inner_enum.size();
    cx_mat A(n, n);
    if(mod == 'E') {
        calcMatrixE(g, k, A);
    } else {
        calcMatrixM(g, k, A);
    }
    spdlog::info("Saving matrix");
    A.save("./logs/matrix.txt", arma::arma_ascii);
    /*
     * Calcing free vector f
     */
    spdlog::info("Starting calculation of right side");
    cx_vec f(n);
    if(mod == 'E') {
        calcFE(g, k, {0, 1, 0}, {-1, 0, 0}, f);
    } else {
        calcFM(g, k, {0, 1, 0}, {-1, 0, 0}, f);
    }
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
        if(mod == 'E') {
            sigma[i] = calcSigmaE(g, j, k, tau);
        } else {
            sigma[i] = calcSigmaM(g, j, k, tau);
        }
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