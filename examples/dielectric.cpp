#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>

#include "Grid.h"
#include "Integrator.h"
#include "operatorK.h"
#include "operatorR.h"

void calcMatrixD(const Grid &g, double k1, double k2, double eps1, double eps2, cx_mat &M){

}
void calcFD(const Grid &g, double k, vec3 Eplr, vec3 v0, cx_vec& f) {
    int i = 0;
    for(auto [e, v]: g.edges){
        //f[i] = intEdge_e_Einc(g,e,v,k,Eplr,v0);
        ++i;
    }
}

void dielecrtic_proceed(shared_ptr<spdlog::logger> &logger, Grid &g, double k1, double k2, double eps1, double eps2){
    /*
     * Matrix calculation
     */
    size_t n = g.edges_inner_enum.size();
    cx_mat A(n * 2, n * 2);
    calcMatrixD(g,k1,k2,eps1,eps2, A);
    spdlog::info("Saving matrix");
    A.save("./logs/matrix.txt", arma::arma_ascii);
    /*
     * Calcing free vector f
     */
    spdlog::info("Starting calculation of right side");
    cx_vec f(n * 2);
    calcFD(g,k2, {0,1,0}, {-1, 0, 0}, f);
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
    /*int parts = 360;
    vector<double> sigma(parts), x(parts);
    for(int i = 0; i < parts; ++i){
        double alpha = M_PI * i / parts;
        vec3 tau({cos(alpha), sin(alpha), 0});
        sigma[i] = calcSigmaD(g, j, k1, k2, tau);
        x[i] = i / 2.;
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
     */
    /*
     * Calculating vector flow
     */
    //calcTotalFlow(g,j,"./logs/points.txt", "./logs/greal.txt", "./logs/gimag.txt", "./logs/norms.txt");*/
}

