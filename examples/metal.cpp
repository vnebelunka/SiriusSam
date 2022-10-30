#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>

#include "Grid.h"
#include "Integrator.h"
#include "operatorK.h"
#include "operatorR.h"
#include "progressbar.h"

void calcFE(const Grid &g, double k, vec3 Eplr, vec3 v0, cx_vec &f) {
    int i = 0;
    progressbar p(g.edges_inner_enum.size());
    for(auto [e, v]: g.edges){
        if(v.second == -1){
            continue;
        }
        f[i] = -intEdge_en_Einc(g, e, v, k, Eplr, v0);
        ++i;
        p.update();
    }
    std::cerr<<std::endl;
}

void calcMatrixE(const Grid &g, double k, cx_mat &M) {
    int i = 0, j;
    progressbar p(int(g.edges_inner_enum.size()));
    for(auto [e1, v1]: g.edges){
        j = 0;
        if(v1.second == -1){
            continue;
        }
        for(auto [e2, v2]: g.edges){
            if(v2.second == -1){
                continue;
            }
            M(i,j) = intEdge_en_Re(g, e1, e2, v1, v2, k) - 0.5 * intEdge_e1_e2(g, e1, e2, v1, v2);
            ++j;
        }
        ++i;
        p.update();
    }
    std::cerr << std::endl;
}

void calcFM(const Grid &g, double k, vec3 Eplr, vec3 v0, cx_vec& f) {
    progressbar p(g.edges_inner_enum.size());
    int i = 0;
    for(auto [e, v]: g.edges){
        f[i] = intEdge_e_Einc(g,e,v,k,Eplr,v0);
        ++i;
        p.update();
    }
    std::cerr << std::endl;
}

void calcMatrixM(const Grid &g, double k, cx_mat &M) {
    int i = 0, j;
    progressbar p(g.edges.size());
    for(auto [e1, v1]: g.edges){
        j = 0;
        if(v1.second == -1){
            continue;
        }
        for(auto [e2, v2]: g.edges){
            if(v2.second == -1){
                continue;
            }
            M(i,j) = intEdge_e_Ke(g, e1, e2, v1, v2, k);
            ++j;
        }
        ++i;
        p.update();
    }
    std::cerr << std::endl;
}

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