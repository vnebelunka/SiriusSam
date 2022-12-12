#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>

#include "Grid.h"
#include "Integrator.h"
#include "operatorK.h"
#include "operatorR.h"
#include "progressbar.h"

void calcMatrixD(const Grid &g, double w, double k1, double k2, double eps1, double eps2, cx_mat &M){
    int i = 0;
    int n = g.edges_inner_enum.size();
    progressbar p(int(g.edges_inner_enum.size()));
    for(auto [e1, v1]: g.edges){
        int k = 0;
        for(auto [e2, v2]: g.edges){
//            // equations in 1 region
//            auto temp = intEdge_e_Ke(g, e1, e2, v1, v2, k1) / (eps1 * w);
//            M(i, k) = temp;
//            temp = -intEdge_e_Ren(g, e1, e2, v1, v2, k1);
//            M(i, k + n) = temp + 0.5 * intEdge_e1_e2(g,e1,e2,v1,v2);
//            // equations in 2nd region
//            temp = intEdge_e_Ke(g,e1,e2,v1,v2, k2) / (w * eps2);
//            M(i + n, k) = temp;
//            temp = -intEdge_e_Ren(g, e1, e2, v1, v2, k2);
//            M(i + n, k + n) = temp - 0.5 * intEdge_e1_e2(g,e1,e2,v1,v2);
//            ++k;
            assert(i < g.edges.size() && k < g.edges.size());
            auto temp = 1 / (w * eps1) * intEdge_e_Ke(g,e1,e2,v1,v2,k1) - 1 / (w * eps2) * intEdge_e_Ke(g,e1,e2,v1,v2,k2);
            M(i,k) = temp;
            temp = 1 / (w * eps1) * intEdge_e_Ke(g,e1,e2,v1,v2,k1) + 1 / (w * eps2) * intEdge_e_Ke(g,e1,e2,v1,v2,k2);
            M(i+n,k) = temp;
            temp = intEdge_e_Ren(g,e1,e2,v1,v2,k2) - intEdge_e_Ren(g,e1,e2,v1,v2,k1) + intEdge_e1_e2(g,e1,e2,v1,v2);
            M(i,n+k) = temp;
            temp = -intEdge_e_Ren(g,e1,e2,v1,v2,k2) -intEdge_e_Ren(g,e1,e2,v1,v2,k1);
            M(i+n,n+k) = temp;
            ++k;
        }
        ++i;
        p.update();
    }
}
void calcFD(const Grid &g, double k, vec3 Eplr, vec3 v0, cx_vec& f) {
    int i = 0;
    int n = g.edges_inner_enum.size();
    for(auto [e, v]: g.edges){
        f[i] = -intEdge_e_Einc(g,e,v,k,Eplr,v0);
        f[i + n] = 0;
        ++i;
    }
}

void dielecrtic_proceed(shared_ptr<spdlog::logger> &logger, Grid &g, double w, double k1, double k2, double eps1, double eps2,
                        bool load = false){
    size_t n = g.edges_inner_enum.size();
    cx_vec j(2 * n);
    if(!load) {
        spdlog::info("Starting calculation of Matrix coefficients");
        /*
         * Matrix calculation
         */
        cx_mat A(n * 2, n * 2);
        calcMatrixD(g, w, k1, k2, eps1, eps2, A);
        spdlog::info("Saving matrix");
        A.save("./logs/matrix.txt", arma::arma_ascii);
        /*
         * Calcing free vector f
         */
        spdlog::info("Starting calculation of right side");
        cx_vec f(n * 2);
        calcFD(g, k2, {0, 1, 0}, {-1, 0, 0}, f);
        spdlog::info("Saving right side");
        f.save("./logs/f.txt", arma::arma_ascii);
        /*
         * Solving System
         */
        spdlog::info("Solving linear system");
        arma::solve(j, A, f);
        spdlog::info(   "Saving system solution");
        j.save("./logs/j.txt", arma::arma_ascii);
    } else {
        j.load("./logs/j.txt", arma::arma_ascii);
    }

    cx_vec je(n), jm(n);
    for(int i = 0; i < n; ++i){
        je(i) = j(i);
        jm(i) = j(i + n);
    }

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
        sigma[i] = calcSigmaEM(g, je, jm, k2, w,eps2, tau);
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

    calcTotalFlow(g,je,"./logs/points.txt", "./logs/grealE.txt", "./logs/gimagE.txt", "./logs/norms.txt", 'e');
    calcTotalFlow(g,jm,"./logs/points.txt", "./logs/grealM.txt", "./logs/gimagM.txt", "./logs/norms.txt", 'n');
}

