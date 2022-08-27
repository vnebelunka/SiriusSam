#include <iostream>
#include <complex>
#include <chrono>
#include <armadillo>

#include "include/Grid.h"
#include "include/Integrator.h"
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
    int i = 0, j;
    for(auto [e1, v1]: g.edges){
        j = 0;
        for(auto [e2, v2]: g.edges){
            M(i,j) = intEdge(g, e1, e2, v1, v2, k);
            ++j;
        }
        ++i;
        cout << i << '/' << n << '\n';
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
    char *dname = argv[1];
    double k = atof(argv[2]);
    Grid g(dname);
    g.read();
    g.get_unique_edges();
    g.getMarkedTriangles();
    size_t n = g.edges.size();
    string s;
    cin >> s;
    if(s == "calc") {
        ofstream outA("./calcs/matrix.txt");
        cx_mat A(n, n);
        calcMatrix(g, k, A);
        for(int i = 0; i < n; ++i){
            for(int j = 0; j < n; ++j){
                outA << A(i, j) << " ";
            }
            outA << '\n';
        }
        cx_vec f = calcF(g, k, {0, 1, 0}, {-1, 0, 0});
        cout << "f ready" << endl;
        ofstream outF("./calcs/f.txt");
        for (int i = 0; i < n; ++i) {
            outF << f(i) << '\n';
        }
        cx_vec j = arma::solve(A, f);
        ofstream outJ("./calcs/j.txt");
        for(int i = 0; i < n; ++i){
            outJ <<j(i) << '\n';
        }
        int parts = 360;
        vector<double> sigma(parts), x(parts);
        for(int i = 0; i < parts; ++i){
            double alpha = M_PI * i / parts;
            vec3 tau({cos(alpha), sin(alpha), 0});
            sigma[i] = calcSigma(g, j, k, tau);
            x[i] = i / 2.;
            cout << i  << '/' << parts << '\n';
        }
        ofstream out("./calcs/sigma.txt");
        for(int i = 0; i < parts; ++i){
            out << x[i] << " ";
        }
        out << '\n';
        for(int i = 0; i < parts; ++i){
            out << sigma[i] << " ";
        }
    } else {
        double d = g.diametr_grid();
        cout << "diameter of grid = " << d << '\n';
        cout << "lambda =" << 2. * M_PI / k <<  '\n';
        cout << "d/lambda = " << d / (2. * M_PI / k) << '\n';
    }
}