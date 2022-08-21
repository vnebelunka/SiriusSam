#include <iostream>
#include <complex>
#include <chrono>

#include "include/Grid.h"
#include "include/Integrator.h"
//#include "matplot/matplot.h"

//TODO: сделать namespace
//TODO: документация

using namespace std::chrono;

complex<double> intEdge(const Grid &g, const pair<int, int>& e1, const pair<int, int>& e2, const pair<int, int>& v1,
                        const pair<int, int>& v2, double k);

vector<vector<complex<double>>> calcMatrix(const Grid &g, double k){
    size_t n = g.edges.size();
    vector<vector<complex<double>>> M(n, vector<complex<double>>(n));
    std::ofstream outM("../matrix10_512.txt");
    int i = 0, j;
    for(auto [e1, v1]: g.edges){
        j = 0;
        for(auto [e2, v2]: g.edges){
            M[i][j] = intEdge(g, e1, e2, v1, v2, k);
            outM << M[i][j] << " ";
            ++j;
        }
        outM << "\n";
        ++i;
        cout << i << '/' << n << '\n';
    }
    return M;
}

complex<double> intEdge(const Grid &g, const pair<int, int>& e1, const pair<int, int>& e2, const pair<int, int>& v1,
                       const pair<int, int>& v2, double k){
    MarkedTriangle txPlus(g.points[e1.first], g.points[e1.second], g.points[v1.first]);
    MarkedTriangle txMinus(g.points[e1.first], g.points[e1.second], g.points[v1.second]);
    MarkedTriangle tyPlus(g.points[e2.first], g.points[e2.second], g.points[v2.first]);
    MarkedTriangle tyMinus(g.points[e2.first], g.points[e2.second], g.points[v2.second]);
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

vector<complex<double>> calcF(const Grid &g, double k, vec3 Eplr, vec3 v0){
    vector<complex<double>> f(g.edges.size());
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


int main(){
    Grid g("../data/sphere_5cm_reg.obj");
    g.read();
    g.get_unique_edges();
    int n = g.edges.size();
    double k = 10;
    string s;
    cin >> s;
    if(s == "calc") {
        auto A = calcMatrix(g, k);
        auto f = calcF(g, k, {0, 1, 0}, {-1, 0, 0});
        cout << "f ready\n";
        ofstream outF("../f10_512.txt");
        for (int i = 0; i < n; ++i) {
            outF << f[i] << '\n';
        }
    } else if(s == "plot"){
        FILE *fd = fopen("../results/k=10/j10_3012.txt", "r");
        vector<complex<double>> j(n);
        for(int i = 0; i < n; ++i){
            double a, b;
            fscanf(fd, " (%lf", &a);
            fscanf(fd, "%lfj)", &b);
            //cout << a <<" "<< b << '\n';
            j[i] = complex<double>(a, b);
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
        ofstream out("res10_3012.txt");
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