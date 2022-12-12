#include <iostream>
#include <complex>
#include <armadillo>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>

#include "Grid.h"

const double eps0 = 8.85418782e-12;
const double mu0 = 1.25663706e10-6;
const double c0 = 299792458;


void metal_proceed(shared_ptr<spdlog::logger> &logger, Grid &g, double k, char mod='M', bool load=false);
void dielecrtic_proceed(shared_ptr<spdlog::logger> &logger, Grid &g, double w, double k1, double k2, double eps1, double eps2,
                        bool load=false);

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
    double d = g.diametr_grid();
    logger->info("Diameter of grid (max Edge length) = d = {}, lambda/d = {}\n", d, d / (2. * M_PI / k));
    logger->info("Num of Triangles: {}, Num of Edges: {}", g.triangles.size(), g.edges.size());

    double eps1 = 1.;
    double eps2 = 1.;
    double mu1 = 1.;
    double  mu2 = 1.;

    double w = 30;

    double v = w / (2 * M_PI);

    double n1 = sqrt(eps1 * eps0 * mu1 * mu0);
    double n2 = sqrt(eps2 * eps0 * mu2 * mu0);


    double k1 = n1 * w;
    double k2 = n2 * w;


    //dielecrtic_proceed(logger,g, w, k1, k2, eps1, eps2, false);

    metal_proceed(logger,g, k1, 'E', false);
}