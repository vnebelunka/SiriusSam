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

void metal_proceed(shared_ptr<spdlog::logger> &logger, Grid &g, double k, char mod='M');

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
    spdlog::info("Starting calculation of Matrix coefficients");
    metal_proceed(logger, g, k, 'E');
}