cmake_minimum_required(VERSION 3.22)
project(Sirius)

set(CMAKE_CXX_STANDARD 20)


set(GCC_COVERAGE_COMPILE_FLAGS "-O2 -pg -larmadillo")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

add_executable(Sirius main.cpp include/Grid.h test/test.cpp include/Geometry.h src/Geometry.cpp include/Integrator.h src/Integrator.cpp src/Grid.cpp src/Integrator.tpp)

target_include_directories(Sirius PRIVATE include)