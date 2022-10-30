cmake_minimum_required(VERSION 3.22)
project(Sirius)

set(CMAKE_CXX_STANDARD 20)

if(NOT TARGET spdlog)
    # Stand-alone build
    find_package(spdlog REQUIRED)
endif()


find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

set(GCC_COVERAGE_COMPILE_FLAGS "-O2 -larmadillo")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

add_executable(Sirius main.cpp include/Grid.h test/test.cpp include/Geometry.h
        src/Geometry.cpp include/Integrator.h src/Integrator.cpp src/Grid.cpp
        src/Integrator.tpp include/operatorK.h src/operators/operatorK.cpp
        src/operators/operatorR.cpp include/operatorR.h examples/metal.cpp
        examples/dielectric.cpp)

target_include_directories(Sirius PRIVATE include)
target_link_libraries(Sirius PRIVATE spdlog::spdlog ${ARMADILLO_LIBRARIES})