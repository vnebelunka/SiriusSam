cmake_minimum_required(VERSION 3.22)
project(Sirius)

set(CMAKE_CXX_STANDARD 20)

#add_subdirectory(matplotplusplus)

add_executable(Sirius main.cpp include/Grid.h test/test.cpp include/Geometry.h src/Geometry.cpp include/Integrator.h src/Integrator.cpp src/Grid.cpp)

#target_link_libraries(Sirius PUBLIC matplot)

#target_link_libraries(Sirius ${Eigen_LIBRARIES})