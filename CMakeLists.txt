cmake_minimum_required(VERSION 3.0.0)
project(TriangulationProject VERSION 0.1.0)

# include(CTest)
# enable_testing()

add_executable(TriangulationProject Test.cc Triangulation.cc)

# add library files
SET(Triangulation.h)

find_package(Eigen3 3.4 REQUIRED NO_MODULE)
target_link_libraries(TriangulationProject Eigen3::Eigen)

# include_directories(/usr/include/eigen3)


set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
