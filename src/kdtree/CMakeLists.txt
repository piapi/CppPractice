project(kdtree)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "-std=c++14 -O3")

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)



find_package(PCL REQUIRED QUIET)


# Eigen
include_directories("/usr/include/eigen3")

include_directories(${PCL_INCLUDE_DIRS})
add_executable(kdt KD_Tree.cpp)
