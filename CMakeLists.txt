cmake_minimum_required(VERSION 3.19)

PROJECT(LOWRANK)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")
find_package(Eigen3 REQUIRED)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
add_executable(${PROJECT_NAME} LRM.cpp )
