cmake_minimum_required(VERSION 3.19)

PROJECT(LOWRANK)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")
enable_testing()
find_package(GTEST REQUIRED)
find_package(Eigen3 REQUIRED)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include)

add_executable(test_FDA test/test_FDA.cpp)
target_link_libraries(test_FDA ${GTEST_LIBRARIES})
add_test(NAME test_FDA COMMAND test_FDA)

# add_executable(FDA check.cpp )
