cmake_minimum_required(VERSION 3.24)
project(traveling_salesperson)

set(CMAKE_CXX_STANDARD 17)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

add_executable(traveling_salesperson main.cpp tsp.cpp)
target_link_libraries(traveling_salesperson PUBLIC MPI::MPI_CXX)
target_link_libraries(traveling_salesperson PRIVATE OpenMP::OpenMP_CXX)
