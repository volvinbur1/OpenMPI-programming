cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(lab1)

find_package(MPI REQUIRED)

add_executable(lab1 lab1/main.cxx)

target_link_libraries(lab1 PRIVATE MPI::MPI_CXX)