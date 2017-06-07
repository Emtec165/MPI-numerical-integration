cmake_minimum_required(VERSION 3.6)
project(ar_proj_3)

set(CMAKE_CXX_STANDARD 11)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set(SOURCE_FILES main.cpp)
add_executable(ar_proj_2 ${SOURCE_FILES})
target_link_libraries(ar_proj_2 ${MPI_LIBRARIES})