cmake_minimum_required(VERSION 3.15)
project(cmd_rtree_primer)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add your source files here
add_executable(cmd_rtree_primer 
    cmd_rtree_primer.cpp)

# Find Boost (adjust components as needed)
find_package(Boost REQUIRED COMPONENTS geometry)

target_link_libraries(cmd_rtree_primer PRIVATE Boost::geometry)
