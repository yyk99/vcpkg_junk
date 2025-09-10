cmake_minimum_required(VERSION 3.15)
project(cmd_boost_rtree_traversal)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add your source files here
add_executable(cmd_boost_rtree_traversal
    cmd_boost_rtree_traversal.cpp)

# Find Boost (adjust components as needed)
find_package(Boost REQUIRED geometry)

target_link_libraries(cmd_boost_rtree_traversal PRIVATE Boost::geometry)