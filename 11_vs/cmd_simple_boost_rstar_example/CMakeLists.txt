cmake_minimum_required(VERSION 3.15)
project(cmd_simple_boost_rstar_example)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add your source files here
add_executable(cmd_simple_boost_rstar_example 
    cmd_simple_boost_rstar_example.cpp)

# Find Boost (adjust components as needed)
find_package(Boost REQUIRED COMPONENTS geometry)

target_link_libraries(cmd_simple_boost_rstar_example PRIVATE Boost::geometry)
