cmake_minimum_required(VERSION 3.15)
project(cmd_pcl_octree_primer)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add your source files here
add_executable(cmd_pcl_octree_primer 
    cmd_pcl_octree_primer.cpp)

# Find Boost (adjust components as needed)
find_package(PCL REQUIRED)

target_include_directories(cmd_pcl_octree_primer PRIVATE ${PCL_LIBRARIES})
target_link_libraries(cmd_pcl_octree_primer PRIVATE ${PCL_LIBRARIES})
