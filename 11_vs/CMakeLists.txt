cmake_minimum_required(VERSION 3.22)

if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()

option(BUILD_TESTING "Build tests" ON)
if(BUILD_TESTING)
    list(APPEND VCPKG_MANIFEST_FEATURES "tests")
endif()

project(11_vs_cmake)

# Prevent a "command line is too long" failure in Windows.
set(CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "Force Ninja to use response files.")

set_property(GLOBAL PROPERTY USE_FOLDERS TRUE)

find_package(Boost REQUIRED COMPONENTS
    geometry
)

find_package(PCL REQUIRED
)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(cmd_app1)
add_subdirectory(cmd_boost_rtree_primer)
add_subdirectory(cmd_boost_rtree_traversal)
add_subdirectory(cmd_boost_rtree_view_example)
add_subdirectory(cmd_cesium_3d_tiles_pointcloud)
add_subdirectory(cmd_mem_allocate)
add_subdirectory(cmd_octree_primer)
add_subdirectory(cmd_pcl_octree_primer)
add_subdirectory(cmd_rtree_primer)
add_subdirectory(cmd_simple_boost_rstar_example)

set_target_properties(
    cmd_app1
    cmd_boost_rtree_primer
    cmd_boost_rtree_traversal
    cmd_boost_rtree_view_example
    cmd_cesium_3d_tiles_pointcloud
    cmd_mem_allocate
    cmd_octree_primer
    cmd_pcl_octree_primer
    cmd_rtree_primer
    cmd_simple_boost_rstar_example
    PROPERTIES FOLDER "Apps 11_vs")

# add_definitions(${Boost_LIB_DIAGNOSTIC_DEFINITIONS})

# add_executable(
#     cmd_getargs_03
#     cmd_getargs.cpp
#     ../.clang-format
#     README.md
# )
# target_link_libraries(cmd_getargs_03 PRIVATE Boost::program_options)

# add_executable(
#     voronoi_basic_tutorial
#     voronoi_basic_tutorial.cpp
#     ../.clang-format
#     README.md
# )
# target_link_libraries(voronoi_basic_tutorial PRIVATE Boost::geometry Boost::polygon)

# add_executable(
#     ConstrainedDelaunay2D
#     ConstrainedDelaunay2D.cxx
#     ../.clang-format
#     README.md
# )
# target_link_libraries(ConstrainedDelaunay2D PRIVATE ${VTK_LIBRARIES} Boost::geometry Boost::polygon)
# vtk_module_autoinit(TARGETS ConstrainedDelaunay2D MODULES ${VTK_LIBRARIES})

# add_executable(
#     Triangulate
#     Triangulate.cxx
#     ../.clang-format
#     README.md
# )
# target_link_libraries(Triangulate PRIVATE ${VTK_LIBRARIES} Boost::geometry Boost::polygon)
# vtk_module_autoinit(TARGETS Triangulate MODULES ${VTK_LIBRARIES})

# set_target_properties(
#     Triangulate ConstrainedDelaunay2D voronoi_basic_tutorial cmd_getargs_03
#     PROPERTIES FOLDER "Apps")

# add_executable(
#     GLTFImporter
#     GLTFImporter.cxx
#     ../.clang-format
#     README.md
# )
# target_link_libraries(GLTFImporter PRIVATE ${VTK_LIBRARIES})
# vtk_module_autoinit(TARGETS GLTFImporter MODULES ${VTK_LIBRARIES})
# set_target_properties(GLTFImporter PROPERTIES FOLDER "Apps")

# if(BUILD_TESTING)
#     enable_testing()

#     include(GoogleTest)

#     list(APPEND VCPKG_MANIFEST_FEATURES "tests")

#     find_package(GTest REQUIRED)

#     add_executable(test_boost_03
#         test_boost.cpp
#     )
#     target_link_libraries(
#         test_boost_03
#         PRIVATE
#         Boost::program_options
#         Boost::filesystem
#         Boost::geometry
#         GTest::gtest
#         GTest::gtest_main
#     )

#     gtest_discover_tests(test_boost_03)

#     set_target_properties(test_boost_03 PROPERTIES FOLDER "Tests")
# endif()
