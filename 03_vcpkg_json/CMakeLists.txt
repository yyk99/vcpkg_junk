#
# See more details
# https://learn.microsoft.com/en-us/vcpkg/users/buildsystems/cmake-integration#settings-reference
#

cmake_minimum_required(VERSION 3.25)

if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()

option(BUILD_TESTING "Build tests" ON)
if(BUILD_TESTING)
    list(APPEND VCPKG_MANIFEST_FEATURES "tests")
endif()

project(03_vcpkg_json)

# Prevent a "command line is too long" failure in Windows.
set(CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "Force Ninja to use response files.")

set_property(GLOBAL PROPERTY USE_FOLDERS TRUE)

find_package(Boost REQUIRED COMPONENTS
    program_options
    filesystem
    geometry
    polygon
)

find_package(VTK COMPONENTS
    CommonColor
    CommonCore
    CommonDataModel
    FiltersCore
    FiltersSources
    InteractionStyle
    IOImport
    RenderingContextOpenGL2
    RenderingCore
    RenderingFreeType
    RenderingGL2PSOpenGL2
    RenderingOpenGL2
)

if(NOT VTK_FOUND)
    message(FATAL_ERROR "Unable to find the VTK build folder.")
endif()


add_definitions(${Boost_LIB_DIAGNOSTIC_DEFINITIONS})

add_executable(
    cmd_getargs_03
    cmd_getargs.cpp
    ../.clang-format
    README.md
)
target_link_libraries(cmd_getargs_03 PRIVATE Boost::program_options)

add_executable(
    voronoi_basic_tutorial
    voronoi_basic_tutorial.cpp
    ../.clang-format
    README.md
)
target_link_libraries(voronoi_basic_tutorial PRIVATE Boost::geometry Boost::polygon)

add_executable(
    ConstrainedDelaunay2D
    ConstrainedDelaunay2D.cxx
    ../.clang-format
    README.md
)
target_link_libraries(ConstrainedDelaunay2D PRIVATE ${VTK_LIBRARIES} Boost::geometry Boost::polygon)
vtk_module_autoinit(TARGETS ConstrainedDelaunay2D MODULES ${VTK_LIBRARIES})

add_executable(
    Triangulate
    Triangulate.cxx
    ../.clang-format
    README.md
)
target_link_libraries(Triangulate PRIVATE ${VTK_LIBRARIES} Boost::geometry Boost::polygon)
vtk_module_autoinit(TARGETS Triangulate MODULES ${VTK_LIBRARIES})

set_target_properties(
    Triangulate ConstrainedDelaunay2D voronoi_basic_tutorial cmd_getargs_03
    PROPERTIES FOLDER "Apps")

add_executable(
    GLTFImporter
    GLTFImporter.cxx
    ../.clang-format
    README.md
)
target_link_libraries(GLTFImporter PRIVATE ${VTK_LIBRARIES})
vtk_module_autoinit(TARGETS GLTFImporter MODULES ${VTK_LIBRARIES})
set_target_properties(GLTFImporter PROPERTIES FOLDER "Apps")

if(BUILD_TESTING)
    enable_testing()

    include(GoogleTest)

    list(APPEND VCPKG_MANIFEST_FEATURES "tests")

    find_package(GTest REQUIRED)

    add_executable(test_boost_03
        test_boost.cpp
    )
    target_link_libraries(
        test_boost_03
        PRIVATE
        Boost::program_options
        Boost::filesystem
        Boost::geometry
        GTest::gtest
        GTest::gtest_main
    )

    gtest_discover_tests(test_boost_03)

    set_target_properties(test_boost_03 PROPERTIES FOLDER "Tests")
endif()
