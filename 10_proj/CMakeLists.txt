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

project(10_proj)

# Prevent a "command line is too long" failure in Windows.
set(CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "Force Ninja to use response files.")

set_property(GLOBAL PROPERTY USE_FOLDERS TRUE)

find_package(Boost REQUIRED COMPONENTS
    # program_options
    # filesystem
    geometry
    # polygon
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
    message(FATAL_ERROR "ConstrainedDelaunay2D: Unable to find the VTK build folder.")
endif()

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

    add_executable(test_10
        test_boost.cpp
    )
    target_link_libraries(
        test_10
        PRIVATE
        Boost::geometry
        GTest::gtest
        GTest::gtest_main
    )

    gtest_discover_tests(test_10)

    set_target_properties(test_10 PROPERTIES FOLDER "Tests")
endif()
