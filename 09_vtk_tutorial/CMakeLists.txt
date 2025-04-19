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

project(09_vtk_tutorial)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

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
    CommonTransforms
    FiltersCore
    FiltersSources
    InteractionStyle
    InteractionWidgets
    RenderingContextOpenGL2
    RenderingCore
    RenderingFreeType
    RenderingGL2PSOpenGL2
    RenderingOpenGL2
)

if(NOT VTK_FOUND)
    message(FATAL_ERROR "ConstrainedDelaunay2D: Unable to find the VTK build folder.")
endif()

add_definitions(${Boost_LIB_DIAGNOSTIC_DEFINITIONS})

macro(add_vtk_tutorial tutorial_name)
    add_executable(${tutorial_name} MACOSX_BUNDLE ${tutorial_name}.cxx)
    target_link_libraries(${tutorial_name} PRIVATE ${VTK_LIBRARIES})
    vtk_module_autoinit(
        TARGETS ${tutorial_name}
        MODULES ${VTK_LIBRARIES}
    )
    set_target_properties(${tutorial_name} PROPERTIES FOLDER "VTK Tutorials")
endmacro()

add_vtk_tutorial(Tutorial_Step1)
add_vtk_tutorial(Tutorial_Step2)
add_vtk_tutorial(Tutorial_Step3)
add_vtk_tutorial(Tutorial_Step4)
add_vtk_tutorial(Tutorial_Step5)
add_vtk_tutorial(Tutorial_Step6)

if(BUILD_TESTING)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(JSONCPP jsoncpp)

    enable_testing()

    include(GoogleTest)

    list(APPEND VCPKG_MANIFEST_FEATURES "tests")

    find_package(GTest REQUIRED)

    add_executable(test_boost_09
        test_boost.cpp
        test_vtk.cpp
    )
    target_link_libraries(
        test_boost_09
        PRIVATE
        Boost::program_options
        Boost::filesystem
        Boost::geometry
        GTest::gtest
        GTest::gtest_main
        ${VTK_LIBRARIES}
        ${JSONCPP_LIBRARIES}
    )

    gtest_discover_tests(test_boost_09)

    set_target_properties(test_boost_09 PROPERTIES FOLDER "Tests")
endif()
