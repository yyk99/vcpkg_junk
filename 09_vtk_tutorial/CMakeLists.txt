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

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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
    #    cli11
    CommonColor
    CommonCore
    CommonDataModel
    CommonTransforms
    FiltersCore
    FiltersSources
    fmt
    InteractionStyle
    InteractionWidgets
    IOCesium3DTiles
    IOGeometry
    IOImport
    RenderingAnnotation
    RenderingContextOpenGL2
    RenderingCore
    RenderingFreeType
    RenderingGL2PSOpenGL2
    RenderingLabel
    RenderingOpenGL2
)

if(NOT VTK_FOUND)
    message(FATAL_ERROR "Unable to find the VTK build folder.")
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

macro(add_vtk_example tutorial_name)
    add_executable(${tutorial_name} MACOSX_BUNDLE ${tutorial_name}.cxx)
    target_link_libraries(${tutorial_name}
        PRIVATE ${VTK_LIBRARIES} Boost::program_options
    )
    vtk_module_autoinit(
        TARGETS ${tutorial_name}
        MODULES ${VTK_LIBRARIES}
    )
    set_target_properties(${tutorial_name} PROPERTIES FOLDER "VTK Examples")
endmacro()

add_vtk_tutorial(Tutorial_Step1)
add_vtk_tutorial(Tutorial_Step2)
add_vtk_tutorial(Tutorial_Step3)
add_vtk_tutorial(Tutorial_Step4)
add_vtk_tutorial(Tutorial_Step5)
add_vtk_tutorial(Tutorial_Step6)
add_vtk_tutorial(Tutorial_Step7_build_mesh)
add_vtk_example(Polygon)
add_vtk_example(Pyramid)
add_vtk_example(PolyLine1)
add_vtk_example(PolyLine1_from_boost)
add_vtk_example(LinearCellsDemo)
if (NOT TARGET GLTFImporter)
    add_vtk_example(GLTFImporter)
endif()
add_vtk_example(Axis)

if(BUILD_TESTING)
    # this is heuristically generated, and may not be correct
    find_package(jsoncpp CONFIG REQUIRED)
    # target_link_libraries(main PRIVATE JsonCpp::JsonCpp)

    enable_testing()

    include(GoogleTest)

    list(APPEND VCPKG_MANIFEST_FEATURES "tests")

    find_package(GTest REQUIRED)

    add_executable(test_09
        test_boost.cpp
        test_vtk.cpp
        README.md
    )
    target_link_libraries(
        test_09
        PRIVATE
        Boost::program_options
        Boost::filesystem
        Boost::geometry
        GTest::gtest
        GTest::gtest_main
        ${VTK_LIBRARIES}
        JsonCpp::JsonCpp
    )

    gtest_discover_tests(test_09)

    set_target_properties(test_09 PROPERTIES FOLDER "Tests")
endif()
