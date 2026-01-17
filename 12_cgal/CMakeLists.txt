cmake_minimum_required(VERSION 3.22)

if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()

option(BUILD_TESTING "Build tests" ON)
if(BUILD_TESTING)
    list(APPEND VCPKG_MANIFEST_FEATURES "tests")
endif()

project(12_cgal)

# Prevent a "command line is too long" failure in Windows.
set(CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "Force Ninja to use response files.")

set_property(GLOBAL PROPERTY USE_FOLDERS TRUE)

find_package(Boost REQUIRED COMPONENTS
    geometry
)

set (CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)
find_package(CGAL REQUIRED)
message(STATUS "CGAL version: ${CGAL_VERSION}")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


macro (add_cmd_cgal_example _cmd_name)
    add_executable (
        ${_cmd_name}
        ${_cmd_name}.cpp 
    )

    target_link_libraries(
        ${_cmd_name}
        PRIVATE
        CGAL::CGAL
    )
    set_target_properties(${_cmd_name} PROPERTIES FOLDER "Apps")
endmacro(add_cmd_cgal_example)

add_cmd_cgal_example(cmd_cgal_example1)
add_cmd_cgal_example(cmd_cgal_2)
add_cmd_cgal_example(cmd_cgal_3)
add_cmd_cgal_example(cmd_cgal_4)
add_cmd_cgal_example(cmd_cgal_5)
add_cmd_cgal_example(cmd_cgal_6)

if(BUILD_TESTING)
    enable_testing()

    include(GoogleTest)

    list(APPEND VCPKG_MANIFEST_FEATURES "tests")

    find_package(GTest REQUIRED)

    add_executable(test_12
        test_12.cpp
    )
    target_link_libraries(
        test_12
        PRIVATE
        CGAL::CGAL
        Boost::geometry
        GTest::gtest
        GTest::gtest_main
    )

    gtest_discover_tests(test_12)

    set_target_properties(test_12 PROPERTIES FOLDER "Tests")
endif()
