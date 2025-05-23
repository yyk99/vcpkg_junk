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

find_package(PROJ CONFIG REQUIRED)

find_package(Eigen3 CONFIG REQUIRED)

# set_target_properties(GLTFImporter PROPERTIES FOLDER "Apps")

if(BUILD_TESTING)
    enable_testing()

    include(GoogleTest)

    list(APPEND VCPKG_MANIFEST_FEATURES "tests")

    find_package(GTest REQUIRED)

    add_executable(test_10
        test_boost.cpp
        test_proj.cpp
        lv95_globalmapper.h
        lv95_good.h
        lv95_original.h
        test_eigen.cpp
        README.md
    )

    target_compile_definitions(test_10
        PRIVATE
        "PROJ_DIR=\"${PROJ_DIR}\""
    )

    target_link_libraries(
        test_10
        PRIVATE
        Boost::geometry
        PROJ::proj
        Eigen3::Eigen
        GTest::gtest
        GTest::gtest_main
    )

    gtest_discover_tests(test_10)

    set_target_properties(test_10 PROPERTIES FOLDER "Tests")
endif()
