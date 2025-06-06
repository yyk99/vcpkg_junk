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

project(assimp_junk)

set_property(GLOBAL PROPERTY USE_FOLDERS TRUE)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Boost REQUIRED COMPONENTS 
    program_options
    filesystem
    geometry
)

find_package(assimp REQUIRED)

find_package(Stb REQUIRED)

find_package(json-c CONFIG REQUIRED)

find_package(Eigen3 CONFIG REQUIRED)

add_definitions(${Boost_LIB_DIAGNOSTIC_DEFINITIONS})

add_executable(
    cmd_getargs_08
    cmd_getargs.cpp
    ../.clang-format
    README.md
)
target_link_libraries(cmd_getargs_08 PRIVATE Boost::program_options)
set_target_properties(cmd_getargs_08 PROPERTIES FOLDER "Apps")

if(BUILD_TESTING)
    enable_testing()
    
    include(GoogleTest)

    list(APPEND VCPKG_MANIFEST_FEATURES "tests")

    find_package(GTest REQUIRED)

    add_executable(test_08
        TilesetJson.h
        assimp_aux.h
        meshtoolbox.h
        test_assimp.cpp
        test_boost.cpp

        .editorconfig
    )
    target_include_directories(test_08 PRIVATE ${Stb_INCLUDE_DIR})
    target_link_libraries(
        test_08
        PRIVATE
        Boost::program_options
        Boost::filesystem
        Boost::geometry
        assimp::assimp
        json-c::json-c
        Eigen3::Eigen

        GTest::gtest
        GTest::gtest_main
    )
    set_target_properties(test_08 PROPERTIES FOLDER "Tests")
    gtest_discover_tests(test_08)

endif()
