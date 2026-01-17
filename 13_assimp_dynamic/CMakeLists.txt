#
# See more details
# https://learn.microsoft.com/en-us/vcpkg/users/buildsystems/cmake-integration#settings-reference
#

cmake_minimum_required(VERSION 3.22)

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

find_package(assimp REQUIRED)
find_package(Stb REQUIRED)
find_package(draco REQUIRED)

if(BUILD_TESTING)
    enable_testing()
    include(GoogleTest)
    find_package(GTest REQUIRED)

    add_executable(test_13
        test_assimp.cpp
        test_draco.cpp
        test_draco.h

        .editorconfig
    )
    target_include_directories(test_13 PRIVATE ${Stb_INCLUDE_DIR})
    target_link_libraries(test_13 PRIVATE
        assimp::assimp
        draco::draco

        GTest::gtest
        GTest::gtest_main
    )
    set_target_properties(test_13 PROPERTIES FOLDER "Tests")
    gtest_discover_tests(test_13)

endif()
