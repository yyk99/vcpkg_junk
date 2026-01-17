cmake_minimum_required(VERSION 3.15)
project(cmd_tree_from_inside_predicate)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(BUILD_TESTING)
     enable_testing()

     include(GoogleTest)

     list(APPEND VCPKG_MANIFEST_FEATURES "tests")

     find_package(GTest REQUIRED)

    # Add your source files here
    add_executable(cmd_tree_from_inside_predicate 
        cmd_tree_from_inside_predicate.cpp)

    # Find Boost (adjust components as needed)
    find_package(Boost REQUIRED COMPONENTS geometry)

    target_link_libraries(cmd_tree_from_inside_predicate PRIVATE GTest::gtest)
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

endif()



