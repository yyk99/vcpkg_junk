#
# See more details
# https://learn.microsoft.com/en-us/vcpkg/users/buildsystems/cmake-integration#settings-reference
#

if(POLICY CMP0167)
	cmake_policy(SET CMP0167 NEW)
endif()

cmake_minimum_required(VERSION 3.25)
project(assimp_junk)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Boost REQUIRED COMPONENTS 
	program_options
	filesystem
	geometry
)

find_package(assimp REQUIRED)

add_definitions(${Boost_LIB_DIAGNOSTIC_DEFINITIONS})

add_executable(
	cmd_getargs_08
	cmd_getargs.cpp
	../.clang-format
	README.md
)
target_link_libraries(cmd_getargs_08 PRIVATE Boost::program_options)

option(BUILD_TESTING "Build tests" ON)
if(BUILD_TESTING)
	enable_testing()
	
	include(GoogleTest)

	list(APPEND VCPKG_MANIFEST_FEATURES "tests")

	find_package(GTest REQUIRED)

	add_executable(test_boost_08
		test_boost.cpp
		test_assimp.cpp
	)
	target_link_libraries(
		test_boost_08
		PRIVATE
		Boost::program_options
		Boost::filesystem
		Boost::geometry
		assimp::assimp
		GTest::gtest
		GTest::gtest_main
	)

	gtest_discover_tests(test_boost_08)

endif()
