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

project(cmake_with_boost)

find_package(Boost REQUIRED COMPONENTS 
	program_options
	filesystem
	geometry
)

add_definitions(${Boost_LIB_DIAGNOSTIC_DEFINITIONS})

add_executable(
	cmd_getargs_03
	cmd_getargs.cpp
	../.clang-format
	README.md
)
target_link_libraries(cmd_getargs_03 PRIVATE Boost::program_options)

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

endif()
