#
# See more details
# https://learn.microsoft.com/en-us/vcpkg/users/buildsystems/cmake-integration#settings-reference
#

cmake_minimum_required(VERSION 3.25)

option(BUILD_TESTING "Build tests" ON)
if(BUILD_TESTING)
  list(APPEND VCPKG_MANIFEST_FEATURES "tests")
endif()

project(04_cgal_tbb)

find_package(Boost REQUIRED COMPONENTS 
	program_options
	filesystem
)

find_package(TBB REQUIRED)
find_package(CGAL REQUIRED)
include (CGAL_TBB_support)

add_executable(
	cmd_getargs_04
	cmd_getargs.cpp
	../.clang-format
	README.txt
)
target_link_libraries(cmd_getargs_04 PRIVATE Boost::program_options)

if(BUILD_TESTING)
	enable_testing()
	include(GoogleTest)

	list(APPEND VCPKG_MANIFEST_FEATURES "tests")

	find_package(GTest REQUIRED)

	add_executable(test_boost_04
		test_boost.cpp
	)
	target_link_libraries(
		test_boost_04
		PRIVATE 
		Boost::program_options 
		Boost::filesystem

		CGAL::CGAL 
		CGAL::TBB_support

		GTest::gtest 
		GTest::gtest_main
	)
	if (FALSE)
		# DEPRECATED
		CGAL_target_use_TBB(test_boost)
	endif()
	gtest_discover_tests(test_boost_04)

endif()
