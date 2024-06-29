#
# See more details
# https://learn.microsoft.com/en-us/vcpkg/users/buildsystems/cmake-integration#settings-reference
#

cmake_minimum_required(VERSION 3.28)
project(05_manifest_features)

option(BUILD_WITH_VCPKG_CGAL "Build with CGAL from VCPKG" OFF)
option(BUILD_TESTING "Build tests" OFF)
set(VCPKG_MANIFEST_FEATURES "" CACHE STRING "Manifest features")

if (BUILD_WITH_VCPKG_CGAL)
	list(APPEND VCPKG_MANIFEST_FEATURES "cgal-vcpkg")
endif()

if (BUILD_TESTING)
	list(APPEND VCPKG_MANIFEST_FEATURES "tests")
endif()

find_package(Boost REQUIRED COMPONENTS 
	program_options
	filesystem
)

find_package(TBB REQUIRED)
find_package(CGAL REQUIRED)
include (CGAL_TBB_support)

if (CGAL_FOUND)
	message(STATUS "CGAL found at CGAL_DIR=${CGAL_DIR}")
endif()

add_executable(
	cmd_getargs 
	cmd_getargs.cpp
	../.clang-format
	README.txt
)
target_link_libraries(cmd_getargs PRIVATE Boost::program_options)

if(BUILD_TESTING)
	enable_testing()

	find_package(GTest REQUIRED)

	add_executable(test_boost 
		test_boost.cpp
	)
	target_link_libraries(
		test_boost
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
	add_test(test_boost test_boost)

endif()
