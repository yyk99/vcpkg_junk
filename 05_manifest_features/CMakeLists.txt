#
# See more details
# https://learn.microsoft.com/en-us/vcpkg/users/buildsystems/cmake-integration#settings-reference
#

cmake_minimum_required(VERSION 3.28)

option(BUILD_WITH_VCPKG_CGAL "Build with CGAL from VCPKG" OFF)
option(BUILD_TESTING "Build tests" OFF)
# set(VCPKG_MANIFEST_FEATURES "" CACHE STRING "Manifest features")

if (BUILD_WITH_VCPKG_CGAL)
	list(APPEND VCPKG_MANIFEST_FEATURES "cgal-vcpkg")
else()
	if (EXISTS ${CMAKE_SOURCE_DIR}/cgal)
		message(STATUS "---------- set GDAL_DIR ----------")
		set(GDAL_DIR "${CMAKE_SOURCE_DIR}/cgal" CACHE PATH "")
		list(APPEND VCPKG_MANIFEST_FEATURES "cgal-local")
	endif()
endif()

if (BUILD_TESTING)
	list(APPEND VCPKG_MANIFEST_FEATURES "tests")
endif()


project(05_manifest_features)

find_package(Boost REQUIRED COMPONENTS 
	program_options
	filesystem
)

message (STATUS "*********** VCPKG_MANIFEST_FEATURES=${VCPKG_MANIFEST_FEATURES} **************")

find_package(TBB REQUIRED)
find_package(CGAL REQUIRED)
include (CGAL_TBB_support)

if (CGAL_FOUND)
	message(STATUS "CGAL found at CGAL_DIR=${CGAL_DIR}")
endif()

add_executable(
	cmd_getargs_05
	cmd_getargs.cpp
	../.clang-format
	README.txt
	vcpkg.json
)
target_link_libraries(cmd_getargs_05 PRIVATE Boost::program_options)

if(BUILD_TESTING)
	enable_testing()

	find_package(GTest REQUIRED)

	add_executable(test_boost_05
		test_boost.cpp
	)
	target_link_libraries(
		test_boost_05
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
		CGAL_target_use_TBB(test_boost_05)
	endif()
	add_test(test_boost_05 test_boost_05)

endif()
