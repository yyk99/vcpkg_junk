#
# See more details
# https://learn.microsoft.com/en-us/vcpkg/users/buildsystems/cmake-integration#settings-reference
#

cmake_minimum_required(VERSION 3.25)
project(04_cgal_tbb)

find_package(Boost REQUIRED COMPONENTS 
	program_options
	filesystem
)

find_package(TBB REQUIRED)
find_package(CGAL REQUIRED)
include (CGAL_TBB_support)

add_executable(
	cmd_getargs_06
	cmd_getargs.cpp
	../.clang-format
	README.txt
)
target_link_libraries(cmd_getargs_06 PRIVATE Boost::program_options)

option(BUILD_TESTING "Build tests" ON)
if(BUILD_TESTING)
	enable_testing()
	include(GoogleTest)

	set (BUILD_GMOCK OFF CACHE BOOL "")
	list(APPEND VCPKG_MANIFEST_FEATURES "tests")

	find_package(GTest)
	if (NOT GTEST_FOUND)
		message(STATUS "Build gtest from source")
		include(FetchContent)
		FetchContent_Declare(
  			googletest
  			# Specify the commit you depend on and update it regularly.
  			URL https://github.com/google/googletest/archive/5376968f6948923e2411081fd9372e71a59d8e77.zip
		)
		# For Windows: Prevent overriding the parent project's compiler/linker settings
		set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
		FetchContent_MakeAvailable(googletest)

		# Now simply link against gtest or gtest_main as needed. Eg
		## add_executable(example example.cpp)
		## target_link_libraries(example gtest_main)
		## add_test(NAME example_test COMMAND example)
	endif()

	add_executable(test_boost_06 
		test_boost.cpp
	)
	target_link_libraries(
		test_boost_06
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
		CGAL_target_use_TBB(test_boost_06)
	endif()
	gtest_discover_tests(test_boost_06)

endif()
