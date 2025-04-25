#
#
#

cmake_minimum_required(VERSION 3.25)
option (BUILD_STATIC "Build with static boost libraries" ON)
option (BUILD_TESTING "Build test app" ON)

project(cmake_with_boost)

set_property(GLOBAL PROPERTY USE_FOLDERS TRUE)

set(Boost_USE_STATIC_LIBS ${BUILD_STATIC})
set(Boost_USE_STATIC_RUNTIME ${BUILD_STATIC})
find_package(Boost REQUIRED COMPONENTS 
	program_options
	filesystem
	geometry
)

#message(STATUS "Boost_LIBRARY_DIRS=${Boost_LIBRARY_DIRS}")
#message(STATUS "Boost_LIB_DIAGNOSTIC_DEFINITIONS=${Boost_LIB_DIAGNOSTIC_DEFINITIONS}")
#message(STATUS "VCPKG_TARGET_TRIPLET=${VCPKG_TARGET_TRIPLET}")
#    
#add_definitions(${Boost_LIB_DIAGNOSTIC_DEFINITIONS})

add_executable(
	cmd_getargs 
	cmd_getargs.cpp
	.clang-format
	README.md
)
target_link_libraries(cmd_getargs PRIVATE Boost::program_options Boost::geometry)
if (BUILD_STATIC)
	set_property(TARGET cmd_getargs PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

if(BUILD_TESTING)
	enable_testing()
	
	include(GoogleTest)

	list(APPEND VCPKG_MANIFEST_FEATURES "tests")

	find_package(GTest REQUIRED)

	add_executable(test_boost_01
		test_boost.cpp
	)
	target_link_libraries(
		test_boost_01
		PRIVATE 
		Boost::program_options 
		Boost::filesystem
		Boost::geometry
		GTest::gtest 
		GTest::gtest_main
	)

	gtest_discover_tests(test_boost_01)

	set_target_properties(test_boost_01 PROPERTIES FOLDER "Tests")
endif()