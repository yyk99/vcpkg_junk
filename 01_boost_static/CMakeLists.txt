#
#
#

cmake_minimum_required(VERSION 3.28)
project(cmake_with_boost)

option (BUILD_STATIC "Build with static boost libraries" ON)

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_STATIC_RUNTIME ON)
find_package(Boost REQUIRED COMPONENTS program_options geometry)

#message(STATUS "Boost_LIBRARY_DIRS=${Boost_LIBRARY_DIRS}")
#message(STATUS "Boost_LIB_DIAGNOSTIC_DEFINITIONS=${Boost_LIB_DIAGNOSTIC_DEFINITIONS}")
#message(STATUS "VCPKG_TARGET_TRIPLET=${VCPKG_TARGET_TRIPLET}")
#    
#add_definitions(${Boost_LIB_DIAGNOSTIC_DEFINITIONS})

add_executable(
	cmd_getargs 
	cmd_getargs.cpp
	.clang-format
	README.txt
)
target_link_libraries(cmd_getargs PRIVATE Boost::program_options Boost::geometry)
if (BUILD_STATIC)
	set_property(TARGET cmd_getargs PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()