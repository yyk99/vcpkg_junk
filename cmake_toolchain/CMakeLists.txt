#
#
#

cmake_minimum_required(VERSION 3.28)
project(cmake_with_boost)

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost REQUIRED COMPONENTS program_options)

message(STATUS "Boost_LIBRARY_DIRS=${Boost_LIBRARY_DIRS}")
message(STATUS "Boost_LIB_DIAGNOSTIC_DEFINITIONS=${Boost_LIB_DIAGNOSTIC_DEFINITIONS}")
    
add_definitions(${Boost_LIB_DIAGNOSTIC_DEFINITIONS})

add_executable(
	cmd_getargs 
	cmd_getargs.cpp
	../.clang-format
	README.txt
)
target_link_libraries(cmd_getargs PRIVATE Boost::program_options)



