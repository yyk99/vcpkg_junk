#
#
#

cmake_minimum_required(VERSION 3.25)
project(cmake_with_boost)

set_property(GLOBAL PROPERTY USE_FOLDERS TRUE)

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost REQUIRED COMPONENTS 
	program_options
	filesystem
)

#message(STATUS "Boost_LIBRARY_DIRS=${Boost_LIBRARY_DIRS}")
#message(STATUS "Boost_LIB_DIAGNOSTIC_DEFINITIONS=${Boost_LIB_DIAGNOSTIC_DEFINITIONS}")
#message(STATUS "VCPKG_TARGET_TRIPLET=${VCPKG_TARGET_TRIPLET}")
#    
#add_definitions(${Boost_LIB_DIAGNOSTIC_DEFINITIONS})

add_executable(
	cmd_getargs_02
	cmd_getargs.cpp
	../.clang-format
	README.txt
)
target_link_libraries(cmd_getargs_02 PRIVATE Boost::program_options)
set_target_properties(cmd_getargs_02 PROPERTIES FOLDER "Apps")


