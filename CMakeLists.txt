#
#
#

cmake_minimum_required(VERSION 3.25)

if(POLICY CMP0167)
	cmake_policy(SET CMP0167 NEW)
endif()

option(BUILD_TESTING "Build tests" ON)
if(BUILD_TESTING)
  list(APPEND VCPKG_MANIFEST_FEATURES "tests")
endif()

project (vcpkg_junk)

if(BUILD_TESTING)
	enable_testing()
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(01_boost_static)
add_subdirectory(02_cmake_toolchain)
add_subdirectory(03_vcpkg_json)
add_subdirectory(04_cgal_tbb)
add_subdirectory(05_manifest_features)
add_subdirectory(06_cgal_tbb_novcpkg)
add_subdirectory(07_cmake_qt)
add_subdirectory(08_assimp_junk)
add_subdirectory(09_vtk_tutorial)

# end of file
