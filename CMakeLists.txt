#
#
#

cmake_minimum_required(VERSION 3.22)

if(POLICY CMP0167)
	cmake_policy(SET CMP0167 NEW)
endif()

option(BUILD_TESTING "Build tests" ON)
if(BUILD_TESTING)
  list(APPEND VCPKG_MANIFEST_FEATURES "tests")
endif()

option(BUILD_01 "Build 01_boost_static" ON)
if(BUILD_01)
  list(APPEND VCPKG_MANIFEST_FEATURES "01-boost-static")
endif()

option(BUILD_02 "Build 02_cmake_toolchain" ON)
if(BUILD_02)
  list(APPEND VCPKG_MANIFEST_FEATURES "02-cmake-toolchain")
endif()

option(BUILD_03 "Build 03_vcpkg_json" ON)
if(BUILD_03)
  list(APPEND VCPKG_MANIFEST_FEATURES "03-vcpkg-json")
endif()

option(BUILD_04 "Build 04_cgal_tbb" ON)
if(BUILD_04)
  list(APPEND VCPKG_MANIFEST_FEATURES "04-cgal-tbb")
endif()

option(BUILD_05 "Build 05_manifest_features" ON)
if(BUILD_05)
  list(APPEND VCPKG_MANIFEST_FEATURES "05-manifest-features")
endif()

option(BUILD_06 "Build 06_cgal_tbb_novcpkg" ON)
if(BUILD_06)
  list(APPEND VCPKG_MANIFEST_FEATURES "06-cgal-tbb-novcpkg")
endif()

option(BUILD_07 "Build 07_cmake_qt" ON)
if(BUILD_07)
  list(APPEND VCPKG_MANIFEST_FEATURES "07-cmake-qt")
endif()

option(BUILD_08 "Build 08_assimp_junk" ON)
if(BUILD_08)
  list(APPEND VCPKG_MANIFEST_FEATURES "08-assimp-junk")
endif()

option(BUILD_09 "Build 09_vtk_tutorial" ON)
if(BUILD_09)
  list(APPEND VCPKG_MANIFEST_FEATURES "09-vtk-tutorial")
endif()

option(BUILD_10 "Build 10_proj" ON)
if(BUILD_10)
  list(APPEND VCPKG_MANIFEST_FEATURES "10-proj")
endif()

option(BUILD_11 "Build 11_vs" ON)
if(BUILD_11)
  list(APPEND VCPKG_MANIFEST_FEATURES "11-vs")
endif()

project (vcpkg_junk)

if(BUILD_TESTING)
	enable_testing()
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if(BUILD_01)
	add_subdirectory(01_boost_static)
endif()
if(BUILD_02)
	add_subdirectory(02_cmake_toolchain)
endif()
if(BUILD_03)
	add_subdirectory(03_vcpkg_json)
endif()
if(BUILD_04)
	add_subdirectory(04_cgal_tbb)
endif()
if(BUILD_05)
	add_subdirectory(05_manifest_features)
endif()
if(BUILD_06)
	add_subdirectory(06_cgal_tbb_novcpkg)
endif()
if(BUILD_07)
	add_subdirectory(07_cmake_qt)
endif()
if(BUILD_08)
	add_subdirectory(08_assimp_junk)
endif()
if(BUILD_09)
	add_subdirectory(09_vtk_tutorial)
endif()
if(BUILD_10)
	message("Include 10_proj...")
	add_subdirectory(10_proj)
endif()
if(BUILD_11)
	message("Include 11_vs...")
	add_subdirectory(11_vs)
endif()

# end of file
