#
#
#

cmake_minimum_required(VERSION 3.20)
project (vcpkg_junk)

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

enable_testing()

# end of file
