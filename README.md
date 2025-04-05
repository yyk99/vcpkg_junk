# Collection of vcpkg related project

## How to Build in Linux

	cmake --preset kestrel
	cmake --build build-linux

## How to build in Windows

	cmake --preset default -DCMAKE_POLICY_DEFAULT_CMP0167=NEW
	cmake --build build-vs2022

## EXAMPLE command line(s)

	cmake -S 02_cmake_toolchain/ -B build-ubuntu-vcpkg -DCMAKE_TOOLCHAIN_FILE=/home/yyk/src/vcpkg/scripts/buildsystems/vcpkg.cmake
