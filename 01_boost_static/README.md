# 01_boost_static

## Preconditions

- install vcpkg
- install boots required packages...

In Windows:

    vcpkg install boost-filesystem --triplet x64-windows-static
    vcpkg install boost-geometry --triplet x64-windows-static
    vcpkg install boost-program_options --triplet x64-windows-static
    vcpkg install gtest --triplet x64-windows-static

---------------- cut here ----------------------
Total install time: 4.8 min
The package boost-program-options is compatible with built-in CMake targets of FindBoost.cmake:

    find_package(Boost REQUIRED COMPONENTS program_options)
    target_link_libraries(main PRIVATE Boost::program_options)

or the generated cmake configs via:

    find_package(boost_program_options REQUIRED CONFIG)
    target_link_libraries(main PRIVATE Boost::program_options)
---------------- cut here ----------------------


## How to configure

Windows

    cmake -S . -B build-vs2022-x64 -DCMAKE_TOOLCHAIN_FILE=E:/opt/vcpkg/scripts/buildsystems/vcpkg.cmake -DVCPKG_TARGET_TRIPLET=x64-windows-static

Linux (x64 or arm64)

It seems that vcpkg in Linux is static (by default). Nevertheless:

    cmake -S . -B build-linux -DCMAKE_TOOLCHAIN_FILE=/home/yyk/src/vcpkg/scripts/buildsystems/vcpkg.cmake


## How to build

    cmake --build build-vs2022-x64
