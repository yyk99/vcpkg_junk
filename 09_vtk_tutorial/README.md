# General info

Display the use of vcpkg.json file

- kestrel preset uses "central" location to install required packages

## Preconditions

Install vcpkg

    git clone https://github.com/microsoft/vcpkg.git
    cd vcpkg
    ./bootstrap-vcpkg.sh # in Unix/Linux

## How to configure

In Linux

    cmake --preset kestrel
    cmake --build build-linux
    cmake --build build-linux --target test

In Windows (WIP)

    cmake --preset default
    cmake --build build-vs2022
    cmake --build build-vs2022 --target RUN_TESTS

## Setting environment

In Windows:

    In the Command Line (CMD)

    set PROJ_DATA=vcpkg_installed/x64-windows/share/proj/
    Debug\test_09.exe

If Visual Studio is used, then set Debugging / Environment property:

    PROJ_DATA=$(ProjectDir)vcpkg_installed/x64-windows/share/proj/

## REFERENCE

See more details
https://learn.microsoft.com/en-us/vcpkg/users/buildsystems/cmake-integration#settings-reference

VTX Tutorial Step1 (and so on)
https://examples.vtk.org/site/Cxx/Tutorial/Tutorial_Step1/


