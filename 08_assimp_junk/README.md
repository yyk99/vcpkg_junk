# General info

Display the use of vcpkg.json file

## Preconditions

Install vcpkg

    git clone .... blah blah blah

## How to configure

In Linux

    cmake --preset kestrel
    cmake --build build-linux
    cmake --build build-linux --target test

In Windows (WIP)

    cmake --preset default -DCMAKE_POLICY_DEFAULT_CMP0167=NEW
    cmake --build build-vs2022
    cmake --build build-vs2022 --target RUN_TESTS

## REFERENCE

See more details
https://learn.microsoft.com/en-us/vcpkg/users/buildsystems/cmake-integration#settings-reference


