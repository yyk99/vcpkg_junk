# =========================================================================
#   Unity - A Test Framework for C
#   ThrowTheSwitch.org
#   Copyright (c) 2007-24 Mike Karlesky, Mark VanderVoord, & Greg Williams
#   SPDX-License-Identifier: MIT
# =========================================================================

cmake_minimum_required(VERSION 3.20)
project(unity_example_2 C)

# Set C standard to C99 (equivalent to -std=c99)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# Unity root directory (relative to this CMakeLists.txt)
set(UNITY_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../Unity)

# Compiler warning flags (compiler-dependent)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    # GCC/Clang warning flags (matching makefile)
    set(WARNING_FLAGS
        -Wall
        -Wextra
        -Wpointer-arith
        -Wcast-align
        -Wwrite-strings
        -Wswitch-default
        -Wunreachable-code
        -Winit-self
        -Wmissing-field-initializers
        -Wno-unknown-pragmas
        -Wstrict-prototypes
        -Wundef
        -Wold-style-definition
    )
elseif(CMAKE_C_COMPILER_ID MATCHES "MSVC")
    # MSVC warning flags
    set(WARNING_FLAGS
        /W4                 # High warning level (equivalent to -Wall -Wextra)
        /wd4068             # Disable unknown pragma warning (equivalent to -Wno-unknown-pragmas)
    )
endif()

# Optional: Add -Werror/WX for CI builds
option(UNITY_WERROR "Treat warnings as errors" OFF)
if(UNITY_WERROR)
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
        list(APPEND WARNING_FLAGS -Werror)
    elseif(CMAKE_C_COMPILER_ID MATCHES "MSVC")
        list(APPEND WARNING_FLAGS /WX)
    endif()
endif()

# Test executable: all_tests (uses Unity fixture extension)
add_executable(ex2_all_tests
    ${UNITY_ROOT}/src/unity.c
    ${UNITY_ROOT}/extras/fixture/src/unity_fixture.c
    src/ProductionCode.c
    src/ProductionCode2.c
    test/TestProductionCode.c
    test/TestProductionCode2.c
    test/test_runners/TestProductionCode_Runner.c
    test/test_runners/TestProductionCode2_Runner.c
    test/test_runners/all_tests.c
)

target_include_directories(ex2_all_tests PRIVATE
    src
    ${UNITY_ROOT}/src
    ${UNITY_ROOT}/extras/fixture/src
)

target_compile_definitions(ex2_all_tests PRIVATE
    UNITY_FIXTURE_NO_EXTRAS
)

target_compile_options(ex2_all_tests PRIVATE ${WARNING_FLAGS})

# Register test with CTest
add_test(NAME ex2_all_tests COMMAND ex2_all_tests -v)
