# =========================================================================
#   Unity - A Test Framework for C
#   ThrowTheSwitch.org
#   Copyright (c) 2007-24 Mike Karlesky, Mark VanderVoord, & Greg Williams
#   SPDX-License-Identifier: MIT
# =========================================================================

cmake_minimum_required(VERSION 3.20)
project(unity_example_5 C)

# Set C standard to C89 (equivalent to -std=c89)
set(CMAKE_C_STANDARD 90)
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

# Test executable: TestProductionCode
add_executable(ex5_test1
    ${UNITY_ROOT}/src/unity.c
    src/ProductionCode.c
    test/TestProductionCode.c
    test/test_runners/TestProductionCode_Runner.c
)

target_include_directories(ex5_test1 PRIVATE
    src
    ${UNITY_ROOT}/src
)

target_compile_definitions(ex5_test1 PRIVATE
    UNIT_TESTING
)

target_compile_options(ex5_test1 PRIVATE ${WARNING_FLAGS})

# Force include unity_detail_config.h for all source files
# This enables Unity's Details Stack feature for enhanced error reporting
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(ex5_test1 PRIVATE
        -include ${CMAKE_CURRENT_SOURCE_DIR}/test/unity_detail_config.h
    )
elseif(CMAKE_C_COMPILER_ID MATCHES "MSVC")
    target_compile_options(ex5_test1 PRIVATE
        /FI${CMAKE_CURRENT_SOURCE_DIR}/test/unity_detail_config.h
    )
endif()

# Register test with CTest
add_test(NAME ex5_test1 COMMAND ex5_test1)
