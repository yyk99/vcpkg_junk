# =========================================================================
#   Unity - A Test Framework for C
#   ThrowTheSwitch.org
#   Copyright (c) 2007-24 Mike Karlesky, Mark VanderVoord, & Greg Williams
#   SPDX-License-Identifier: MIT
# =========================================================================

cmake_minimum_required(VERSION 3.20)
project(unity_example_1 C)

# Set C standard to C89 (equivalent to -std=c89)
set(CMAKE_C_STANDARD 90)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# Unity root directory (relative to this CMakeLists.txt)
set(UNITY_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../Unity)

# Compiler warning flags (matching makefile)
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

# Optional: Add -Werror for CI builds
option(UNITY_WERROR "Treat warnings as errors" OFF)
if(UNITY_WERROR)
    list(APPEND WARNING_FLAGS -Werror)
endif()

# Test executable 1: TestProductionCode
add_executable(test1
    ${UNITY_ROOT}/src/unity.c
    src/ProductionCode.c
    test/TestProductionCode.c
    test/test_runners/TestProductionCode_Runner.c
)

target_include_directories(test1 PRIVATE
    src
    ${UNITY_ROOT}/src
)

target_compile_options(test1 PRIVATE ${WARNING_FLAGS})

# Test executable 2: TestProductionCode2
add_executable(test2
    ${UNITY_ROOT}/src/unity.c
    src/ProductionCode2.c
    test/TestProductionCode2.c
    test/test_runners/TestProductionCode2_Runner.c
)

target_include_directories(test2 PRIVATE
    src
    ${UNITY_ROOT}/src
)

target_compile_options(test2 PRIVATE ${WARNING_FLAGS})

# Register tests with CTest
add_test(NAME test1 COMMAND test1)
add_test(NAME test2 COMMAND test2)

# Optional: Custom commands to regenerate test runners using Ruby
# Uncomment if you want to regenerate runners automatically
# find_program(RUBY_EXECUTABLE ruby)
# if(RUBY_EXECUTABLE)
#     add_custom_command(
#         OUTPUT test/test_runners/TestProductionCode_Runner.c
#         COMMAND ${RUBY_EXECUTABLE} ${UNITY_ROOT}/auto/generate_test_runner.rb
#                 ${CMAKE_CURRENT_SOURCE_DIR}/test/TestProductionCode.c
#                 ${CMAKE_CURRENT_SOURCE_DIR}/test/test_runners/TestProductionCode_Runner.c
#         DEPENDS test/TestProductionCode.c
#         COMMENT "Generating TestProductionCode_Runner.c"
#     )
#
#     add_custom_command(
#         OUTPUT test/test_runners/TestProductionCode2_Runner.c
#         COMMAND ${RUBY_EXECUTABLE} ${UNITY_ROOT}/auto/generate_test_runner.rb
#                 ${CMAKE_CURRENT_SOURCE_DIR}/test/TestProductionCode2.c
#                 ${CMAKE_CURRENT_SOURCE_DIR}/test/test_runners/TestProductionCode2_Runner.c
#         DEPENDS test/TestProductionCode2.c
#         COMMENT "Generating TestProductionCode2_Runner.c"
#     )
# endif()
