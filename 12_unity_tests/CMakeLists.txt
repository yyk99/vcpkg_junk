cmake_minimum_required(VERSION 3.20)
project(unity_example C)

set(CMAKE_C_STANDARD 11)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Add Unity test framework from submodule
add_subdirectory(Unity)

# Add the calculator library
add_library(calculator calculator.c calculator.h)

# Find Ruby for generating test runners
find_program(RUBY_EXECUTABLE ruby REQUIRED)

# Generate test runner for test_calculator.c
set(TEST_CALCULATOR_RUNNER ${CMAKE_BINARY_DIR}/test_calculator_Runner.c)
add_custom_command(
    OUTPUT ${TEST_CALCULATOR_RUNNER}
    COMMAND ${RUBY_EXECUTABLE} ${CMAKE_SOURCE_DIR}/Unity/auto/generate_test_runner.rb
            ${CMAKE_SOURCE_DIR}/test_calculator.c
            ${TEST_CALCULATOR_RUNNER}
    DEPENDS test_calculator.c
    COMMENT "Generating test runner for test_calculator.c"
    VERBATIM
)

# Add the test executable with generated runner
add_executable(test_calculator
    test_calculator.c
    ${TEST_CALCULATOR_RUNNER}
)

# Add include directories so the generated runner can find headers
target_include_directories(test_calculator PRIVATE ${CMAKE_SOURCE_DIR})

# Link the calculator library and Unity
target_link_libraries(test_calculator PRIVATE calculator unity::framework)

# Enable testing
enable_testing()
add_test(NAME calculator_tests COMMAND test_calculator)

# Add examples subdirectory
add_subdirectory(examples)
