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

# Add the test executable
add_executable(test_calculator test_calculator.c)

# Link the calculator library and Unity
target_link_libraries(test_calculator PRIVATE calculator unity::framework)

# Enable testing
enable_testing()
add_test(NAME calculator_tests COMMAND test_calculator)
