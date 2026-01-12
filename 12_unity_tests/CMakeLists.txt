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

# Add Unity example_1 targets
# Production code libraries
add_library(production_code Unity/examples/example_1/src/ProductionCode.c Unity/examples/example_1/src/ProductionCode.h)
target_include_directories(production_code PUBLIC Unity/examples/example_1/src)

add_library(production_code2 Unity/examples/example_1/src/ProductionCode2.c Unity/examples/example_1/src/ProductionCode2.h)
target_include_directories(production_code2 PUBLIC Unity/examples/example_1/src)

# Test executables for example_1
add_executable(test_production_code
    Unity/examples/example_1/test/TestProductionCode.c
    Unity/examples/example_1/test/test_runners/TestProductionCode_Runner.c
)
target_link_libraries(test_production_code PRIVATE production_code unity::framework)

add_executable(test_production_code2
    Unity/examples/example_1/test/TestProductionCode2.c
    Unity/examples/example_1/test/test_runners/TestProductionCode2_Runner.c
)
target_link_libraries(test_production_code2 PRIVATE production_code2 unity::framework)

# Register tests with CTest
add_test(NAME production_code_tests COMMAND test_production_code)
add_test(NAME production_code2_tests COMMAND test_production_code2)
