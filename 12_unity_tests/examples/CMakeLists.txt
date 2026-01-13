# Unity example_1 targets
# Production code libraries
add_library(production_code example_1/src/ProductionCode.c example_1/src/ProductionCode.h)
target_include_directories(production_code PUBLIC example_1/src)

add_library(production_code2 example_1/src/ProductionCode2.c example_1/src/ProductionCode2.h)
target_include_directories(production_code2 PUBLIC example_1/src)

# Test executables for example_1
add_executable(test_production_code
    example_1/test/TestProductionCode.c
    example_1/test/test_runners/TestProductionCode_Runner.c
)
target_link_libraries(test_production_code PRIVATE production_code unity::framework)

add_executable(test_production_code2
    example_1/test/TestProductionCode2.c
    example_1/test/test_runners/TestProductionCode2_Runner.c
)
target_link_libraries(test_production_code2 PRIVATE production_code2 unity::framework)

# Register tests with CTest
add_test(NAME production_code_tests COMMAND test_production_code)
add_test(NAME production_code2_tests COMMAND test_production_code2)
