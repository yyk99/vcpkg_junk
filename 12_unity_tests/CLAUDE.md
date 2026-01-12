# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Unity test framework example project that demonstrates how to use the Unity C testing framework (from ThrowTheSwitch/Unity) with CMake. The Unity framework is included as a Git submodule in the `Unity/` directory.

## Build System

The project uses CMake with custom presets defined in `CMakePresets.json`:

**Configure:**
```bash
# Debug build (default)
cmake --preset debug

# Release build
cmake --preset release
```

**Build:**
```bash
# Build debug
cmake --build --preset debug

# Build release
cmake --build --preset release
```

**Build outputs** are organized as follows:
- Executables: `build/{debug|release}/bin/`
- Libraries: `build/{debug|release}/lib/`

## Testing

**Run all tests via CTest:**
```bash
ctest --preset debug
```

**Run test executable directly:**
```bash
./build/debug/bin/test_calculator
```

## Project Architecture

### Unity Framework Integration

Unity is included as a Git submodule and integrated via `add_subdirectory(Unity)` in the root CMakeLists.txt. The Unity library is linked using the `unity::framework` target.

### Test Structure

Tests follow this pattern (see `test_calculator.c` as reference):
- Include `unity.h` and the header of the code under test
- Implement `setUp()` and `tearDown()` functions (run before/after each test)
- Write test functions with `void test_name(void)` signature
- Use Unity assertions like `TEST_ASSERT_EQUAL(expected, actual)`
- In `main()`, wrap tests with `UNITY_BEGIN()` and `UNITY_END()`
- Register tests with `RUN_TEST(test_function_name)`

### Adding New Tests

When adding a new module to test:
1. Create the production code as a library in CMakeLists.txt
2. Create a test file (e.g., `test_module.c`)
3. Add test executable and link it with the library and `unity::framework`
4. Register the test with `add_test()` in CMakeLists.txt
