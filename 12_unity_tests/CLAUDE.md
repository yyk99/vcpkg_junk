# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Unity test framework example project that demonstrates how to use the Unity C testing framework (from ThrowTheSwitch/Unity) with CMake. The Unity framework is included as a Git submodule in the `Unity/` directory.

## Build System

The project uses CMake with custom presets defined in `CMakePresets.json`.

**Note:** Unix Makefiles presets (debug/release) can be used on Windows if `make` is available (e.g., from Cygwin, MinGW, or MSYS2). For native Windows development, use the vs2022 preset.

**Configure:**
```bash
# Debug build (Unix Makefiles)
cmake --preset debug

# Release build (Unix Makefiles)
cmake --preset release

# Visual Studio 2022 (multi-config generator)
cmake --preset vs2022
```

**Build:**
```bash
# Build debug (Unix Makefiles)
cmake --build --preset debug

# Build release (Unix Makefiles)
cmake --build --preset release

# Build with Visual Studio 2022 (specify config with --config)
cmake --build --preset vs2022 --config Debug
cmake --build --preset vs2022 --config Release
```

**Build outputs** are organized as follows:
- Unix Makefiles: Executables in `build/{debug|release}/bin/`, libraries in `build/{debug|release}/lib/`
- Visual Studio: Multi-config output in `build/vs2022/bin/{Debug|Release}/` and `build/vs2022/lib/{Debug|Release}/`

## Testing

**Run all tests via CTest:**
```bash
# Unix Makefiles
ctest --preset debug

# Visual Studio 2022 (uses Debug config by default)
ctest --preset vs2022

# Visual Studio 2022 with specific configuration
ctest --preset vs2022 -C Release
```

**Run test executable directly:**
```bash
# Unix Makefiles
./build/debug/bin/test_calculator

# Visual Studio 2022
./build/vs2022/bin/Debug/test_calculator.exe
```

**Note:** Visual Studio is a multi-configuration generator, so you must specify the configuration (Debug/Release) when building and can override it when testing using `-C <config>`.

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
