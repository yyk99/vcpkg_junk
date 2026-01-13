# Unity test framework examples

This project demonstrates how to use the Unity test framework (from ThrowTheSwitch/Unity) with CMake.

## Prerequisites

Before building this project, ensure you have:

- **CMake** 3.20 or higher
- **C compiler** (GCC, Clang, or MSVC)
- **Ruby interpreter** - Required for auto-generating test runners
  - Ubuntu/Debian: `sudo apt install ruby`
  - macOS: `brew install ruby` (or use system Ruby)
  - Windows: Download from [ruby-lang.org](https://www.ruby-lang.org/en/downloads/)

Verify Ruby installation:
```bash
ruby --version
```

## Project Structure

- `calculator.h` / `calculator.c` - Simple calculator functions to test
- `test_calculator.c` - Unity tests for the calculator (test runners auto-generated)
- `CMakeLists.txt` - CMake build configuration
- `CMakePresets.json` - CMake presets for easy configuration
- `Unity/` - Unity test framework (Git submodule)
- `examples/` - Unity example tests from the framework

## Building and Running

### Using Unix Makefiles (Linux/MinGW/Cygwin)

**Note:** Unix Makefiles preset can be used on Windows if `make` is available (e.g., from Cygwin, MinGW, or MSYS2).

```bash
# Configure with CMake preset
cmake --preset debug

# Build
cmake --build --preset debug

# Run tests
ctest --preset debug

# Or run the test executable directly
./build/debug/bin/test_calculator
```

### Using Clang Compiler

```bash
# Configure with Clang (Debug)
cmake --preset clang-debug

# Build
cmake --build --preset clang-debug

# Run tests
ctest --preset clang-debug

# Or run the test executable directly
./build/clang-debug/bin/test_calculator
```

### Using Visual Studio 2022 (Windows)

```bash
# Configure with Visual Studio 2022
cmake --preset vs2022

# Build (Debug configuration)
cmake --build --preset vs2022 --config Debug

# Build (Release configuration)
cmake --build --preset vs2022 --config Release

# Run tests (uses Debug configuration by default)
ctest --preset vs2022

# Run tests with specific configuration
ctest --preset vs2022 -C Release
```

## Test Runner Generation

This project uses Unity's Ruby-based test runner generator. During the build:

1. CMake invokes `Unity/auto/generate_test_runner.rb`
2. The script scans test files for functions matching `test_*` pattern
3. A test runner file is generated in the build directory (e.g., `test_calculator_Runner.c`)
4. The runner automatically calls all discovered tests with proper `setUp()`/`tearDown()`

**Benefits:**
- No manual test registration needed
- Just write test functions starting with `test_` and they're automatically included
- Test runners regenerate when test files change

## Example Output

When all tests pass, you should see:
```
test_calculator.c:12:test_add_positive_numbers:PASS
test_calculator.c:17:test_add_negative_numbers:PASS
test_calculator.c:23:test_subtract:PASS
test_calculator.c:30:test_multiply:PASS

-----------------------
4 Tests 0 Failures 0 Ignored
OK
```
