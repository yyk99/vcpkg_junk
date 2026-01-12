# Unity test framework examples

This project demonstrates how to use the Unity test framework (from ThrowTheSwitch/Unity) with CMake and vcpkg.

## Project Structure

- `calculator.h` / `calculator.c` - Simple calculator functions to test
- `test_calculator.c` - Unity tests for the calculator
- `CMakeLists.txt` - CMake build configuration
- `CMakePresets.json` - CMake presets for easy configuration
- `vcpkg.json` - vcpkg dependencies

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
