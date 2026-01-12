# Unity test framework examples

This project demonstrates how to use the Unity test framework (from ThrowTheSwitch/Unity) with CMake and vcpkg.

## Project Structure

- `calculator.h` / `calculator.c` - Simple calculator functions to test
- `test_calculator.c` - Unity tests for the calculator
- `CMakeLists.txt` - CMake build configuration
- `CMakePresets.json` - CMake presets for easy configuration
- `vcpkg.json` - vcpkg dependencies

## Building and Running

```bash
# Configure with CMake preset
cmake --preset default

# Build
cmake --build --preset default

# Run tests
ctest --preset default

# Or run the test executable directly
./build/test_calculator
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
