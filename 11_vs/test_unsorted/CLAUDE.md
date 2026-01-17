# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a C++ research/testing project focused on exploring virtual base class behavior and inheritance patterns in C++. The project uses vcpkg for dependency management and GoogleTest for unit testing.

## Build System

### Dependencies
- vcpkg (package manager)
- CMake 3.15+
- C++17 compiler
- GoogleTest (via vcpkg)

### Configuration

**Using CMake Presets (Recommended):**
```bash
# Windows - use one of: default, windows-developmentyk, cherry
cmake --preset=default

# Linux - use kestrel (GCC) or kestrel-clang (Clang)
cmake --preset=kestrel
cmake --preset=kestrel-clang
```

Available presets in CMakeUserPresets.json:
- `default`: Windows (G:\opt\vcpkg)
- `windows-developmentyk`: Windows (C:\opt\vcpkg)
- `cherry`: Windows (E:\opt\vcpkg)
- `kestrel`: Linux with GCC (/home/yyk/src/vcpkg)
- `kestrel-clang`: Linux with Clang (/home/yyk/src/vcpkg)

**Manual Configuration (Alternative):**

Windows:
```bash
cmake -S . -B build-vs2022-x64 -DCMAKE_TOOLCHAIN_FILE=E:/opt/vcpkg/scripts/buildsystems/vcpkg.cmake -DBUILD_TESTING=ON
```

Linux (x64 or arm64):
```bash
cmake -S . -B build-linux -DCMAKE_TOOLCHAIN_FILE=/home/yyk/src/vcpkg/scripts/buildsystems/vcpkg.cmake
```

Linux with Clang:
```bash
cmake -S . -B build-linux -DCMAKE_TOOLCHAIN_FILE=/home/yyk/src/vcpkg/scripts/buildsystems/vcpkg.cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++
```

### Building
```bash
# After configuring with presets:
cmake --build --preset=kestrel         # Linux (GCC)
cmake --build --preset=kestrel-clang   # Linux (Clang)
cmake --build --preset=default         # Windows

# Or manually specify build directory:
cmake --build build-linux
cmake --build build-linux-clang
cmake --build build-vs2022
```

### Running Tests
```bash
cd build-linux && ctest
# or run the test executable directly:
./build-linux/test_unsorted
```

## Code Architecture

### Core Components

**Widget and Data Base Classes (virtual_base_test.cpp)**
- `Widget`: Simple handle wrapper class
- `Data`: Base class with virtual functions and `property1` field

**Inheritance Pattern Testing**
The project compares three different inheritance approaches through parallel implementations:

1. **Plain Multiple Inheritance** (namespace `plain`, virtual_base_test_plain.h)
   - `Dialog`: inherits from both `Widget` and `Data` (non-virtual)
   - `DialogMore`: inherits from both `Dialog` and `DataMore` (non-virtual)
   - Results in diamond problem: two separate `Data` subobjects in `DialogMore`
   - MSVC and Clang behave differently when accessing ambiguous base members

2. **Virtual Inheritance** (namespace `virt`, virtual_base_test_virt.h)
   - `Dialog`: inherits from `Widget` and `virtual Data`
   - `DataMore`: inherits from `virtual Data`
   - `DialogMore`: inherits from both `Dialog` and `DataMore`
   - Resolves diamond problem: single shared `Data` subobject in `DialogMore`

3. **Pointer-Based Design** (namespace `ptr`, in virtual_base_test.cpp)
   - `Dialog` contains `std::unique_ptr<Data>`
   - `DialogMore` replaces pointer with `DataMore` instance
   - Uses `dynamic_cast` and virtual `data()` accessor
   - Avoids multiple inheritance complexity entirely

### Test Structure

All tests are in the `VirtualBaseF` test fixture:
- `t1`: Tests plain multiple inheritance patterns and demonstrates ambiguity
- `t2`: Tests virtual inheritance resolving the diamond problem
- `t3`: Demonstrates virtual base class dominance with weak/dominant pattern
- `t4`: Tests the pointer-based design alternative

### Key Macros

Defined in test_unsorted.h:
- `CONSOLE(x)`: Debug-only console output with function name
- `CONSOLE_EVAL(x)`: Evaluates and prints expression with its value

## Platform-Specific Notes

The code includes platform-specific checks:
- MSVC-specific expectations in tests (e.g., ambiguous base access behavior)
- Tests account for different behavior between MSVC and Clang (`#if _MSC_VER && !defined(__clang__)`)

## Branch Structure

- Main branch: `master`
- Development branch: `devel/yurik42/11_vs/test_unsorted`
