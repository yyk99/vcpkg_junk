# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a CGAL (Computational Geometry Algorithms Library) testing project that demonstrates various CGAL features including polygon operations, boolean set operations, and geometric computations. It's part of a larger vcpkg_junk_bmg repository containing multiple vcpkg-related test projects.

## Build System

This project uses CMake with vcpkg for dependency management. Multiple CMake presets are configured for different platforms and environments.

### Building on Linux

```bash
# Configure with the kestrel preset (Linux/Unix Makefiles)
cmake --preset kestrel

# Build all targets
cmake --build build-linux

# Build a specific example
cmake --build build-linux --target cmd_cgal_example1

# Run tests
cmake --build build-linux --target test_12
ctest --test-dir build-linux
```

### Building on Windows

```bash
# Configure with default preset (Visual Studio 2022)
cmake --preset default

# Build all targets
cmake --build build-vs2022
```

### CMake Presets

- **kestrel**: Linux build using Unix Makefiles with vcpkg root at `/home/yyk/src/vcpkg`
- **default**: Windows build inheriting from vcpkg preset with VCPKG_ROOT at `g:\opt\vcpkg`
- **windows-developmentyk**: Windows alternative with VCPKG_ROOT at `c:\opt\vcpkg`
- **cherry**: Windows alternative with VCPKG_ROOT at `e:\opt\vcpkg`

Note: VCPKG_ROOT environment variable must be set correctly for your system before building.

## Project Structure

### Executables

All example executables are built using the `add_cmd_cgal_example()` macro, which links against CGAL::CGAL and outputs to `build-*/bin/`:

- **cmd_cgal_example1**: Basic polygon operations (area, convexity, point-in-polygon testing)
- **cmd_cgal_2**: Polygons with holes demonstration
- **cmd_cgal_3**: Boolean set operations (union, intersection, difference)
- **cmd_cgal_4, cmd_cgal_5, cmd_cgal_6**: Additional CGAL examples (WIP)

### Tests

- **test_12**: GTest-based unit test that uses CGAL and Boost.Geometry
- Tests are only built when `BUILD_TESTING=ON` (default)
- Uses GoogleTest discovery via `gtest_discover_tests()`

## Dependencies (vcpkg.json)

Core dependencies:
- **cgal**: The main CGAL library
- **fmt**: Formatting library
- **boost-geometry**: Boost geometry library

Test dependencies (feature "tests"):
- **gtest**: Google Test framework

## Key CGAL Patterns

All examples use the `Exact_predicates_inexact_constructions_kernel` which provides:
- Exact geometric predicates (reliable geometric tests)
- Inexact constructions (fast but approximate coordinates)
- Common typedefs: `K`, `Point_2`, `Polygon_2`, `Polygon_with_holes_2`

## Running Examples

After building, run examples from the build directory:

```bash
# Linux
./build-linux/bin/cmd_cgal_example1
./build-linux/bin/cmd_cgal_2

# Windows
.\build-vs2022\bin\Debug\cmd_cgal_example1.exe
```

## CMake Features

- Output directories are unified: binaries go to `bin/`, libraries to `lib/`
- Targets are organized in IDE folders: "Apps" for examples, "Tests" for test executables
- Ninja response files are forced on to prevent command line length issues on Windows
- Policy CMP0167 is set to NEW if available
