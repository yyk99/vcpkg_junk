# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a CMake-based C++ project exploring spatial data structures and algorithms, with a primary focus on R-trees, octrees, and point cloud processing. The project uses Boost Geometry for spatial indexing, PCL (Point Cloud Library) for 3D point cloud operations, and fmt for formatted output. It consists of 12 command-line applications demonstrating different spatial computing concepts, from low-level custom implementations to high-level library usage:

- **cmd_fmt_example**: Comprehensive fmt library demonstration including text formatting, number formatting (binary/hex/octal), string alignment, custom type formatting, and locale-aware output
- **cmd_boost_rtree_primer**: Comprehensive Boost R-tree spatial indexing primer demonstrating 2D/3D operations, k-nearest neighbor searches, range queries, custom predicates, performance comparisons, and visitor pattern traversal
- **cmd_boost_rtree_traversal**: Multiple R-tree traversal methods including iterator-based, query-based, visitor pattern, nearest neighbor traversal, and depth-first simulation with tree statistics analysis
- **cmd_boost_rtree_view_example**: R-tree view utilities featuring tree structure visualization, memory layout analysis, spatial distribution analysis, and ASCII grid visualization with tree analyzer for height/branching estimation
- **cmd_cesium_3d_tiles_pointcloud**: Point cloud to Cesium 3D Tiles conversion with PLY file loading, spatial partitioning, hierarchical tile generation, and geometric error computation
- **cmd_mem_allocate**: Memory allocation stress testing using exponential doubling strategy to determine system memory limits (C implementation)
- **cmd_octree_primer**: Custom octree data structure implementation with adaptive subdivision based on point density, range queries, nearest point finding, and spatial search optimization
- **cmd_pcl_octree_primer**: PCL octree examples including voxel-based searches, k-nearest neighbor queries, radius-based range searches, and change detection experiments
- **cmd_rtree_primer**: Generic R-tree implementation from scratch with custom node structures, split algorithms, bounding box operations, and performance analysis
- **cmd_simple_boost_rstar_example**: Simplified Boost R* tree usage demonstrating bulk insertion, range queries, k-NN searches, intersection queries, custom predicates, and point removal
- **cmd_sombrero_elevation_generator**: Mathematical sombrero function elevation data generator with multiple PNG output formats (grayscale, color-mapped, 16-bit heightmap) and CSV export capabilities
- **cmd_tree_from_inside_predicate**: Hierarchical tree construction using containment predicates with automatic parent-child relationship establishment and tree reparenting (includes gtest integration)

## Build System

The project uses **vcpkg** in manifest mode for dependency management with Visual Studio integration.

### Build Commands

**Visual Studio (Windows):**
```bash
# Configure with vcpkg preset
cmake --preset vcpkg
# Build all targets
cmake --build build-vs2022
# Or open in Visual Studio
start 11_vs.sln
```

**Linux/GCC:**
```bash
# Manual configure with vcpkg
mkdir build-linux && cd build-linux
cmake -DCMAKE_TOOLCHAIN_FILE=$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake ..
# Build
cmake --build .
```

### Running Individual Applications
```bash
# From build directory
./bin/cmd_fmt_example
./bin/cmd_boost_rtree_primer
# etc.
```

## Dependencies

Main dependencies managed through `vcpkg.json`:
- **fmt**: String formatting library
- **pcl**: Point Cloud Library for 3D processing
- **stb**: Image processing utilities
- **boost-geometry**: Spatial algorithms and data structures
- **gtest**: Testing framework (optional, enabled with BUILD_TESTING)

## Architecture Notes

- Each `cmd_*` directory contains a self-contained application with its own CMakeLists.txt
- All applications use C++17 standard
- Build outputs go to `build-*/bin/` for executables and `build-*/lib/` for libraries
- Project uses folder organization in IDE with "Apps 11_vs" folder grouping
- MSVC builds include `/utf-8` compiler flag for Unicode support

### Key Technical Themes

The applications demonstrate several core concepts:
- **Spatial Indexing**: R-tree variants (Linear, Quadratic, R*) for efficient spatial queries
- **Point Cloud Processing**: 3D point data manipulation, partitioning, and format conversion
- **Query Optimization**: K-nearest neighbor, range queries, intersection tests, and custom predicates
- **Hierarchical Structures**: Octrees, tile hierarchies, and containment-based tree building
- **Data Export**: Multiple output formats including PNG (grayscale, color-mapped, 16-bit), CSV, and 3D Tiles
- **Performance Analysis**: Comparative benchmarks between custom implementations and library-based solutions

## Development Workflow

1. Make changes to source files in respective `cmd_*` directories
2. Build specific target: `cmake --build build-vs2022 --target <target_name>`
3. Run executable from build directory
4. For new applications, add subdirectory to root CMakeLists.txt and include in set_target_properties

## vcpkg Integration

Initial setup requires:
```bash
cd cmd_fmt_example  # or any subdirectory with vcpkg usage
vcpkg integrate install
```

The project uses vcpkg baseline `6ecbbbdf31cba47aafa7cf6189b1e73e10ac61f8` for reproducible builds.