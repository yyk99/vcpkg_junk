# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Vulkan graphics programming primer project containing three example applications that demonstrate Vulkan API usage with GLFW for window management. The project showcases basic Vulkan initialization, classic triangle rendering with a complete graphics pipeline, and advanced 3D cube rendering.

## Build System

### Dependencies
- vcpkg (package manager)
- CMake 3.10+
- C++17 compiler
- Vulkan SDK
- GLFW3 (via vcpkg)
- GLM (via vcpkg)

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

**Manual Configuration:**
```bash
cmake -S . -B build-vs2022-x64 -DCMAKE_TOOLCHAIN_FILE=c:/opt/vcpkg/scripts/buildsystems/vcpkg.cmake
```

### Building
```bash
# After configuring with presets:
cmake --build --preset=kestrel         # Linux (GCC)
cmake --build --preset=kestrel-clang   # Linux (Clang)
cmake --build --preset=default         # Windows

# Or manually specify build directory:
cmake --build build-linux
cmake --build build-vs2022
```

### Running the Applications
```bash
# Basic Vulkan window application
./build-linux/vulkan_app

# Classic triangle with vertex colors
./build-linux/vulkan_triangle

# 3D rotating cube with lighting
./build-linux/vulkan_cube
```

## Code Architecture

### Applications

**vulkan_app.cpp** - Basic Vulkan Application
- Demonstrates fundamental Vulkan initialization pipeline
- Creates a window with GLFW
- Sets up Vulkan instance, physical device, logical device
- Creates window surface and swap chain
- Minimal example showing core Vulkan setup steps

**vulkan_triangle.cpp** - Triangle Rendering (Hello Triangle)
- Classic introductory Vulkan example with complete graphics pipeline
- Renders a single triangle with per-vertex colors (red, green, blue)
- Implements vertex buffer with position and color attributes
- Uses embedded SPIR-V shader bytecode (no external shader files needed)
- Demonstrates full render pass, pipeline, framebuffer, and command buffer setup
- Simpler than cube: no transformations, no depth buffer, just 3 vertices
- Good starting point for understanding the Vulkan rendering pipeline

**vulkan_cube.cpp** - 3D Cube Rendering
- Full Vulkan graphics pipeline implementation
- Renders a rotating 3D cube with colored faces (red, green, blue, yellow, magenta, cyan)
- Implements vertex and index buffers
- Uses uniform buffers for MVP (Model-View-Projection) matrices
- Includes depth buffering for proper 3D rendering
- Demonstrates descriptor sets for shader uniforms
- Implements double buffering with synchronization (MAX_FRAMES_IN_FLIGHT = 2)

### Key Data Structures

**Vertex** (vulkan_cube.cpp:15-42)
- Position (vec3) and color (vec3) attributes
- Provides binding and attribute descriptions for Vulkan pipeline

**UniformBufferObject** (vulkan_cube.cpp:44-48)
- Contains model, view, and projection matrices
- Updated per-frame for cube rotation animation

**Cube Geometry** (vulkan_cube.cpp:51-97)
- 24 vertices (4 per face, 6 faces)
- 36 indices (2 triangles × 3 vertices per face × 6 faces)
- Each face has a distinct color

### Vulkan Pipeline Setup

All applications follow similar initialization patterns:
1. Window creation (GLFW)
2. Vulkan instance creation
3. Surface creation for window rendering
4. Physical device selection
5. Logical device and queue creation
6. Swap chain setup

The triangle application additionally implements:
7. Render pass configuration
8. Graphics pipeline with embedded SPIR-V shaders
9. Framebuffer creation
10. Vertex buffer (3 vertices with position + color)
11. Command buffer recording
12. Synchronization objects (semaphores and fences)

The cube application extends this with:
13. Index buffer for efficient geometry
14. Uniform buffers for MVP matrices
15. Descriptor sets for shader resources
16. Depth buffering for proper 3D occlusion

## Code Style

The project includes a `.clang-format` file for consistent code formatting. Use `clang-format` to format code before committing changes.
