# Shader Compilation Instructions for vulkan_cube

## Problem
The `vulkan_cube.exe` requires compiled SPIR-V shader files (`vert.spv` and `frag.spv`) but they are not included in the repository.

## Solution

### Option 1: Install Vulkan SDK and Compile Shaders (Recommended)

1. **Install Vulkan SDK:**
   - Download from: https://vulkan.lunarg.com/sdk/home
   - Install to default location (e.g., `C:\VulkanSDK\1.x.x.x\`)

2. **Compile shaders:**
   - Run the provided batch file:
     ```batch
     compile_shaders.bat
     ```
   - Or manually compile:
     ```batch
     glslc shader.vert -o vert.spv
     glslc shader.frag -o frag.spv
     ```

3. **Copy compiled shaders to build directory:**
   ```batch
     copy vert.spv build-vs2022-x64\Release\
   copy frag.spv build-vs2022-x64\Release\
   ```

4. **Run the application:**
   ```batch
   cd build-vs2022-x64\Release
   vulkan_cube.exe
   ```

### Option 2: Use Online SPIR-V Compiler

If you cannot install the Vulkan SDK:

1. Visit: https://shader-playground.timjones.io/
2. Paste the contents of `shader.vert` (vertex shader)
3. Select "GLSL" as input and "SPIR-V" as output
4. Compile and download as `vert.spv`
5. Repeat for `shader.frag` → `frag.spv`
6. Place both `.spv` files in `build-vs2022-x64\Release\`

### Option 3: Use Embedded SPIR-V (Future Enhancement)

Similar to `vulkan_triangle.cpp`, the shaders could be embedded directly in the source code to avoid external file dependencies. This would require valid SPIR-V bytecode arrays.

## Shader Source Files

The GLSL shader source files are included:
- `shader.vert` - Vertex shader with MVP uniform buffer
- `shader.frag` - Fragment shader for color output

## Expected Output

When shaders are correctly compiled and placed, `vulkan_cube.exe` should display a rotating 3D cube with colored faces in a window.
