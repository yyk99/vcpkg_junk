# Vulkan Primer

Three example Vulkan applications demonstrating graphics programming with Vulkan API, GLFW, and GLM.

## Applications

- **vulkan_app** - Basic Vulkan window initialization
- **vulkan_triangle** - Classic "Hello Triangle" with vertex colors (red, green, blue)
- **vulkan_cube** - 3D rotating colored cube with depth buffering

## Prerequisites

- vcpkg installed
- Vulkan SDK
- CMake 3.10+
- C++17 compiler

Ubuntu Linux required more dependencies and tools:

       $ sudo apt install libxinerama-dev libxcursor-dev xorg-dev libglu1-mesa-dev
       $ sudo apt-get install -y glslang-tools


## Building

### Using CMake Presets (Recommended)

Configure:

    cmake --preset=kestrel         # Linux (GCC)
    cmake --preset=kestrel-clang   # Linux (Clang)
    cmake --preset=default         # Windows (G:\opt\vcpkg)
    cmake --preset=cherry          # Windows (E:\opt\vcpkg)

Build:

    cmake --build --preset=kestrel

### Manual Configuration (Alternative)

Windows:

    cmake -S . -B build-vs2022-x64 -DCMAKE_TOOLCHAIN_FILE=c:/opt/vcpkg/scripts/buildsystems/vcpkg.cmake

Linux:

    cmake -S . -B build-linux -DCMAKE_TOOLCHAIN_FILE=/home/yyk/src/vcpkg/scripts/buildsystems/vcpkg.cmake

Build:

    cmake --build build-linux

## Running

    ./build-linux/vulkan_app
    ./build-linux/vulkan_triangle
    ./build-linux/vulkan_cube
