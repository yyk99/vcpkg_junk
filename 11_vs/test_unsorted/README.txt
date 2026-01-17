## How to configure

Using CMake presets (recommended):

    cmake --preset=kestrel         # Linux (GCC)
    cmake --preset=kestrel-clang   # Linux (Clang)
    cmake --preset=default         # Windows (G:\opt\vcpkg)
    cmake --preset=cherry          # Windows (E:\opt\vcpkg)

Build:

    cmake --build --preset=kestrel

