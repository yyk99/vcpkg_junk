cmake_minimum_required(VERSION 3.15)
project(cmd_fmt_example)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add your source files here
add_executable(cmd_fmt_example
    main.cpp)

#The package fmt provides CMake targets:
#
#    find_package(fmt CONFIG REQUIRED)
#    target_link_libraries(main PRIVATE fmt::fmt)
#
#    # Or use the header-only version
#    find_package(fmt CONFIG REQUIRED)
#    target_link_libraries(main PRIVATE fmt::fmt-header-only)

# Find Boost (adjust components as needed)
find_package(fmt REQUIRED)

target_link_libraries(cmd_fmt_example PRIVATE fmt::fmt)
