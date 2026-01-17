cmake_minimum_required(VERSION 3.10)
project(14_vulkan_primer)

set(CMAKE_CXX_STANDARD 17)

find_package(Vulkan REQUIRED)
find_package(glfw3 REQUIRED)

# Find glslc shader compiler
find_program(GLSLC_EXECUTABLE
    NAMES glslc
    HINTS
        ${Vulkan_GLSLC_EXECUTABLE}
        "$ENV{VULKAN_SDK}/Bin"
        "$ENV{VULKAN_SDK}/bin"
        "/usr/bin"
        "/usr/local/bin"
)

if(NOT GLSLC_EXECUTABLE)
    message(WARNING "glslc not found - shaders will not be automatically compiled")
else()
    message(STATUS "Found glslc: ${GLSLC_EXECUTABLE}")

    # Compile vertex shader
    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/vert.spv
        COMMAND ${GLSLC_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/shader.vert -o ${CMAKE_CURRENT_BINARY_DIR}/vert.spv
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/shader.vert
        COMMENT "Compiling vertex shader"
        VERBATIM
    )

    # Compile fragment shader
    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/frag.spv
        COMMAND ${GLSLC_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/shader.frag -o ${CMAKE_CURRENT_BINARY_DIR}/frag.spv
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/shader.frag
        COMMENT "Compiling fragment shader"
        VERBATIM
    )

    # Create a custom target for all shaders
    add_custom_target(shaders ALL
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/vert.spv ${CMAKE_CURRENT_BINARY_DIR}/frag.spv
        COMMENT "Building shaders"
    )
endif()

add_executable(vulkan_app vulkan_app.cpp)
target_link_libraries(vulkan_app Vulkan::Vulkan glfw)

add_executable(vulkan_triangle vulkan_triangle.cpp)
target_link_libraries(vulkan_triangle Vulkan::Vulkan glfw)

find_package(glm REQUIRED)

add_executable(vulkan_cube vulkan_cube.cpp)
target_link_libraries(vulkan_cube Vulkan::Vulkan glfw glm::glm)
