#
#
#
cmake_minimum_required(VERSION 3.16)

project(07_cmake_qt VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Qt6 REQUIRED COMPONENTS Core Widgets Gui)
qt_standard_project_setup()

qt_add_executable(07_helloworld
    main.cpp
)

target_link_libraries(07_helloworld PRIVATE Qt6::Core Qt6::Widgets Qt6::Gui)
set_target_properties(07_helloworld PROPERTIES FOLDER "Apps")