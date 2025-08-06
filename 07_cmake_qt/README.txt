Preconditions

- install vcpkg
- install boots required packages...

---------------- cut here ----------------------
Total install time: 4.8 min
The package boost-program-options is compatible with built-in CMake targets of FindBoost.cmake:

    find_package(Boost REQUIRED COMPONENTS program_options)
    target_link_libraries(main PRIVATE Boost::program_options)

or the generated cmake configs via:

    find_package(boost_program_options REQUIRED CONFIG)
    target_link_libraries(main PRIVATE Boost::program_options)
---------------- cut here ----------------------


How to configure
	cmake -S . -B build-vs2022-x64 -DCMAKE_TOOLCHAIN_FILE=E:/opt/vcpkg/scripts/buildsystems/vcpkg.cmake 

	or

	cmake -S . -B build-linux -DCMAKE_TOOLCHAIN_FILE=/home/yyk/src/vcpkg/scripts/buildsystems/vcpkg.cmake


to be continue...
