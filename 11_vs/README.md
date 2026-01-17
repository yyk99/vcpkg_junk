# Visual Studio with vcpkg in the manifest mode

Apparently

	cd cmd_fmt_example
	vcpkg integrate install

## Build with specific preset

	cmake --build --preset default-debug
	cmake --build --preset cherry-release
	cmake --build --preset kestrel-debug

## CHANGELOG

### Branch devel/yurik42/11_vs_unicode_filenames

	Find out how to use UNICODE filepath in non-UNICODE Windows applications

### Branch devel/yurik42/11_vs_sombrero_elevation_generator

	Add cmd_sombrero_elevation_generator app




