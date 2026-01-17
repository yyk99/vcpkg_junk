#define TINYGLTF_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "tiny_gltf.h"

#include <filesystem>
#include <iostream>
#include <string>

int main() {

    auto test_data = [](const char *relative_path) {
        auto test_data_dir =
            std::filesystem::absolute(__FILE__).parent_path() / "test_data";
        return test_data_dir / relative_path;
    };

    auto input_filename = test_data("tileset2/bounding_boxes.glb").string();

    tinygltf::Model model;
    tinygltf::TinyGLTF loader;
    std::string err;
    std::string warn;

#if 0
    bool ret = loader.LoadASCIIFromFile(&model, &err, &warn, input_filename);
#else
    // Alternatively, for binary glTF (.glb) files:
    bool ret = loader.LoadBinaryFromFile(&model, &err, &warn, input_filename);
#endif

    if (!warn.empty()) {
        std::cout << "WARN: " << warn << std::endl;
    }

    if (!err.empty()) {
        std::cerr << "ERR: " << err << std::endl;
    }

    if (!ret) {
        std::cerr << "Failed to load glTF model." << std::endl;
        return 1;
    }

    std::cout << "Successfully loaded glTF model: " << input_filename
              << std::endl;

    // You can now access the loaded model data, e.g.:
    std::cout << "Number of meshes: " << model.meshes.size() << std::endl;
    std::cout << "Number of nodes: " << model.nodes.size() << std::endl;

    // Further processing and rendering logic would go here,
    // involving iterating through meshes, accessors, bufferViews, etc.
    // and setting up your graphics API to render the geometry.

    return 0;
}