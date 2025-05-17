//
//
//

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string>

namespace fs = std::filesystem;

#include <assimp/Exporter.hpp>
#include <assimp/Importer.hpp>
#include <assimp/scene.h>

#include "../common/DebuggingConsole.h"

/// @brief
class AssimpF : public testing::Test {
protected:
    std::string test_name() const {
        return ::testing::UnitTest::GetInstance()->current_test_info()->name();
    }

    std::string test_case_name() const {
        return ::testing::UnitTest::GetInstance()
            ->current_test_info()
            ->test_case_name();
    }

    fs::path create_ws() {
        auto ws = fs::absolute("out") / test_case_name() / test_name();
        fs::remove_all(ws);
        fs::create_directories(ws);
        return ws;
    }

    void save_as(std::string const &text, fs::path const &location) {
        std::ofstream os(location.string().c_str());
        if (!os.good())
            throw std::runtime_error("Cannot create file: " +
                                     location.string());
        os << text;
        os.close();
    }

    auto test_data(const char *relative_path) {
        auto test_data_dir = fs::absolute(__FILE__).parent_path() / "test_data";
        CONSOLE_EVAL(test_data_dir);

        return test_data_dir / relative_path;
    }
};

/// @brief
/// @param
/// @param
TEST_F(AssimpF, toolkit_create_ws) {
    auto ws = create_ws();

    ASSERT_TRUE(fs::is_directory(ws));
    save_as("hello\n", ws / "foo.txt");
    ASSERT_TRUE(fs::is_regular_file(ws / "foo.txt"));

    auto ws2 = create_ws();

    ASSERT_EQ(ws, ws2);

    ASSERT_TRUE(fs::is_directory(ws2));
    ASSERT_FALSE(fs::is_regular_file(ws2 / "foo.txt"));
    save_as("hello2\n", ws2 / "foo2.txt");
    ASSERT_TRUE(fs::is_regular_file(ws2 / "foo2.txt"));
}

#include <assimp/DefaultLogger.hpp>
#include <assimp/Logger.hpp>

/// @brief
/// @param
/// @param
TEST_F(AssimpF, z1) {
    auto ws = create_ws();

    auto cube_ply =
        (fs::path(__FILE__).parent_path() / "test_data/cube.ply").string();

    std::cout << cube_ply << "\n";

    Assimp::DefaultLogger::create((ws / "logfile.txt").string().c_str(),
                                  Assimp::Logger::VERBOSE);
    ASSIMP_LOG_DEBUG("First message");

    Assimp::Importer sot;
    {
        auto actual_scene = sot.ReadFile("simple.ply", 0);
        ASSERT_EQ(nullptr, actual_scene);
    }

    {
        auto actual_scene = sot.ReadFile(cube_ply.c_str(), 0);
        ASSERT_NE(nullptr, actual_scene);

        Assimp::Exporter exporter;
        auto rc = exporter.Export(actual_scene, "assxml",
                                  (ws / "cube_ply.xml").string());
        ASSERT_EQ(0, rc);
    }
}

std::ostream &operator<<(std::ostream &ss, aiVector3D const &v) {
    ss << "{" << v.x << "," << v.y << "," << v.z << "}";
    return ss;
}

/// @brief Load a textured cube and save it as assxml
/// @param
/// @param
TEST_F(AssimpF, load_textured_cube) {
    auto ws = create_ws();

    // test_data/BoxTextured-glTF/BoxTextured.gltf

    auto BoxTextured_gltf = test_data("BoxTextured-glTF/BoxTextured.gltf");
    ASSERT_TRUE(fs::is_regular_file(BoxTextured_gltf));

    Assimp::Importer sot;
    {
        aiScene const *actual_scene =
            sot.ReadFile(BoxTextured_gltf.string().c_str(), 0);
        ASSERT_NE(nullptr, actual_scene);

        CONSOLE_EVAL(actual_scene->mNumMeshes);
        CONSOLE_EVAL(actual_scene->mNumTextures);
        CONSOLE_EVAL(actual_scene->mNumMaterials);
        CONSOLE_EVAL(actual_scene->mNumSkeletons);
        // CONSOLE_EVAL(actual_scene->mNumTextures);

        ASSERT_EQ(1, actual_scene->mNumMeshes);
        {
            // dump the mesh
            aiMesh const *mp = actual_scene->mMeshes[0];
            CONSOLE_EVAL(mp->mNumUVComponents[0]);
            CONSOLE_EVAL(mp->mNumUVComponents[1]);
            CONSOLE_EVAL(mp->mNumUVComponents[2]);
            CONSOLE_EVAL(mp->mNumUVComponents[3]);

            aiVector3D *textCoord_0 = mp->mTextureCoords[0];
            CONSOLE_EVAL(textCoord_0[0]);
            CONSOLE_EVAL(textCoord_0[1]);
            CONSOLE_EVAL(textCoord_0[2]);
            CONSOLE_EVAL(textCoord_0[3]);
            CONSOLE_EVAL(textCoord_0[4]);
        }

        Assimp::Exporter exporter;
        auto rc = exporter.Export(actual_scene, "assxml",
                                  (ws / "BoxTextured.gltf.xml").string());
        ASSERT_EQ(0, rc);
    }
}

//
// 
//

#include <json-c/json.h>

#include "TilesetJson.h"

TEST_F(AssimpF, json_c_test)
{

}