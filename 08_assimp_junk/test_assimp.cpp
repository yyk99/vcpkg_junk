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

#include "assimp_aux.h"

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
        return test_data_dir / relative_path;
    }

    bool has_failure() const { return ::testing::Test::HasFailure(); }

    bool dump_as_assxml(const char *test_file, const char *out_prefix,
                        fs::path const &ws)
    {
        auto filename_gltf = test_data(test_file);
        if (!fs::is_regular_file(filename_gltf))
            return false;

        Assimp::Importer sot;
        {
            aiScene const *actual_scene =
                sot.ReadFile(filename_gltf.string().c_str(), 0);
            if(!actual_scene)
                return false;

            Assimp::Exporter exporter;
            auto rc = exporter.Export(actual_scene, "assxml",
                (ws / (std::string(out_prefix) + ".xml")).string());
            return rc == AI_SUCCESS;
        }
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

/// @brief Compare different ways of texture definitions
/// @param --gtest_filter=AssimpF.load_textured_cube_variants
/// @param  
TEST_F(AssimpF, load_textured_cube_variants) {
    auto ws = create_ws();

    // clang-format off
    const char *files [] = { 
        "BoxTextured-glTF/BoxTextured.gltf",
        "BoxTextured-glTF-Binary/BoxTextured.glb",
        "BoxTextured-glTF-Embedded/BoxTextured.gltf",
        };
    // clang-format on

    for (auto fp : files) {
        auto prefix = fs ::path(fp).parent_path();
        EXPECT_TRUE(dump_as_assxml(fp, prefix.string().c_str(), ws)) << "File: " << fp;
    }
}


//
//
//

#include <json-c/json.h>

/// @brief a quick json-c sample
/// @param --gtest_filter=AssimpF.json_c_test
TEST_F(AssimpF, json_c_test) {
    json_object *jobj = nullptr;
    const char *str = R"({ "msg-type": [ "0xdeadbeef", "irc log" ],
        "msg-from": { "class": "soldier", "name": "Wixilav" },
        "msg-to": { "class": "supreme-commander", "name": "[Redacted]" },
        "msg-log": [
            "soldier: Boss there is a slight problem with the piece offering to humans",
            "supreme-commander: Explain yourself soldier!",
            "soldier: Well they don't seem to move anymore...",
            "supreme-commander: Oh snap, I came here to see them twerk!"
            ]
        })";
    CONSOLE_EVAL(str);

    jobj = json_tokener_parse(str);
    ASSERT_TRUE(jobj);

    const char *pretty_text = json_object_to_json_string_ext(
        jobj, JSON_C_TO_STRING_SPACED | JSON_C_TO_STRING_PRETTY);
    ASSERT_TRUE(pretty_text);
    CONSOLE_EVAL(pretty_text);

    ASSERT_EQ(1, json_object_put(jobj));
}

TEST_F(AssimpF, json_c_test_array) {

    const char *json_string =
        R"({"name": "John Doe", "age": 30, "hobbies": ["reading", "hiking", "coding"]})";

    json_object *json_obj = json_tokener_parse(json_string);
    ASSERT_TRUE(json_obj) << "Error parsing JSON string";

    json_object *hobbies_array = json_object_object_get(json_obj, "hobbies");
    ASSERT_TRUE(hobbies_array);
    EXPECT_TRUE(json_object_get_type(hobbies_array) == json_type_array)
        << "Error: 'hobbies' is not an array or does not exist";
    if (has_failure()) {
        (void)json_object_put(json_obj); // Clean up json_obj before exiting
        GTEST_FAIL() << "See above for the reason";
    }
    int array_len = json_object_array_length(hobbies_array);
    ASSERT_EQ(3, array_len);
    {
        auto actual = json_object_array_get_idx(hobbies_array, 0);
        EXPECT_STREQ("reading", json_object_get_string(actual));
    }
    {
        auto actual = json_object_array_get_idx(hobbies_array, 1);
        EXPECT_STREQ("hiking", json_object_get_string(actual));
    }
    {
        auto actual = json_object_array_get_idx(hobbies_array, 2);
        EXPECT_STREQ("coding", json_object_get_string(actual));
    }

    json_object_put(json_obj);
}

#include "TilesetJson.h"

/// @brief One more json-c test
/// @param --gtest_filter=AssimpF.json_c_test_2
/// @param
TEST_F(AssimpF, json_c_test_2) {
    auto json_file = test_data("GM13206_test.tb_1_tiles/tileset.json");
    ASSERT_TRUE(fs::is_regular_file(json_file));

    cesiumjs::TilesetJson sot(json_file.string());
    {
        auto actual = sot.asset();
        ASSERT_FALSE(actual.extra_ion_georeferenced.has_value());
        EXPECT_EQ(std::string("1.0"), actual.version);
    }
    {
        auto root = sot.root();
        ASSERT_EQ(std::string("ADD"), root.refine);
    }
}

/// @brief A failure test
/// @param --gtest_filter=AssimpF.json_c_test_3
TEST_F(AssimpF, json_c_test_3) {
    EXPECT_THROW(cesiumjs::TilesetJson sot("/dev/null"), std::runtime_error);
}

#include "meshtoolbox.h"

/// @brief 
/// @param --gtest_filter=AssimpF.meshtoolbox_t0
/// @param  
TEST_F(AssimpF, meshtoolbox_t0) {

    auto ws = create_ws();

    meshtoolbox::Toolbox tb;

    meshtoolbox::box_t b0{{0, 0, 0}, {10, 20, 30}};

    std::vector<meshtoolbox::box_t> boxes{b0};

    auto model = std::unique_ptr<aiScene>(tb.make_boxes(boxes));

    ASSERT_TRUE(model);
    EXPECT_EQ(1, model->mNumMeshes);

    Assimp::Exporter exp;
    auto flags = aiProcess_GenNormals | aiProcess_ValidateDataStructure | 0;
    {
        std::string filename_glb = (ws / "meshtoolbox_t0.glft").string();
        auto err = exp.Export(model.get(), "gltf", filename_glb, flags);
        EXPECT_EQ(AI_SUCCESS, err) << "Failed export to " << filename_glb;
    }
    {
        std::string filename_glb = (ws / "meshtoolbox_t0.glb").string();
        auto err = exp.Export(model.get(), "glb2", filename_glb, flags);
        EXPECT_EQ(AI_SUCCESS, err) << "Failed export to " << filename_glb;
    }
    {
        std::string filename_glb = (ws / "meshtoolbox_t0.xml").string();
        auto err = exp.Export(model.get(), "assxml", filename_glb, flags);
        EXPECT_EQ(AI_SUCCESS, err) << "Failed export to " << filename_glb;
    }
}