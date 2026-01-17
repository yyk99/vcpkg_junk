//
//
//

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string>

namespace fs = std::filesystem;

#include <assimp/DefaultLogger.hpp>
#include <assimp/Exporter.hpp>
#include <assimp/Importer.hpp>
#include <assimp/Logger.hpp>
#include <assimp/scene.h>

#define STB_IMAGE_IMPLEMENTATION

#include <stb_image.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#include "../common/CONSOLE.h"
#include "assimp_aux.h"

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
        auto test_data_dir =
            fs::absolute(__FILE__).parent_path().parent_path() /
            "08_assimp_junk" / "test_data";
        return test_data_dir / relative_path;
    }

    bool has_failure() const { return ::testing::Test::HasFailure(); }

    bool dump_as_assxml(const char *test_file, const char *out_prefix,
                        fs::path const &ws) {
        auto filename_gltf = test_data(test_file);
        if (!fs::is_regular_file(filename_gltf))
            return false;

        Assimp::Importer sot;
        {
            aiScene const *actual_scene =
                sot.ReadFile(filename_gltf.string().c_str(), 0);
            if (!actual_scene)
                return false;

            Assimp::Exporter exporter;
            auto rc = exporter.Export(
                actual_scene, "assxml",
                (ws / (std::string(out_prefix) + ".xml")).string());
            return rc == AI_SUCCESS;
        }
    }

    std::vector<char> read_file(const char *filename) {
        std::ifstream file(filename, std::ios::binary | std::ios::ate);
        if (!file)
            throw std::runtime_error(std::string("Cannot open file: ") +
                                     filename);
        std::streamsize size = file.tellg();
        file.seekg(0, std::ios::beg);
        std::vector<char> buffer(size);
        if (!file.read(buffer.data(), size))
            throw std::runtime_error(std::string("Failed to read file: ") +
                                     filename);
        return buffer;
    }

    /**
     * @brief Finds the longest common substring among all strings in the given
     * vector.
     *
     * This function searches for the longest substring that is present in every
     * string in the input vector. If multiple substrings of the same maximum
     * length exist, the first one found in the first string is returned. If the
     * input vector is empty, an empty string is returned.
     *
     * @param s A vector of strings to search for the longest common substring.
     * @return The longest common substring found in all strings, or an empty
     * string if none exists.
     */
    std::string get_longest_substring(std::vector<std::string> const &s) {
        if (s.empty())
            return "";

        const std::string &first = s[0];
        for (size_t len = first.size(); len > 0; --len) {
            for (size_t start = 0; start + len <= first.size(); ++start) {
                std::string_view candidate(first.data() + start, len);
                bool found = true;
                for (size_t i = 1; i < s.size(); ++i) {
                    if (s[i].find(candidate) == std::string::npos) {
                        found = false;
                        break;
                    }
                }
                if (found)
                    return std::string(candidate);
            }
        }
        return "";
    }

    std::string get_longest_substring(std::string const &s1,
                                      std::string const &s2) {
        if (s1.empty() || s2.empty())
            return "";
        size_t max_len = 0;
        size_t max_pos = 0;
        for (size_t len = s1.size(); len > 0; --len) {
            for (size_t start = 0; start + len <= s1.size(); ++start) {
                std::string_view candidate(s1.data() + start, len);
                if (s2.find(candidate) != std::string::npos) {
                    // Return the first found (longest, leftmost)
                    return std::string(candidate);
                }
            }
        }
        return "";
    }

    /// @brief Compute the bounding volume of a scene with respect to node
    /// transforms
    /// @param scene
    /// @return aiAABB (axis-aligned bounding box in world space)
    aiAABB compute_aabb_with_transform(const aiScene *scene) {
        aiVector3D sceneMin(std::numeric_limits<float>::max(),
                            std::numeric_limits<float>::max(),
                            std::numeric_limits<float>::max());
        aiVector3D sceneMax(std::numeric_limits<float>::lowest(),
                            std::numeric_limits<float>::lowest(),
                            std::numeric_limits<float>::lowest());

        std::function<void(const aiNode *, const aiMatrix4x4 &)> traverse;
        traverse = [&](const aiNode *node, const aiMatrix4x4 &parentTransform) {
            aiMatrix4x4 transform = parentTransform * node->mTransformation;
            for (unsigned int i = 0; i < node->mNumMeshes; ++i) {
                const aiMesh *mesh = scene->mMeshes[node->mMeshes[i]];
                for (unsigned int v = 0; v < mesh->mNumVertices; ++v) {
                    aiVector3D vertex = mesh->mVertices[v];
                    vertex *= transform;
                    sceneMin.x = std::min(sceneMin.x, vertex.x);
                    sceneMin.y = std::min(sceneMin.y, vertex.y);
                    sceneMin.z = std::min(sceneMin.z, vertex.z);
                    sceneMax.x = std::max(sceneMax.x, vertex.x);
                    sceneMax.y = std::max(sceneMax.y, vertex.y);
                    sceneMax.z = std::max(sceneMax.z, vertex.z);
                }
            }
            for (unsigned int i = 0; i < node->mNumChildren; ++i) {
                traverse(node->mChildren[i], transform);
            }
        };

        if (scene && scene->mRootNode) {
            traverse(scene->mRootNode, aiMatrix4x4());
        }
        return aiAABB{sceneMin, sceneMax};
    }

    /// @brief Compute the bounding volume of a scene w/o respect to transforms
    /// @param actual
    /// @return
    aiAABB compute_aabb(aiScene const *actual) {
        aiVector3D sceneMin(std::numeric_limits<float>::max(),
                            std::numeric_limits<float>::max(),
                            std::numeric_limits<float>::max());
        aiVector3D sceneMax(std::numeric_limits<float>::lowest(),
                            std::numeric_limits<float>::lowest(),
                            std::numeric_limits<float>::lowest());

        for (unsigned int i = 0; i < actual->mNumMeshes; ++i) {
            const aiMesh *mesh = actual->mMeshes[i];
            for (unsigned int v = 0; v < mesh->mNumVertices; ++v) {
                const aiVector3D &vertex = mesh->mVertices[v];
                sceneMin.x = std::min(sceneMin.x, vertex.x);
                sceneMin.y = std::min(sceneMin.y, vertex.y);
                sceneMin.z = std::min(sceneMin.z, vertex.z);
                sceneMax.x = std::max(sceneMax.x, vertex.x);
                sceneMax.y = std::max(sceneMax.y, vertex.y);
                sceneMax.z = std::max(sceneMax.z, vertex.z);
            }
        }
        // CONSOLE("Scene AABB min: " << sceneMin << " max: " << sceneMax);
        return aiAABB{sceneMin, sceneMax};
    };
};

#if 0
TEST_F(AssimpF, get_longest_substring_test) {
    {
        std::vector<std::string> v{};
        EXPECT_EQ("", get_longest_substring(v));
    }
    {
        std::vector<std::string> v{
            "export2x2_A0",
        };
        EXPECT_EQ("export2x2_A0", get_longest_substring(v));
    }
    {
        std::vector<std::string> v{
            "export2x2_A0",
            "export2x2_A1",
        };
        EXPECT_EQ("export2x2_A", get_longest_substring(v));
    }
    {
        std::vector<std::string> v{
            "export2x2_A0",
            "export2x2_A1",
            "export2x2_B0",
            "export2x2_B1",
        };
        EXPECT_EQ("export2x2_", get_longest_substring(v));
    }
}

TEST_F(AssimpF, get_longest_substring_test_2) {
    EXPECT_EQ("", get_longest_substring("", ""));
    EXPECT_EQ("export2x2_A0",
              get_longest_substring("export2x2_A0", "export2x2_A0"));
    EXPECT_EQ("export2x2_A",
              get_longest_substring("export2x2_A0", "export2x2_A1"));
}

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

/// @brief Load a textured (jpg) cube and save it as assxml
/// @param --gtest_filter=AssimpF.load_textured_cube_jpg
/// @param
TEST_F(AssimpF, load_textured_cube_jpg) {
    auto ws = create_ws();

    auto BoxTextured_gltf = test_data("BoxTextured-glTF_jpg/BoxTextured.gltf");
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

/// @brief Load a textured cube and save it as assxml
/// @param --gtest_filter=AssimpF.load_textured_cube_embedded
TEST_F(AssimpF, load_textured_cube_embedded) {

    auto BoxTextured_gltf =
        test_data("BoxTextured-glTF-Embedded/BoxTextured.gltf");
    ASSERT_TRUE(fs::is_regular_file(BoxTextured_gltf));

    auto ws = create_ws();

    Assimp::Importer sot;
    {
        aiScene const *actual_scene =
            sot.ReadFile(BoxTextured_gltf.string().c_str(), 0);
        ASSERT_TRUE(actual_scene);
        ASSERT_EQ(1, actual_scene->mNumMeshes);
        ASSERT_EQ(1, actual_scene->mNumTextures);
        ASSERT_EQ(2, actual_scene->mNumMaterials);
        ASSERT_EQ(0, actual_scene->mNumSkeletons);
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

        {
            // dump the texture
            CONSOLE_EVAL(actual_scene->mTextures[0]->achFormatHint);
            CONSOLE_EVAL(actual_scene->mTextures[0]->mFilename);
            CONSOLE_EVAL(actual_scene->mTextures[0]->mHeight);
            CONSOLE_EVAL(actual_scene->mTextures[0]->mWidth);
        }
        {
            // dump the materials
            for (int i = 0; i != actual_scene->mNumMaterials; ++i) {
                CONSOLE("========= i = " << i << " =========");
                CONSOLE_EVAL(actual_scene->mMaterials[i]->GetName());
                CONSOLE_EVAL(actual_scene->mMaterials[i]->mNumProperties);
                auto props = actual_scene->mMaterials[i]->mProperties;
                auto n_props = actual_scene->mMaterials[i]->mNumProperties;
                for (int j = 0; j != n_props; ++j) {
                    CONSOLE_EVAL(props[j]->mKey);
                }
            }
        }

        Assimp::Exporter exporter;
        auto rc = exporter.Export(actual_scene, "assxml",
                                  (ws / "model.xml").string());
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
        EXPECT_TRUE(dump_as_assxml(fp, prefix.string().c_str(), ws))
            << "File: " << fp;
    }
}

/// @brief Load a mesh with external texture and save embedded
/// @param --gtest_filter=AssimpF.save_texture_embedded
/// @param
TEST_F(AssimpF, save_texture_embedded) {
    const char *test_file = "BoxTextured-glTF/BoxTextured.gltf";

    auto filename_gltf = test_data(test_file);
    if (!fs::is_regular_file(filename_gltf))
        GTEST_SKIP();

    auto ws = create_ws();

    Assimp::Importer sot;

    aiScene const *actual_scene =
        sot.ReadFile(filename_gltf.string().c_str(), 0);
    ASSERT_TRUE(actual_scene);

    Assimp::Exporter exporter;
    auto rc =
        exporter.Export(actual_scene, "gltf2", (ws / "model.glb").string(),
                        aiProcess_EmbedTextures);
    ASSERT_EQ(AI_SUCCESS, rc);
    // TODO: verify the export is valid
    //       apparently it is not
}

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
    auto flags =
        /*aiProcess_GenNormals | */ aiProcess_ValidateDataStructure | 0;
    {
        std::string filename_glb = (ws / "model.glft").string();
        auto err = exp.Export(model.get(), "gltf", filename_glb, flags);
        EXPECT_EQ(AI_SUCCESS, err) << "Failed export to " << filename_glb;
    }
    {
        std::string filename_glb = (ws / "model.glb").string();
        auto err = exp.Export(model.get(), "glb2", filename_glb, flags);
        EXPECT_EQ(AI_SUCCESS, err) << "Failed export to " << filename_glb;
    }
    {
        std::string filename_glb = (ws / "model.xml").string();
        auto err = exp.Export(model.get(), "assxml", filename_glb, flags);
        EXPECT_EQ(AI_SUCCESS, err) << "Failed export to " << filename_glb;
    }
}

/// @brief Create a georeferenced "graveyard" scene (the test data)
/// @param
/// --gtest_filter=AssimpF.meshtoolbox_t1
///
TEST_F(AssimpF, meshtoolbox_t1) {

    auto ws = create_ws();

    meshtoolbox::Toolbox tb;

    // A "graveyard"
    ai_real max_X = 2000;
    ai_real max_Z = 1000;
    ai_real step = 50; // 50m interval between "obelisks"

    std::vector<meshtoolbox::box_t> boxes;
    boxes.reserve(size_t(max_X / step + 1) * size_t(max_Z / step + 1));
    for (int i = 0; i < max_X / step; ++i) {
        for (int j = 0; j < max_Z / step; ++j) {
            boxes.emplace_back(meshtoolbox::box_t{{i * step, 0, j * step},
                                                  {step / 5, step, 4}});
        }
    }

    boxes.emplace_back(
        meshtoolbox::box_t{{max_X, 0, max_Z}, {step, step * 10, step}});
    boxes.emplace_back(
        meshtoolbox::box_t{{max_X, 0, 0}, {step, step * 10, step}});
    boxes.emplace_back(
        meshtoolbox::box_t{{0, 0, max_Z}, {step, step * 10, step}});

    auto actual = std::unique_ptr<aiScene>(tb.make_boxes(boxes));
    ASSERT_TRUE((bool)actual);

    {
        Assimp::Exporter exp;
        exp.Export(actual.get(), "glb2", (ws / "actual.glb").string());
        exp.Export(actual.get(), "assxml", (ws / "actual.xml").string());
    }
}

/// @brief Create a GLB model w/o normals
TEST_F(AssimpF, meshtoolbox_wo_normals) {

    auto ws = create_ws();

    meshtoolbox::Toolbox tb;

    // Real world coordinates, Z-Up
    meshtoolbox::box_t b0{{0, 0, 0}, {10, 20, 30}};

    std::vector<meshtoolbox::box_t> boxes{b0};

    auto model = std::unique_ptr<aiScene>(tb.make_boxes_yup(boxes));

    ASSERT_TRUE(model);
    EXPECT_EQ(1, model->mNumMeshes);

    Assimp::Exporter exp;
    auto flags = aiProcess_ValidateDataStructure | 0;
    {
        std::string filename_glb = (ws / "model.glb").string();
        auto err = exp.Export(model.get(), "glb2", filename_glb, flags);
        EXPECT_EQ(AI_SUCCESS, err) << "Failed export to " << filename_glb;
    }
    {
        std::string filename_glb = (ws / "model.xml").string();
        auto err = exp.Export(model.get(), "assxml", filename_glb, flags);
        EXPECT_EQ(AI_SUCCESS, err) << "Failed export to " << filename_glb;
    }
}

/// @brief
/// @param --gtest_filter=AssimpF.meshtoolbox_box_with_axis
/// @param
TEST_F(AssimpF, meshtoolbox_box_with_axis) {

    auto ws = create_ws();

    meshtoolbox::Toolbox tb;

    meshtoolbox::box_t b0{{0, 0, 0}, {10, 20, 30}};
    meshtoolbox::box_t b1{{11, 0, 0}, {1, 1, 1}};
    meshtoolbox::box_t b2{{0, 21, 0}, {1, 1, 1}};
    meshtoolbox::box_t b3{{0, 0, 31}, {1, 1, 1}};

    std::vector<meshtoolbox::box_t> boxes{b0, b1, b2, b3};

    auto model = std::unique_ptr<aiScene>(tb.make_boxes(boxes));

    ASSERT_TRUE(model);
    EXPECT_EQ(4, model->mNumMeshes);

    Assimp::Exporter exp;
    auto flags = aiProcess_ValidateDataStructure | 0;
    {
        std::string filename_glb = (ws / "model.glb").string();
        auto err = exp.Export(model.get(), "glb2", filename_glb, flags);
        EXPECT_EQ(AI_SUCCESS, err) << "Failed export to " << filename_glb;
    }
    {
        std::string filename_glb = (ws / "model.xml").string();
        auto err = exp.Export(model.get(), "assxml", filename_glb, flags);
        EXPECT_EQ(AI_SUCCESS, err) << "Failed export to " << filename_glb;
    }
}

/// @brief Use stb_image.h
/// @param --gtest_filter=AssimpF.meshtoolbox_stb_image
/// @param
TEST_F(AssimpF, meshtoolbox_stb_image) {
    auto logo_png = test_data("BoxTextured-glTF/CesiumLogoFlat.png");
    ASSERT_TRUE(fs::is_regular_file(logo_png));

    int width, height, channels;
    unsigned char *data =
        stbi_load(logo_png.string().c_str(), &width, &height, &channels, 0);
    ASSERT_TRUE(data != nullptr) << "Failed to load image: " << logo_png;

    ASSERT_EQ(211, width);
    ASSERT_EQ(211, height);
    ASSERT_EQ(3, channels);

    // Optionally, check a pixel value
    if (width > 0 && height > 0 && channels >= 3) {
        int idx = 0; // top-left pixel
        unsigned char r = data[idx * channels + 0];
        unsigned char g = data[idx * channels + 1];
        unsigned char b = data[idx * channels + 2];
        CONSOLE_EVAL(unsigned(r));
        CONSOLE_EVAL(unsigned(g));
        CONSOLE_EVAL(unsigned(b));
    }

    // Free the image memory
    stbi_image_free(data);
}

/// @brief Load a large JPG file
/// @param --gtest_filter=AssimpF.meshtoolbox_stb_image_jpg
/// @param
TEST_F(AssimpF, meshtoolbox_stb_image_jpg) {
    const char *filename_jpg =
        R"(C:\home\work\GM-19017\out\scene_0_dense_mesh_material_0_map_kd.jpg)";
    if (!fs::is_regular_file(filename_jpg))
        GTEST_SKIP();

    int width, height, channels;
    auto data = stbi_load(filename_jpg, &width, &height, &channels, 0);
    ASSERT_TRUE(data) << "Failed to load image: " << filename_jpg;
#if 0
    ASSERT_EQ(211, width);
    ASSERT_EQ(211, height);
    ASSERT_EQ(3, channels);

    // Optionally, check a pixel value
    if (width > 0 && height > 0 && channels >= 3) {
        int idx = 0; // top-left pixel
        unsigned char r = data[idx * channels + 0];
        unsigned char g = data[idx * channels + 1];
        unsigned char b = data[idx * channels + 2];
        CONSOLE_EVAL(unsigned(r));
        CONSOLE_EVAL(unsigned(g));
        CONSOLE_EVAL(unsigned(b));
    }
#endif
    // Free the image memory
    stbi_image_free(data);
}

/// @brief Read a small jpg file with stb library
/// @param --gtest_filter=AssimpF.stb_read_jpg
/// @param
TEST_F(AssimpF, stb_read_jpg) {
    auto logo_jpg = test_data("BoxTextured-glTF_jpg/CesiumLogoFlat.jpg");
    ASSERT_TRUE(fs::is_regular_file(logo_jpg));

    int width, height, channels;
    auto stbi_uc_deleter = [](stbi_uc *p) {
        if (p)
            stbi_image_free(p);
    };
    auto data = std::unique_ptr<stbi_uc, decltype(stbi_uc_deleter)>(
        stbi_load(logo_jpg.string().c_str(), &width, &height, &channels, 0),
        stbi_uc_deleter);
#if _WIN32
    ASSERT_TRUE(data) << "Failed: " << logo_jpg
                      << ", reason:" << stbi_failure_reason();

    ASSERT_EQ(211, width);
    ASSERT_EQ(211, height);
    ASSERT_EQ(3, channels);

    // Optionally, check a pixel value
    if (width > 0 && height > 0 && channels >= 3) {
        int idx = 0; // top-left pixel
        unsigned char r = data.get()[idx * channels + 0];
        unsigned char g = data.get()[idx * channels + 1];
        unsigned char b = data.get()[idx * channels + 2];
        CONSOLE_EVAL(unsigned(r));
        CONSOLE_EVAL(unsigned(g));
        CONSOLE_EVAL(unsigned(b));
    }
#else
    // It is expected that JPEG is not supported in Linux
    // (because of the assimp limitations)
    ASSERT_FALSE(data);
    const char *error = stbi_failure_reason();
    ASSERT_TRUE(strstr(error, "unknown") != nullptr) << "Error: " << error;
#endif
}

/// @brief Create a PNG file
/// @param --gtest_filter=AssimpF.meshtoolbox_write_stb_image
TEST_F(AssimpF, meshtoolbox_write_stb_image) {
    auto ws = create_ws();

    int width = 256;
    int height = 128;
    int channels = 3;
    std::vector<unsigned char> image(width * height * channels, 0);

    // Fill with a gradient
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            int idx = (y * width + x) * channels;
            image[idx + 0] =
                static_cast<unsigned char>(x * 255 / (width - 1)); // R
            image[idx + 1] =
                static_cast<unsigned char>(y * 255 / (height - 1)); // G
            image[idx + 2] = 128;                                   // B
        }
    }

    auto out_png = ws / "gradient.png";
    int rc = stbi_write_png(out_png.string().c_str(), width, height, channels,
                            image.data(), width * channels);
    ASSERT_TRUE(rc != 0) << "Failed to write PNG: " << out_png;
    ASSERT_TRUE(fs::is_regular_file(out_png));
}

/// @brief Create a scene with a tile and texture
/// @param --gtest_filter=AssimpF.meshtoolbox_tile0
TEST_F(AssimpF, meshtoolbox_tile0) {

    auto logo_png = test_data("BoxTextured-glTF/CesiumLogoFlat.png");
    ASSERT_TRUE(fs::is_regular_file(logo_png));

    auto ws = create_ws();

    meshtoolbox::Toolbox tb;

    auto tile0 = tb.mesh_tile(211, 211);

    auto model = std::unique_ptr<aiScene>(tb.make_scene({tile0}));

    ASSERT_TRUE(model);
    EXPECT_EQ(1, model->mNumMeshes);

    // Let's attach a "texture"
    {
        auto data = read_file(logo_png.string().c_str());
        ASSERT_EQ(2433, data.size());

        ASSERT_EQ(0, model->mNumMaterials);
        model->mNumMaterials = 1;
        model->mMaterials = new aiMaterial *[model->mNumMaterials];
        model->mMaterials[0] = new aiMaterial{};

        ASSERT_EQ(0, model->mNumTextures);
        model->mNumTextures = 1;
        model->mTextures = new aiTexture *[model->mNumTextures];
        model->mTextures[0] = new aiTexture{};
        auto *tex = model->mTextures[0];
        tex->mWidth = data.size();
        tex->mHeight = 0;
        tex->pcData = (aiTexel *)malloc(data.size());
        memcpy(tex->pcData, data.data(), data.size());
        strcpy(tex->achFormatHint, "png");

        // Set the material to use the texture
        aiString texPath;
        texPath.Set("*0"); // Embedded texture reference
        model->mMaterials[0]->AddProperty(&texPath,
                                          AI_MATKEY_TEXTURE_DIFFUSE(0));

        tile0->mMaterialIndex = 0;
        tile0->mNumUVComponents[0] = 2;
        tile0->mTextureCoords[0] = new aiVector3D[4];
        tile0->mTextureCoords[0][0] = {0, 0, 0};
        tile0->mTextureCoords[0][1] = {1, 0, 0};
        tile0->mTextureCoords[0][2] = {1, 1, 0};
        tile0->mTextureCoords[0][3] = {0, 1, 0};
    }

    Assimp::Exporter exp;
    auto flags = aiProcess_ValidateDataStructure | 0;
    {
        std::string filename = (ws / "model.glb").string();
        auto err = exp.Export(model.get(), "glb2", filename, flags);
        EXPECT_EQ(AI_SUCCESS, err) << "Failed export to " << filename;
        ASSERT_TRUE(fs::is_regular_file(filename)) << filename;
    }
    {
        std::string filename = (ws / "model.xml").string();
        auto err = exp.Export(model.get(), "assxml", filename, flags);
        EXPECT_EQ(AI_SUCCESS, err) << "Failed export to " << filename;
        ASSERT_TRUE(fs::is_regular_file(filename)) << filename;
    }

    // self-check
    {
        std::string filename = (ws / "model.glb").string();
        ASSERT_TRUE(fs::is_regular_file(filename));
        Assimp::Importer imp;
        aiScene const *model = imp.ReadFile(filename.c_str(), 0);
        ASSERT_TRUE(model);

        CONSOLE_EVAL(model->mNumTextures);
        ASSERT_EQ(1, model->mNumTextures);
    }
}

/// @brief Research DRACO compressing algo
/// @param --gtest_filter=AssimpF.meshtoolbox_tile1
TEST_F(AssimpF, meshtoolbox_tile1) {
    auto filename = test_data("draco/2CylinderEngine.gltf");
    ASSERT_TRUE(fs::is_regular_file(filename));

    auto ws = create_ws();

    Assimp::Importer sot;
    aiScene const *model = sot.ReadFile(filename.string().c_str(), 0);
    ASSERT_EQ(nullptr, model);
    EXPECT_STREQ("GLTF: Draco mesh compression not supported.",
                 sot.GetErrorString());
}

/// @brief Create a scene with a cone
/// @param --gtest_filter=AssimpF.meshtoolbox_cone0
///
TEST_F(AssimpF, meshtoolbox_cone0) {

    auto ws = create_ws();

    meshtoolbox::Toolbox tb;

    auto cone0 = tb.mesh_cone(30, 10, 6);

    auto model = std::unique_ptr<aiScene>(tb.make_scene({cone0}));

    ASSERT_TRUE(model);
    EXPECT_EQ(1, model->mNumMeshes);

    Assimp::Exporter exp;
    auto flags = aiProcess_ValidateDataStructure | 0;
    {
        std::string filename = (ws / "model.glb").string();
        auto err = exp.Export(model.get(), "glb2", filename, flags);
        EXPECT_EQ(AI_SUCCESS, err) << "Failed export to " << filename;
        ASSERT_TRUE(fs::is_regular_file(filename)) << filename;
    }
    {
        std::string filename = (ws / "model.xml").string();
        auto err = exp.Export(model.get(), "assxml", filename, flags);
        EXPECT_EQ(AI_SUCCESS, err) << "Failed export to " << filename;
        ASSERT_TRUE(fs::is_regular_file(filename)) << filename;
    }
}

/// @brief Create a scene with a cone and a cylinder
/// @param --gtest_filter=AssimpF.meshtoolbox_cone_and_cylinder
///
TEST_F(AssimpF, meshtoolbox_cone_and_cylinder) {

    auto ws = create_ws();

    meshtoolbox::Toolbox tb;

    auto cone0 = tb.mesh_cone(30, 10, 6);
    auto cyl0 = tb.mesh_cylinder(100, 1, 4, "cyl0");

    auto model = std::unique_ptr<aiScene>(tb.make_scene({cone0, cyl0}));

    for (auto const *np : tb.list_nodes(model.get())) {
        CONSOLE_EVAL(np->mName);
    }

    ASSERT_TRUE(model);
    EXPECT_EQ(2, model->mNumMeshes);

    auto *mp0 = model->mRootNode->FindNode("mesh_0");
    ASSERT_TRUE(mp0);

    auto *mp1 = model->mRootNode->FindNode("mesh_1");
    ASSERT_TRUE(mp1);

    tb.shift_down(mp1, 100);

    Assimp::Exporter exp;
    auto flags = aiProcess_ValidateDataStructure | 0;
    {
        std::string filename = (ws / "model.glb").string();
        auto err = exp.Export(model.get(), "glb2", filename, flags);
        EXPECT_EQ(AI_SUCCESS, err) << "Failed export to " << filename;
        ASSERT_TRUE(fs::is_regular_file(filename)) << filename;
    }
    {
        std::string filename = (ws / "model.xml").string();
        auto err = exp.Export(model.get(), "assxml", filename, flags);
        EXPECT_EQ(AI_SUCCESS, err) << "Failed export to " << filename;
        ASSERT_TRUE(fs::is_regular_file(filename)) << filename;
    }
}

/// @brief Create a scene with three coordinate axis
/// @param --gtest_filter=AssimpF.meshtoolbox_axis
TEST_F(AssimpF, meshtoolbox_axis) {

    auto ws = create_ws();

    meshtoolbox::Toolbox tb;

    auto cone0 = tb.mesh_cone(30, 10, 6);
    tb.paint_mesh(cone0, aiColor4D(0.0f, 1.0f, 0.0f, 1.0f)); // green
    auto cyl0 = tb.mesh_cylinder(100, 1, 4, "cyl0");
    auto cone1 = tb.mesh_cone(30, 10, 6);
    tb.paint_mesh(cone1, aiColor4D(0.0f, 0.0f, 1.0f, 1.0f)); // blue
    auto cyl1 = tb.mesh_cylinder(100, 1, 4, "cyl1");
    auto cone2 = tb.mesh_cone(30, 10, 6);
    tb.paint_mesh(cone2, aiColor4D(1.0f, 0.0f, 0.0f, 1.0f)); // red
    auto cyl2 = tb.mesh_cylinder(100, 1, 4, "cyl2");

    auto model = std::unique_ptr<aiScene>(
        tb.make_scene({cone0, cyl0, cone1, cyl1, cone2, cyl2}));

    for (auto const *np : tb.list_nodes(model.get())) {
        CONSOLE_EVAL(np->mName);
    }

    ASSERT_TRUE(model);
    EXPECT_EQ(6, model->mNumMeshes);

    auto *mp_con0 = model->mRootNode->FindNode("mesh_0");
    ASSERT_TRUE(mp_con0);
    auto *mp_cyl0 = model->mRootNode->FindNode("mesh_1");
    ASSERT_TRUE(mp_cyl0);

    auto *mp_con1 = model->mRootNode->FindNode("mesh_2");
    ASSERT_TRUE(mp_con0);
    auto *mp_cyl1 = model->mRootNode->FindNode("mesh_3");
    ASSERT_TRUE(mp_cyl0);

    auto *mp_con2 = model->mRootNode->FindNode("mesh_4");
    ASSERT_TRUE(mp_con0);
    auto *mp_cyl2 = model->mRootNode->FindNode("mesh_5");
    ASSERT_TRUE(mp_cyl0);

    tb.shift_up(mp_con0, 100);
    tb.shift_up(mp_con1, 100);
    tb.shift_up(mp_con2, 100);

    tb.turn_around_X(mp_con1, 90);
    tb.turn_around_X(mp_cyl1, 90);

    tb.turn_around_Z(mp_con2, -90);
    tb.turn_around_Z(mp_cyl2, -90);

    Assimp::Exporter exp;
    auto flags = aiProcess_ValidateDataStructure | 0;
    {
        std::string filename = (ws / "model.glb").string();
        auto err = exp.Export(model.get(), "glb2", filename, flags);
        EXPECT_EQ(AI_SUCCESS, err) << "Failed export to " << filename;
        ASSERT_TRUE(fs::is_regular_file(filename)) << filename;
    }
    {
        std::string filename = (ws / "model.xml").string();
        auto err = exp.Export(model.get(), "assxml", filename, flags);
        EXPECT_EQ(AI_SUCCESS, err) << "Failed export to " << filename;
        ASSERT_TRUE(fs::is_regular_file(filename)) << filename;
    }
}

void PrintTo(aiVector3D const &v, std::ostream *os) { *os << v; }

class TransF : public AssimpF {
public:
    aiMatrix4x4 translate(ai_real x, ai_real y, ai_real z) const {
        aiMatrix4x4 t;
        aiVector3D v(x, y, z);
        aiMatrix4Translation(&t, &v);
        return t;
    }

    aiMatrix4x4 scaling(ai_real x, ai_real y, ai_real z) const {
        aiMatrix4x4 t;
        aiVector3D v(x, y, z);
        aiMatrix4Scaling(&t, &v);
        return t;
    }

    aiMatrix4x4 scaling(ai_real m) const { return scaling(m, m, m); }

    // generated by co-pilot
    aiMatrix4x4 y_up_too() const {
        aiMatrix4x4 t;
        aiMatrix4x4::RotationX(-AI_MATH_PI / 2,
                               t); // Rotate -90 degrees around X-axis
        return t;
    }

    /// @brief Convert Z-up coordinates to Y-up (my version)
    ///
    /// | 1  0  0  0 |
    /// | 0  0  1  0 |
    /// | 0 -1  0  0 |
    /// | 0  0  0  1 |
    ///
    aiMatrix4x4 y_up() const {
        // clang-format off
        aiMatrix4x4 t = {
            1, 0, 0, 0,
            0, 0, 1, 0,
            0, -1, 0, 0,
            0, 0, 0, 1
        };
        // clang-format on
        return t;
    }

    struct Geodetic {
        double lon_deg, lat_deg, height;
    };

    auto ecef_to_geodetic(double x, double y, double z) {
        // WGS84 constants
        constexpr double a = 6378137.0;           // semi-major axis
        constexpr double f = 1.0 / 298.257223563; // flattening
        constexpr double b = a * (1 - f);         // semi-minor axis
        constexpr double e2 = 1 - (b * b) / (a * a);
        constexpr double ep2 = (a * a - b * b) / (b * b);

        double p = sqrt(x * x + y * y);
        double theta = atan2(z * a, p * b);
        double lon = atan2(y, x);
        double lat = atan2(z + ep2 * b * pow(sin(theta), 3),
                           p - e2 * a * pow(cos(theta), 3));
        double N = a / sqrt(1 - e2 * sin(lat) * sin(lat));
        double height = p / cos(lat) - N;

        Geodetic geo;
        geo.lon_deg = lon * 180.0 / M_PI;
        geo.lat_deg = lat * 180.0 / M_PI;
        geo.height = height;
        return geo;
    };
};

/// @brief test the matrix transforms
/// @param --gtest_filter=TransF.matrix_transforms
/// @param
TEST_F(TransF, matrix_transforms) {
    {
        aiVector3D v0{1, 2, 3};
        aiMatrix4x4 transform;

        aiVector3D off(10, 20, 0);
        aiMatrix4Translation(&transform, &off);
        aiVector3D v1 = transform * v0;
        CONSOLE_EVAL(v1);
        EXPECT_EQ(aiVector3D(11, 22, 3), v1);
    }
    {
        aiVector3D v0{1, 2, 3};
        auto v1 = translate(10, 20, 0) * v0;
        CONSOLE_EVAL(v1);
        EXPECT_EQ(aiVector3D(11, 22, 3), v1);
    }
    {
        aiVector3D v0{1, 2, 3};
        aiMatrix4x4 transform;

        aiVector3D S(2, 2, 2);
        aiMatrix4Scaling(&transform, &S);
        aiVector3D v1 = transform * v0;
        CONSOLE_EVAL(v1);
        EXPECT_EQ(aiVector3D(2, 4, 6), v1);
    }
    {
        aiVector3D v0{1, 2, 3};
        aiVector3D v1 = scaling(2.0) * v0;
        CONSOLE_EVAL(v1);
        EXPECT_EQ(aiVector3D(2, 4, 6), v1);
    }
    {
        // 1. translate, 2. scale
        aiVector3D v0{1, 2, 3};
        aiVector3D v1 = scaling(2.0) * translate(10, 20, 0) * v0;
        CONSOLE_EVAL(v1);
        EXPECT_EQ(aiVector3D(22, 44, 6), v1);
    }
    {
        // 1. scale, 2. translate
        aiVector3D v0{1, 2, 3};
        aiVector3D v1 = translate(10, 20, 0) * scaling(2.0) * v0;
        CONSOLE_EVAL(v1);
        EXPECT_EQ(aiVector3D(12, 24, 6), v1);
    }
    {
        CONSOLE_EVAL(y_up());
        CONSOLE_EVAL(y_up_too());
    }
    {
        // 1. Y-Up
        aiVector3D v0{1, 2, 3};
        aiVector3D v1 = y_up() * v0;
        CONSOLE_EVAL(v1);

        EXPECT_EQ(aiVector3D(1, 3, -2), v1);
    }
    {
        // 1. Y-Up (same, by with co-pilot generated function)
        aiVector3D v0{1, 2, 3};
        aiVector3D v1 = y_up_too() * v0;
        CONSOLE_EVAL(v1);
        // EXPECT_EQ(aiVector3D(1, 3, -2), v1);
        EXPECT_FLOAT_EQ(1, v1.x);
        EXPECT_FLOAT_EQ(3, v1.y);
        EXPECT_FLOAT_EQ(-2, v1.z);
    }
}

/// @brief
/// @param --gtest_filter=TransF.c3dprototype_t0
/// @param
TEST_F(TransF, c3dprototype_t0) {
    ASSERT_TRUE(fs::is_regular_file(test_data("tileset2/tileset2.json")));
    ASSERT_TRUE(fs::is_regular_file(test_data("tileset2/bounding_boxes.glb")));

    auto ws = create_ws();

    auto bounding_boxes_glb = test_data("tileset2/bounding_boxes.glb").string();
    Assimp::Importer importer;
    // clang-format off
    unsigned int postprocess_flags = 0
            | aiProcess_GenBoundingBoxes
            ;
    // clang-format on
    auto actual = importer.ReadFile(bounding_boxes_glb, postprocess_flags);
    ASSERT_TRUE(actual);

    CONSOLE_EVAL(actual->HasMeshes());
    CONSOLE_EVAL(actual->mNumMeshes);

    for (unsigned int i = 0; i < actual->mNumMeshes; ++i) {
        const aiMesh *mesh = actual->mMeshes[i];

        CONSOLE("Mesh " << i << " AABB min: " << mesh->mAABB.mMin
                        << " max: " << mesh->mAABB.mMax);
    }

    // dump mesh[0]
    {
        /*
            EPSG:32613 -- WGS 84 / UTM zone 13N
        */

        aiMesh const *mesh = actual->mMeshes[0];

        for (unsigned int v = 0; v < mesh->mNumVertices; ++v) {
            const aiVector3D &vertex = mesh->mVertices[v];
            std::cout << "Vertex " << v << ": (" << vertex.x << ", " << vertex.y
                      << ", " << vertex.z << ")\n";
        }
    }

    auto bounding_box = compute_aabb(actual);
    CONSOLE("Scene AABB min: " << bounding_box.mMin
                               << " max: " << bounding_box.mMax);

    Assimp::Exporter exporter;

    auto status =
        exporter.Export(actual, "assxml", (ws / "actual.xml").string());
    EXPECT_EQ(0, status);
}

/// @brief
/// @param --gtest_filter=TransF.c3dprototype_bb
/// @param
TEST_F(TransF, c3dprototype_bb) {
    auto actual_bb = aiAABB{{496805.969, 4420721.5, 1625.33777},
                            {496866.031, 4420784, 1636.51453}};

    aiVector3D longitude_latitude_min = {-105.0373839395, 39.9366031861, 0.0};

    // echo "496805.969 4420721.5" | cs2cs +proj=utm +zone=13 +datum=WGS84
    // +units=m +to +proj=latlong +datum=WGS84 -f "%.10f" echo "496866.031
    // 4420784" | cs2cs +proj=utm +zone=13 +datum=WGS84 +units=m +to
    // +proj=latlong +datum=WGS84 -f "%.10f"

    /*
        Cesium 3D Tiles boundingVolume.box format:
        [centerX, centerY, centerZ, halfAxisXx, halfAxisXy, halfAxisXz,
       halfAxisYx, halfAxisYy, halfAxisYz, halfAxisZx, halfAxisZy, halfAxisZz]
        where:
        - center: midpoint of AABB
        - half axes: half the size along each axis (for axis-aligned box, these
       are half extents along X, Y, Z)
    */

    aiVector3D center = (actual_bb.mMin + actual_bb.mMax) * 0.5f;
    aiVector3D half_extent = (actual_bb.mMax - actual_bb.mMin) * 0.5f;

    // Axis-aligned box: half axes are just half extents along X, Y, Z
    std::array<double, 12> cesium_box = {
        center.x, center.y, center.z, half_extent.x, 0, 0, 0, half_extent.y,
        0,        0,        0,        half_extent.z};

    CONSOLE("Cesium boundingVolume.box: ["
            << std::setprecision(16) << cesium_box[0] << ", " << cesium_box[1]
            << ", " << cesium_box[2] << ", " << cesium_box[3] << ", "
            << cesium_box[4] << ", " << cesium_box[5] << ", " << cesium_box[6]
            << ", " << cesium_box[7] << ", " << cesium_box[8] << ", "
            << cesium_box[9] << ", " << cesium_box[10] << ", " << cesium_box[11]
            << "]");

    /*
        To transform from local coordinates (e.g., UTM) to ECEF (Earth-Centered,
       Earth-Fixed), you need to:
        1. Convert local coordinates to geodetic (longitude, latitude, height).
        2. Convert geodetic to ECEF.

        This requires a geodetic library (e.g., proj, GeographicLib).
        Here is a simple ECEF conversion for WGS84 ellipsoid.
    */

    struct Geodetic {
        double lon_deg;
        double lat_deg;
        double height;
    };

    struct ECEF {
        double x;
        double y;
        double z;
    };

    auto geodetic_to_ecef = [](const Geodetic &geo) -> ECEF {
        // WGS84 constants
        constexpr double a = 6378137.0;           // semi-major axis
        constexpr double f = 1.0 / 298.257223563; // flattening
        constexpr double b = a * (1 - f);         // semi-minor axis
        constexpr double e2 = 1 - (b * b) / (a * a);

        double lon = geo.lon_deg * M_PI / 180.0;
        double lat = geo.lat_deg * M_PI / 180.0;
        double N = a / sqrt(1 - e2 * sin(lat) * sin(lat));
        double x = (N + geo.height) * cos(lat) * cos(lon);
        double y = (N + geo.height) * cos(lat) * sin(lon);
        double z = ((b * b) / (a * a) * N + geo.height) * sin(lat);
        return {x, y, z};
    };

    // Example: convert bounding box center from UTM to ECEF
    // You need to convert UTM to lon/lat first (use proj or similar).
    // Here, we use the provided longitude_latitude_min as an example.

    Geodetic geo_center;
    geo_center.lon_deg = longitude_latitude_min.x; // longitude
    geo_center.lat_deg = longitude_latitude_min.y; // latitude
    geo_center.height = center.z;                  // use Z as height

    auto ecef_center = geodetic_to_ecef(geo_center);

    CONSOLE("ECEF center: x=" << ecef_center.x << " y=" << ecef_center.y
                              << " z=" << ecef_center.z);
}

/// @brief
/// @param --gtest_filter=TransF.c3dprototype_triangulate
/// @param
TEST_F(TransF, c3dprototype_triangulate) {
    ASSERT_TRUE(fs::is_regular_file(test_data("tileset2/tileset2.json")));
    ASSERT_TRUE(fs::is_regular_file(test_data("tileset2/bounding_boxes.glb")));

    auto ws = create_ws();

    auto bounding_boxes_glb = test_data("tileset2/bounding_boxes.glb").string();
    Assimp::Importer importer;
    // clang-format off
        unsigned int postprocess_flags = 0
            // | aiProcess_GenBoundingBoxes
            // | aiProcess_ValidateDataStructure
            // | aiProcess_CalcTangentSpace
            | aiProcess_Triangulate
            // | aiProcess_GenNormals
            // | aiProcess_JoinIdenticalVertices
            ;
    // clang-format on
    auto actual = importer.ReadFile(bounding_boxes_glb, postprocess_flags);
    ASSERT_TRUE(actual);

    CONSOLE_EVAL(actual->HasMeshes());
    CONSOLE_EVAL(actual->mNumMeshes);

    Assimp::Exporter exporter;

    auto status = exporter.Export(actual, "glb", (ws / "actual.glb").string());
    EXPECT_EQ(0, status);
}

/// @brief Load the original GLB model and translate the coordinates to (0,0,0)
/// @param --gtest_filter=TransF.c3dprototype_translate_coordinates
/// @param
TEST_F(TransF, c3dprototype_translate_coordinates) {
    ASSERT_TRUE(fs::is_regular_file(test_data("tileset2/tileset2.json")));
    ASSERT_TRUE(fs::is_regular_file(test_data("tileset2/bounding_boxes.glb")));

    auto ws = create_ws();
    {
        auto bounding_boxes_glb =
            test_data("tileset2/bounding_boxes.glb").string();
        Assimp::Importer importer;
        // clang-format off
        unsigned int postprocess_flags = 0
            //| aiProcess_GenBoundingBoxes
            //| aiProcess_ValidateDataStructure
            //| aiProcess_CalcTangentSpace
            //| aiProcess_Triangulate
            //| aiProcess_GenNormals
            //| aiProcess_JoinIdenticalVertices
            ;
        // clang-format on

        // The original coordinates are in UTM 13N
        auto actual = importer.ReadFile(bounding_boxes_glb, postprocess_flags);
        ASSERT_TRUE(actual);

        CONSOLE(std::setprecision(15));

        CONSOLE_EVAL(actual->mMeshes[0]->mAABB);
        CONSOLE_EVAL(actual->mRootNode->mTransformation);

        EXPECT_EQ(aiMatrix4x4(1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 1),
                  actual->mRootNode->mTransformation);

        auto bounding_box = compute_aabb(actual);
        CONSOLE_EVAL(bounding_box.mMin);
        CONSOLE_EVAL(bounding_box.mMax);
        // mMin : {496805.969, 4420721.5, 1625.33777}
        // mMax : {496866.031, 4420784, 1636.51453}

        aiMatrix4x4 translate;
        aiMatrix4x4::Translation(-bounding_box.mMin, translate);
        CONSOLE_EVAL(translate);
        {
            // Verify sure the transform works as expected - moves min_point to
            // (0,0,0)
            aiVector3D min_point{496805.969, 4420721.5, 1625.33777};
            EXPECT_EQ(min_point, bounding_box.mMin);
            CONSOLE_EVAL(translate * min_point);
            EXPECT_TRUE(aiVector3D(0, 0, 0) == translate * min_point);
            aiVector3D max_point{496866.031, 4420784, 1636.51453};
            CONSOLE_EVAL(translate * max_point);
            EXPECT_TRUE(aiVector3D(60.0625, 62.5, 11.1767578) ==
                        translate * max_point);
            EXPECT_TRUE(max_point - min_point == translate * max_point);
        }

        /* update the root transform to map min to (0,0,0) */
        actual->mRootNode->mTransformation =
            translate * actual->mRootNode->mTransformation;
        CONSOLE_EVAL(actual->mRootNode->mTransformation);

        Assimp::Exporter exporter;
        auto actual_status =
            exporter.Export(actual, "glb2", (ws / "actual.glb").string());
        ASSERT_EQ(0, actual_status);
        ASSERT_EQ(
            0, exporter.Export(actual, "assxml", (ws / "actual.xml").string()));
        ASSERT_EQ(0,
                  exporter.Export(actual, "ply", (ws / "actual.ply").string()));
    }
    // verify the saved model
    {
#if 0
        // Create default logger that outputs to console
        Assimp::DefaultLogger::create("", Assimp::Logger::VERBOSE,
                                      aiDefaultLogStream_STDOUT);

#endif
        Assimp::Importer importer;
        auto model = importer.ReadFile((ws / "actual.glb").string(), 0);
        ASSERT_TRUE(model);
        CONSOLE_EVAL(model->mRootNode->mTransformation);
        CONSOLE_EVAL(model->mNumMeshes);
        EXPECT_EQ(9, model->mNumMeshes);
    }
    {
        // verify with tinygltf::TinyGLTF
        {
            tinygltf::TinyGLTF loader;
            tinygltf::Model model;
            std::string err, warn;

            bool success = loader.LoadBinaryFromFile(
                &model, &err, &warn, (ws / "actual.glb").string());

            ASSERT_TRUE(success) << "Failed to load GLTF: " << err;
            EXPECT_TRUE(warn.empty()) << "Warning: " << warn;
            CONSOLE_EVAL(model.meshes.size());
            EXPECT_EQ(9, model.meshes.size());
            EXPECT_EQ(10, model.nodes.size());

            auto vector2str = [](std::vector<double> const &v) {
                std::ostringstream ss;
                ss << "{";
                auto vv = v.begin();
                if (vv != v.end())
                    ss << *vv++;
                while (vv != v.end())
                    ss << ", " << *vv++;
                ss << "}";
                return ss.str();
            };

            for (auto const &np : model.nodes) {
                CONSOLE_EVAL(np.extensions_json_string);
                CONSOLE_EVAL(np.children.size());
                CONSOLE_EVAL(vector2str(np.matrix));
                CONSOLE("\n");
            }
        }
    }
}

/*
    echo "EASTING NORTHING" | cs2cs +proj=utm +zone=13 +ellps=WGS84 +to
   +proj=latlong +datum=WGS84

    in decimal format

    echo "EASTING NORTHING" | cs2cs +proj=utm +zone=13 +ellps=WGS84 +to
   +proj=latlong +ellps=WGS84 -f "%.8f"

    To produce a tileset json

    npx 3d-tiles-tools createTilesetJson -f -i cropped/model.glb -o
   cropped/model.json --cartographicPositionDegrees -105.03697644 39.93681616
   0.00000000
*/

static const char *cropped_model_glb = R"({
  "asset": {
    "version": "1.1"
  },
  "geometricError": 4096,
  "root": {
    "boundingVolume": {
      "box": [
        496840.7947318554,
        4420745.123109758,
        1628.3355102539062,
        1.4763270616531372,
        0,
        0,
        0,
        1.3503844738006592,
        0,
        0,
        0,
        1.38446044921875
      ]
    },
    "geometricError": 512,
    "content": {
      "uri": "model.glb"
    },
    "refine": "ADD",
    "transform": [
      0.9657585935309204,
      -0.2594423616551052,
      0,
      0,
      0.16654706558975313,
      0.6199614388126924,
      0.7667528215329492,
      0,
      -0.1989281628242237,
      -0.7404981265095258,
      0.6419424512144849,
      0,
      -1270544.8042035468,
      -4729526.648222642,
      4072608.8590562753,
      1
    ]
  }
})";

/// @brief Load a GLB file and generate tileset.json
/// @param --gtest_filter=TransF.c3dprototype_make_tileset_json
/// @param
TEST_F(TransF, c3dprototype_make_tileset_json) {
    auto ws = create_ws();
    {
        // Assimp::DefaultLogger::create("", Assimp::Logger::VERBOSE,
        // aiDefaultLogStream_STDOUT);
        auto model_glb =
            test_data("CesiumTilesF/GM20993_export_cropped_as_cesium3d_tile/"
                      "model.glb")
                .string();
        ASSERT_TRUE(fs::is_regular_file(model_glb));

        Assimp::Importer importer;
        // clang-format off
        unsigned int postprocess_flags = 0
            //| aiProcess_GenBoundingBoxes
            //| aiProcess_ValidateDataStructure
            //| aiProcess_CalcTangentSpace
            //| aiProcess_Triangulate
            //| aiProcess_GenNormals
            //| aiProcess_JoinIdenticalVertices
            ;
        // clang-format on

        // The original coordinates are in UTM 13N
        auto actual = importer.ReadFile(model_glb, postprocess_flags);
        ASSERT_TRUE(actual) << importer.GetErrorString();

        CONSOLE_EVAL(compute_aabb(actual));
        CONSOLE_EVAL(compute_aabb_with_transform(actual));
    }
    {
        double b[] = {496840.7947318554,
                      4420745.123109758,
                      1628.3355102539062,
                      1.4763270616531372,
                      0,
                      0,
                      0,
                      1.3503844738006592,
                      0,
                      0,
                      0,
                      1.38446044921875};

        double t[] = {0.9657585935309204,
                      -0.2594423616551052,
                      0,
                      0,
                      0.16654706558975313,
                      0.6199614388126924,
                      0.7667528215329492,
                      0,
                      -0.1989281628242237,
                      -0.7404981265095258,
                      0.6419424512144849,
                      0,
                      -1270544.8042035468,
                      -4729526.648222642,
                      4072608.8590562753,
                      1};
        aiVector3t<double> center{t[0], t[1], t[2]};
        aiMatrix4x4t<double> transform{t[0], t[4],  t[8],  t[12], t[1],  t[5],
                                       t[9], t[13], t[2],  t[6],  t[10], t[14],
                                       t[3], t[7],  t[11], t[15]};

        auto center_transformed = transform * center;
        CONSOLE_EVAL(center_transformed);

        // original: -105.03697644 39.93681616 0.00000000
        // Geodetic: lon=-105.036965 lat=39.9368138 height=7.72997737e-08

        auto geo = ecef_to_geodetic(center_transformed.x, center_transformed.y,
                                    center_transformed.z);
        CONSOLE("Geodetic: lon=" << geo.lon_deg << " lat=" << geo.lat_deg
                                 << " height=" << geo.height);
    }
}

/// @brief
/// @param --gtest_filter=TransF.c3dprototype_box_utm_13N_coordinates
/// @param
TEST_F(TransF, c3dprototype_box_utm_13N_coordinates) {

    aiVector3d p0{496840, 0, -4420750}; // Y-UP coordinates

    {

        auto ws = create_ws();

        meshtoolbox::Toolbox tb;

        meshtoolbox::box_t b0{p0, {10, 20, 30}};

        std::vector<meshtoolbox::box_t> boxes{b0};

        auto model = std::unique_ptr<aiScene>(tb.make_boxes(boxes));

        ASSERT_TRUE(model);
        EXPECT_EQ(1, model->mNumMeshes);

        Assimp::Exporter exp;
        auto flags =
            /*aiProcess_GenNormals | */ aiProcess_ValidateDataStructure | 0;
        {
            std::string filename_glb = (ws / "model.glft").string();
            auto err = exp.Export(model.get(), "gltf", filename_glb, flags);
            EXPECT_EQ(AI_SUCCESS, err) << "Failed export to " << filename_glb;
        }
        {
            std::string filename_glb = (ws / "model.glb").string();
            auto err = exp.Export(model.get(), "glb2", filename_glb, flags);
            EXPECT_EQ(AI_SUCCESS, err) << "Failed export to " << filename_glb;
        }
        {
            std::string filename_glb = (ws / "model.xml").string();
            auto err = exp.Export(model.get(), "assxml", filename_glb, flags);
            EXPECT_EQ(AI_SUCCESS, err) << "Failed export to " << filename_glb;
        }
    }
    {
        double t[] = {0.965758551283516,
                      -0.259442518918477,
                      0,
                      0,
                      0.16654731910137624,
                      0.6199619795783732,
                      0.7667523292285515,
                      0,
                      -0.19892815568166483,
                      -0.7404976186690275,
                      0.641943039235251,
                      0,
                      -1270544.7618038643,
                      -4729523.416653709,
                      4072612.5998952053,
                      1};

        static_assert(sizeof(t) / sizeof(*t) == 16);

        aiMatrix4x4t<double> transform{t[0], t[4],  t[8],  t[12], t[1],  t[5],
                                       t[9], t[13], t[2],  t[6],  t[10], t[14],
                                       t[3], t[7],  t[11], t[15]};

        aiVector3d p0_zup{496840, 4420750, 0}; // Y-UP coordinates

        auto p0_transformed = transform * p0_zup;
        CONSOLE_EVAL(p0_transformed);

        auto geo = ecef_to_geodetic(p0_transformed.x, p0_transformed.y,
                                    p0_transformed.z);
        CONSOLE("Geodetic: lon=" << geo.lon_deg << " lat=" << geo.lat_deg
                                 << " height=" << geo.height);
    }
}
#endif // if 0

/// @brief Read a small jpg file with stb library
/// @param --gtest_filter=AssimpF.stb_read_jpg_too
/// @param
TEST_F(AssimpF, stb_read_jpg_too) {
    auto logo_jpg = test_data("ctzb3.jpg");
    ASSERT_TRUE(fs::is_regular_file(logo_jpg));

    int width, height, channels;
    auto stbi_uc_deleter = [](stbi_uc *p) {
        if (p)
            stbi_image_free(p);
    };
    auto data = std::unique_ptr<stbi_uc, decltype(stbi_uc_deleter)>(
        stbi_load(logo_jpg.string().c_str(), &width, &height, &channels, 0),
        stbi_uc_deleter);

    ASSERT_TRUE(data) << "Failed: " << logo_jpg
                      << ", reason:" << stbi_failure_reason();

    ASSERT_EQ(600, width);
    ASSERT_EQ(418, height);
    ASSERT_EQ(3, channels);

    // Optionally, check a pixel value
    if (width > 0 && height > 0 && channels >= 3) {
        int idx = 0; // top-left pixel
        unsigned char r = data.get()[idx * channels + 0];
        unsigned char g = data.get()[idx * channels + 1];
        unsigned char b = data.get()[idx * channels + 2];
        CONSOLE_EVAL(unsigned(r));
        CONSOLE_EVAL(unsigned(g));
        CONSOLE_EVAL(unsigned(b));
    }
}
