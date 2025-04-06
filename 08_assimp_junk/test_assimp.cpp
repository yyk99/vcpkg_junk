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
};

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

    auto cube_ply = (fs::path(__FILE__).parent_path() / "test_data/cube.ply").string();

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
