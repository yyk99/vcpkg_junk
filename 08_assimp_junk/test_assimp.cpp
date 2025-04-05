//
//
//

#include <gtest/gtest.h>

#include <assimp/Importer.hpp>
#include <assimp/Exporter.hpp>
#include <assimp/scene.h>

#include <filesystem>
#include <string>
#include <fstream>

namespace fs = std::filesystem;

class AssimpF : public testing::Test {
protected:
    std::string test_name() const {
        return ::testing::UnitTest::GetInstance()->current_test_info()->name();
    }

    std::string test_case_name() const {
        return ::testing::UnitTest::GetInstance()->current_test_info()->test_case_name();
    }

    fs::path create_ws() {
        auto ws = fs::absolute("out") / test_case_name() / test_name();
        fs::remove_all(ws);
        fs::create_directories(ws);
        return ws;
    }

    void save_as(std::string const &text, fs::path const &location)
    {
        std::ofstream os(location.string().c_str());
        if(!os.good())
            throw std::runtime_error("Cannot create file: " + location.string());
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

TEST_F(AssimpF, z1) {
    Assimp::Exporter exp;

    std::cout << fs::absolute(".") << "\n";
    std::cout << test_name() << "\n";
    std::cout << test_case_name() << "\n";
}
