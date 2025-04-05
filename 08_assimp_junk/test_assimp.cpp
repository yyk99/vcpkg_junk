//
//
//

#include <gtest/gtest.h>

#include <assimp/Importer.hpp>
#include <assimp/Exporter.hpp>
#include <assimp/scene.h>

#include <filesystem>
#include <string>

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
};

TEST_F(AssimpF, toolkit_create_ws) {
    auto ws = create_ws();
    ASSERT_TRUE(fs::is_directory(ws));
    std::cout << ws << "\n";
}

TEST_F(AssimpF, z1) {
    Assimp::Exporter exp;

    std::cout << fs::absolute(".") << "\n";
    std::cout << test_name() << "\n";
    std::cout << test_case_name() << "\n";
}
