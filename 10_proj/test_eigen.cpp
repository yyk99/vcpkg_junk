//
//
//

#include <gtest/gtest.h>

#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <proj.h>

#include "../common/DebuggingConsole.h"

class EigenF : public testing::Test {};

TEST_F(EigenF, version) {
    CONSOLE(EIGEN_WORLD_VERSION);
    CONSOLE(EIGEN_MAJOR_VERSION);
    CONSOLE(EIGEN_MINOR_VERSION);
    EXPECT_EQ(3, EIGEN_WORLD_VERSION);
    EXPECT_EQ(4, EIGEN_MAJOR_VERSION);
}

TEST_F(EigenF, t0) {
    // clang-format off
    Eigen::Matrix4d yup2zup;
    yup2zup << 
        1.0, 0.0, 0.0, 0.0,
        0.0, 0.0, -1.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 1.0;
    // clang-format on

    Eigen::Vector4d pos = {1, 2, 3, 1};
    CONSOLE("pos:\n" << pos);

    pos = yup2zup * pos;
    CONSOLE("pos:\n" << pos);
}

TEST_F(EigenF, translate_transform) {
    Eigen::Vector3f vv{1.0f, 2.0f, 3.0f};
    Eigen::Translation<float, 3> T{1.0f, 2.0f, 3.0f};

    CONSOLE_EVAL(vv.transpose());
    vv = T * vv;
    CONSOLE_EVAL(vv.transpose());

    CONSOLE_EVAL(Eigen::Vector4f::Identity());

    CONSOLE_EVAL(Eigen::Matrix4f::Identity());
}