//
//
//

#include <gtest/gtest.h>

#include <Eigen/Sparse>
#include <proj.h>

#include "../common/DebuggingConsole.h"

class EigenF : public testing::Test {};

TEST_F(EigenF, t0) {
    CONSOLE("Hello...");

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