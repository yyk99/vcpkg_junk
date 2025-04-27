//
//
//

#include <gtest/gtest.h>

#include <proj.h>
#include <stdio.h>

#include "../common/DebuggingConsole.h"

class ProjF : public testing::Test {};

TEST_F(ProjF, version) {
    // #define PROJ_VERSION_MAJOR 9
    // #define PROJ_VERSION_MINOR 4
    // #define PROJ_VERSION_PATCH 0

    CONSOLE("VERSION: " << PROJ_VERSION_MAJOR << "." << PROJ_VERSION_MINOR
                        << "." << PROJ_VERSION_PATCH);
    EXPECT_EQ(9, PROJ_VERSION_MAJOR);
}

/// @brief A c-api primer
/// @param --gtest_filter=ProjF.c_api_primer
/// @param
TEST_F(ProjF, c_api_primer) {
    PJ_CONTEXT *C;
    PJ *P;
    PJ *norm;
    PJ_COORD a, b;

    C = proj_context_create();
    ASSERT_TRUE(C);

    P = proj_create_crs_to_crs(
        C, "EPSG:4326", "+proj=utm +zone=32 +datum=WGS84", /* or EPSG:32632 */
        NULL);

    ASSERT_TRUE(P) << "Failed to create transformation object";
}