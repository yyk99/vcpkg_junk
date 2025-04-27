//
//
//

#include <gtest/gtest.h>

#include <math.h>
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

    /* or you may set C=PJ_DEFAULT_CTX if you are sure you will     */
    /* use PJ objects from only one thread                          */
    C = proj_context_create();

    P = proj_create_crs_to_crs(
        C, "EPSG:4326", "+proj=utm +zone=32 +datum=WGS84", /* or EPSG:32632 */
        NULL);

    ASSERT_TRUE(P) << "Failed to create transformation object";

    /* This will ensure that the order of coordinates for the input CRS */
    /* will be longitude, latitude, whereas EPSG:4326 mandates latitude, */
    /* longitude */
    norm = proj_normalize_for_visualization(C, P);
    ASSERT_TRUE(norm) << "Failed to normalize transformation object";

    proj_destroy(P);
    P = norm;

    /* a coordinate union representing Copenhagen: 55d N, 12d E */
    /* Given that we have used proj_normalize_for_visualization(), the order */
    /* of coordinates is longitude, latitude, and values are expressed in */
    /* degrees. */
    a = proj_coord(12, 55, 0, 0);

    /* transform to UTM zone 32, then back to geographical */
    b = proj_trans(P, PJ_FWD, a);
    printf("easting: %.3f, northing: %.3f\n", b.enu.e, b.enu.n);

    b = proj_trans(P, PJ_INV, b);
    printf("longitude: %g, latitude: %g\n", b.lp.lam, b.lp.phi);
    EXPECT_DOUBLE_EQ(12, b.lp.lam);
    EXPECT_DOUBLE_EQ(55, b.lp.phi);

    /* Clean up */
    proj_destroy(P);
    proj_context_destroy(C); /* may be omitted in the single threaded case */
}

/// @brief
/// @param
/// @param
TEST_F(ProjF, c_api_primer_2) {

    /* Create the context. */
    /* You may set C=PJ_DEFAULT_CTX if you are sure you will     */
    /* use PJ objects from only one thread                       */
    PJ_CONTEXT *C = proj_context_create();

    /* Create a projection. */
    PJ *P = proj_create(C, "+proj=utm +zone=32 +datum=WGS84 +type=crs");

    ASSERT_TRUE(P) << "Failed to create transformation object";

    /* Get the geodetic CRS for that projection. */
    PJ *G = proj_crs_get_geodetic_crs(C, P);

    /* Create the transform from geodetic to projected coordinates.*/
    PJ_AREA *A = NULL;
    const char *const *options = NULL;
    PJ *G2P = proj_create_crs_to_crs_from_pj(C, G, P, A, options);

    /* Longitude and latitude of Copenhagen, in degrees. */
    double lon = 12.0, lat = 55.0;

    /* Prepare the input */
    PJ_COORD c_in;
    c_in.lpzt.z = 0.0;
    c_in.lpzt.t = HUGE_VAL; // important only for time-dependent projections
    c_in.lp.lam = lon;
    c_in.lp.phi = lat;
    printf("Input longitude: %g, latitude: %g (degrees)\n", c_in.lp.lam,
           c_in.lp.phi);

    /* Compute easting and northing */
    PJ_COORD c_out = proj_trans(G2P, PJ_FWD, c_in);
    printf("Output easting: %g, northing: %g (meters)\n", c_out.enu.e,
           c_out.enu.n);

    /* Apply the inverse transform */
    PJ_COORD c_inv = proj_trans(G2P, PJ_INV, c_out);
    printf("Inverse applied. Longitude: %g, latitude: %g (degrees)\n",
           c_inv.lp.lam, c_inv.lp.phi);
    EXPECT_DOUBLE_EQ(12, c_inv.lp.lam);
    EXPECT_DOUBLE_EQ(55, c_inv.lp.phi);

    /* Clean up */
    proj_destroy(P);
    proj_destroy(G);
    proj_destroy(G2P);
    proj_context_destroy(C); /* may be omitted in the single threaded case */
}