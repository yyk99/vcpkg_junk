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

/// @brief Set the location of proj.db
/// @param --gtest_filter=ProjF.c_api_primer_find_proj_db
/// @param
TEST_F(ProjF, c_api_primer_find_proj_db) {
    PJ_CONTEXT *C;
    PJ *P;
    PJ *norm;
    PJ_COORD a, b;

    /* or you may set C=PJ_DEFAULT_CTX if you are sure you will     */
    /* use PJ objects from only one thread                          */
    C = proj_context_create();

    /* */
#ifdef PROJ_DIR
    char const *data_path[] = {PROJ_DIR};
    proj_context_set_search_paths(C, 1, data_path);
#endif
    P = proj_create_crs_to_crs(
        C, "EPSG:4326", "+proj=utm +zone=32 +datum=WGS84", /* or EPSG:32632 */
        NULL);

    ASSERT_TRUE(P) << "Failed to create transformation object";

    /* Clean up */
    proj_destroy(P);
    proj_context_destroy(C); /* may be omitted in the single threaded case */
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
    /* */
#ifdef PROJ_DIR
    char const *data_path[] = {PROJ_DIR};
    proj_context_set_search_paths(C, 1, data_path);
#endif

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
    /* */
#ifdef PROJ_DIR
    char const *data_path[] = {PROJ_DIR};
    proj_context_set_search_paths(C, 1, data_path);
#endif

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

#include <cassert>
#include <cmath>   // for HUGE_VAL
#include <iomanip> // for std::setprecision()

#include <iostream>

#include "proj/coordinateoperation.hpp"
#include "proj/crs.hpp"
#include "proj/io.hpp"
#include "proj/util.hpp" // for nn_dynamic_pointer_cast

using namespace NS_PROJ::crs;
using namespace NS_PROJ::io;
using namespace NS_PROJ::operation;
using namespace NS_PROJ::util;

/// @see https://proj.org/en/stable/development/quickstart_cpp.html

///
TEST_F(ProjF, cpp_api) {
    auto dbContext = DatabaseContext::create(); 

    // Instantiate a generic authority factory, that is not tied to a particular
    // authority, to be able to get transformations registered by different
    // authorities. This can only be used for CoordinateOperationContext.
    auto authFactory = AuthorityFactory::create(dbContext, std::string());

    // Create a coordinate operation context, that can be customized to amend
    // the way coordinate operations are computed. Here we ask for default
    // settings.
    auto coord_op_ctxt =
        CoordinateOperationContext::create(authFactory, nullptr, 0.0);

    // Instantiate a authority factory for EPSG related objects.
    auto authFactoryEPSG = AuthorityFactory::create(dbContext, "EPSG");

    // Instantiate source CRS from EPSG code
    auto sourceCRS = authFactoryEPSG->createCoordinateReferenceSystem("4326");

    // Instantiate target CRS from PROJ.4 string (commented out, the equivalent
    // from the EPSG code)
#if 1
    auto targetCRS = authFactoryEPSG->createCoordinateReferenceSystem("32631");
#else
    auto targetCRS =
        NN_CHECK_THROW(nn_dynamic_pointer_cast<CRS>(createFromUserInput(
            "+proj=utm +zone=31 +datum=WGS84 +type=crs", dbContext)));
#endif
    // List operations available to transform from EPSG:4326
    // (WGS 84 latitude/longitude) to EPSG:32631 (WGS 84 / UTM zone 31N).
    auto list = CoordinateOperationFactory::create()->createOperations(
        sourceCRS, targetCRS, coord_op_ctxt);
    CONSOLE_EVAL(list.size());

    // Check that we got a non-empty list of operations
    // The list is sorted from the most relevant to the less relevant one.
    // Cf
    // https://proj.org/operations/operations_computation.html#filtering-and-sorting-of-coordinate-operations
    // for more details on the sorting of those operations.
    // For a transformation between a projected CRS and its base CRS, like
    // we do here, there will be only one operation.
    ASSERT_TRUE(!list.empty());

    // Create an execution context (must only be used by one thread at a time)
    PJ_CONTEXT *ctx = proj_context_create();

    // Create a coordinate transformer from the first operation of the list
    auto transformer = list[0]->coordinateTransformer(ctx);

    // Perform the coordinate transformation.
    PJ_COORD c = {{
        49.0,    // latitude in degree
        2.0,     // longitude in degree
        0.0,     // z ordinate. unused
        HUGE_VAL // time ordinate. unused
    }};
    c = transformer->transform(c);

    // Display result
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Easting: " << c.v[0] << std::endl;  // should be 426857.988
    std::cout << "Northing: " << c.v[1] << std::endl; // should be 5427937.523

    ASSERT_NEAR(426857.988, c.v[0], 0.001);
    ASSERT_NEAR(5427937.523, c.v[1], 0.001);

    ASSERT_FLOAT_EQ(426857.988, c.v[0]);
    ASSERT_FLOAT_EQ(5427937.523, c.v[1]);

    // Destroy execution context
    proj_context_destroy(ctx);
}

#include <proj.h>
#include <stdio.h>

TEST_F(ProjF, c_error_handling) {
    PJ_CONTEXT *c;
    PJ *p;
    int err;
    const char *errstr;

    c = proj_context_create();
    p = proj_create_crs_to_crs(c, "EPSG:4326", "EPSG:XXXX", NULL);
    /* it is expected to fail */
    ASSERT_TRUE(p == 0);
    /* Something is wrong, let's try to get details ... */
    err = proj_context_errno(c);
    ASSERT_NE(0, err) << "Failed to create transformation, reason unknown";

    errstr = proj_context_errno_string(c, err);
    CONSOLE("Failed to create transformation: " << errstr);

    proj_context_destroy(c);
}

///
///
///

// EPSG:4326 is a common code in GIS that stands for the WGS84 geographic
// coordinate system. It's the coordinate system used by GPS and Google Earth,
// representing the Earth as a three-dimensional ellipsoid with latitude and
// longitude coordinates. It's often used when accuracy over large areas is
// important, and it doesn't have the distortions of projected systems like Web
// Mercator (EPSG:3857).

// EPSG:32631 (WGS 84 / UTM zone 31N).

// EPSG:4978 - WGS 84 - earth centered CS

// Portland, Maine is located in UTM zone 19 North.
// This zone is used to describe the location of North Deering,
// Portland, ME using the UTM coordinate system.

// 19N	WGS84 / UTM zone 19N	EPSG:32619

TEST_F(ProjF, frost_hill_rd_68) {
    PJ_CONTEXT *ctx = proj_context_create();
    PJ *p = proj_create_crs_to_crs(ctx, "EPSG:4326", "EPSG:32619", NULL);
    ASSERT_TRUE(p);

    PJ_COORD coord = {{
        43.717449894149915, // latitude
        -70.28796806594934, // longitude
        0,                  // elevation
        HUGE_VAL            // hz...
    }};

    PJ *p_norm = proj_normalize_for_visualization(ctx, p);

    auto a = proj_coord(12, 55, 0, 0);

    /* transform to UTM zone 32, then back to geographical */
    auto b = proj_trans(p_norm, PJ_FWD, a);
    printf("easting: %.3f, northing: %.3f\n", b.enu.e, b.enu.n);

    b = proj_trans(p_norm, PJ_INV, b);
    printf("longitude: %g, latitude: %g\n", b.lp.lam, b.lp.phi);
    EXPECT_NEAR(12, b.lp.lam, 1E-12);
    EXPECT_DOUBLE_EQ(55, b.lp.phi);

    proj_context_destroy(ctx);
}

//
// proj4 pipelines
//
// This code defines a pipeline that performs a coordinate transformation
// from ED50/UTM32 to ETRS89/UTM33. It initializes a PJ_CONTEXT and a
// PJ object with the pipeline definition.
// Then, it defines an input coordinate a and applies the transformation using
// proj_trans to get the output coordinate b. Finally, it prints the transformed
// coordinates and cleans up the allocated resources.
//
TEST_F(ProjF, pipelines) {
    PJ_CONTEXT *C;
    PJ *P;
    PJ_COORD a, b;

    C = proj_context_create();
    ASSERT_TRUE(C);

    P = proj_create(C,
                    "+proj=pipeline "
                    "+step +inv +proj=utm +zone=32 +ellps=intl "
                    "+step +proj=cart +ellps=intl "
                    "+step +proj=helmert +x=-81.0703 +y=-89.3603 +z=-115.7526 "
                    "+rx=-0.48488 +ry=-0.02436 +rz=-0.41321 +s=-0.540645 "
                    "   +convention=coordinate_frame "
                    "+step +inv +proj=cart +ellps=GRS80 "
                    "+step +proj=utm +zone=33 +ellps=GRS80");

    ASSERT_TRUE(P) << "Error creating pipeline: "
                   << proj_errno_string(proj_context_errno(C));

    a = proj_coord(488200, 5385800, 0, 0);
    b = proj_trans(P, PJ_FWD, a);

    CONSOLE(std::setprecision(16));
    CONSOLE("Easting: " << b.enu.e << " Northing: " << b.enu.n);

    EXPECT_DOUBLE_EQ(46077.984547946835, b.enu.e);
    EXPECT_DOUBLE_EQ(5403933.1579267327, b.enu.n);

    proj_destroy(P);
    proj_context_destroy(C);
}

// EPSG:21781 - swiss coordinate system LV03
// EPSG:4978 - WGS 84 - earth centered CS
// EPSG:4326 - the WGS84 geographic coordinate system.
TEST_F(ProjF, lv03) {
    const char wkt_string[] = R"WKT(
    PROJCS["CH1903 / LV03",
        GEOGCS["CH1903",
        DATUM["CH1903",SPHEROID["Bessel 1841",6377397.155,299.1528128],
            TOWGS84[674.374,15.056,405.346,0,0,0,0]],
            PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],
            UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],
            AUTHORITY["EPSG","4149"]],PROJECTION["Hotine_Oblique_Mercator_Azimuth_Center"],
            PARAMETER["latitude_of_center",46.9524055555556],
            PARAMETER["longitude_of_center",7.43958333333333],
            PARAMETER["azimuth",90],
            PARAMETER["rectified_grid_angle",90],
            PARAMETER["scale_factor",1],
            PARAMETER["false_easting",600000],
            PARAMETER["false_northing",200000],
            UNIT["metre",1,
            AUTHORITY["EPSG","9001"]],
            AXIS["Easting",EAST],
            AXIS["Northing",NORTH],
            AUTHORITY["EPSG","21781"]])WKT";

    PJ_CONTEXT *ctx = proj_context_create();
    {
        PJ *p = proj_create_crs_to_crs(ctx, "EPSG:4326", "EPSG:21781", NULL);
        ASSERT_TRUE(p);

        {
            PJ_COORD coord = {{
                46.9524055555556, // latitude
                7.43958333333333, // longitude
                0,                // elevation
                HUGE_VAL          // hz...
            }};

            PJ *p_norm = proj_normalize_for_visualization(ctx, p);

            PJ_COORD a = proj_trans(p, PJ_FWD, coord);

            printf("easting: %.3f, northing: %.3f\n", a.enu.e, a.enu.n);
            EXPECT_DOUBLE_EQ(600072.38970262511, a.enu.e);
            EXPECT_DOUBLE_EQ(200147.05558247434, a.enu.n);
        }
        {
            // (Sidlerstrasse 5 - 46°57'3.9" N, 7°26'19.1" E).
            PJ_COORD coord = {{
                46.0 + 57 / 60. + 3.9 / 3600., // latitude
                7.0 + 26 / 60. + 19.1 / 3600., // longitude
                0,                             // elevation
                HUGE_VAL                       // hz...
            }};

            PJ *p_norm = proj_normalize_for_visualization(ctx, p);

            PJ_COORD a = proj_trans(p, PJ_FWD, coord);

            printf("easting: %.3f, northing: %.3f\n", a.enu.e, a.enu.n);
            EXPECT_DOUBLE_EQ(600000.49294705561, a.enu.e);
            EXPECT_DOUBLE_EQ(200000.0635578575, a.enu.n);
        }
    }
    {
        PJ *p = proj_create_crs_to_crs(ctx, "EPSG:4326", wkt_string, NULL);
        ASSERT_TRUE(p);

        {
            PJ_COORD coord = {{
                46.9524055555556, // latitude
                7.43958333333333, // longitude
                0,                // elevation
                HUGE_VAL          // hz...
            }};
            printf("easting: %.16f, northing: %.16f\n", coord.lpzt.lam,
                   coord.lpzt.phi);

            PJ *p_norm = proj_normalize_for_visualization(ctx, p);

            PJ_COORD a = proj_trans(p, PJ_FWD, coord);

            printf("easting: %.3f, northing: %.3f\n", a.enu.e, a.enu.n);
        }
        {
            // (Sidlerstrasse 5 - 46°57'3.9" N, 7°26'19.1" E).
            PJ_COORD coord = {{
                46.0 + 57 / 60. + 3.9 / 3600., // latitude
                7.0 + 26 / 60. + 19.1 / 3600., // longitude
                0,                             // elevation
                HUGE_VAL                       // hz...
            }};

            printf("easting: %.16f, northing: %.16f\n", coord.lpzt.lam,
                   coord.lpzt.phi);

            PJ *p_norm = proj_normalize_for_visualization(ctx, p);

            PJ_COORD a = proj_trans(p, PJ_FWD, coord);

            printf("easting: %.3f, northing: %.3f\n", a.enu.e, a.enu.n);
        }
    }

    proj_context_destroy(ctx);
}

// https://epsg.io/2056
// +proj=somerc +lat_0=46.9524055555556 +lon_0=7.43958333333333 +k_0=1
// +x_0=2600000 +y_0=1200000 +ellps=bessel
// +towgs84=674.374,15.056,405.346,0,0,0,0 +units=m +no_defs +type=crs

// EPSG:21781 - swiss coordinate system LV03
// EPSG:4978 - WGS 84 - earth centered CS
// EPSG:4326 - the WGS84 geographic coordinate system.
TEST_F(ProjF, lv95) {
    const char proj4_text[] =
        R"(+proj=somerc +lat_0=46.9524055555556 +lon_0=7.43958333333333 +k_0=1
    +x_0=2600000 +y_0=1200000 +ellps=bessel
    +towgs84=674.374,15.056,405.346,0,0,0,0 +units=m +no_defs +type=crs)";

    PJ_CONTEXT *ctx = proj_context_create();
    PJ *p = proj_create_crs_to_crs(ctx, "EPSG:4326", proj4_text, NULL);
    ASSERT_TRUE(p);

    // project same coordinates as in lv03 test
    {
        PJ_COORD coord = {{
            46.9524055555556, // latitude
            7.43958333333333, // longitude
            0,                // elevation
            HUGE_VAL          // hz...
        }};

        PJ_COORD a = proj_trans(p, PJ_FWD, coord);
        EXPECT_NEAR(2600072.390, a.enu.e, 5E-4);
        EXPECT_NEAR(1200147.056, a.enu.n, 5E-4);
    }
    {
        // (Sidlerstrasse 5 - 46°57'3.9" N, 7°26'19.1" E).
        PJ_COORD coord = {{
            46.0 + 57 / 60. + 3.9 / 3600., // latitude
            7.0 + 26 / 60. + 19.1 / 3600., // longitude
            0,                             // elevation
            HUGE_VAL                       // hz...
        }};

        printf("central (lam,phi) = %.13f, %.13f\n", coord.lpzt.lam,
               coord.lpzt.phi);

        PJ_COORD a = proj_trans(p, PJ_FWD, coord);
        EXPECT_NEAR(2600000.493, a.enu.e, 5E-4);
        EXPECT_NEAR(1200000.064, a.enu.n, 5E-4);
    }
    {
        // x_0=2600000 +y_0=1200000
        PJ_COORD xy_origin = proj_coord(2600000, 1200000, 0, 0);
        PJ_COORD a = proj_trans(p, PJ_INV, xy_origin);
        EXPECT_NEAR(46.9511, a.lp.lam, 5E-5);
        EXPECT_NEAR(7.43863, a.lp.phi, 5E-6);
    }
    proj_context_destroy(ctx);
}

/// @brief Test related to the GM-20040 case
/// @param --gtest_filter=ProjF.lv95_GM20040
TEST_F(ProjF, lv95_GM20040) {
#if 0
#include "lv95_good.h"
#else
#include "lv95_globalmapper.h"
#endif
    PJ_CONTEXT *ctx = proj_context_create();
    PJ *p = proj_create_crs_to_crs(ctx, "EPSG:4326", LV95, NULL);
    ASSERT_TRUE(p);

    {
        PROJ_STRING_LIST out_warnings{};
        PROJ_STRING_LIST out_grammar_errors{};
        PJ *p_verify = proj_create_from_wkt(ctx, LV95, NULL, &out_warnings,
                                            &out_grammar_errors);
        ASSERT_TRUE(p_verify);
        EXPECT_FALSE(out_warnings);
        EXPECT_FALSE(out_grammar_errors);
        if (out_grammar_errors) {
            for (char **lp = out_grammar_errors; *lp; ++lp)
                CONSOLE(*lp);
        }
    }

    // project same coordinates as in lv03 test
    {
        PJ_COORD coord = {{
            46.9524055555556, // latitude
            7.43958333333333, // longitude
            0,                // elevation
            HUGE_VAL          // hz...
        }};

        printf("central (lam,phi) = %.16f, %.16f\n", coord.lpzt.lam,
               coord.lpzt.phi);

        PJ_COORD a = proj_trans(p, PJ_FWD, coord);

        EXPECT_NEAR(2600072, a.enu.e, 0.5);
        EXPECT_NEAR(1200147, a.enu.n, 0.5);
    }
    {
        // (Sidlerstrasse 5 - 46°57'3.9" N, 7°26'19.1" E).
        PJ_COORD coord = {{
            46.0 + 57 / 60. + 3.9 / 3600., // latitude
            7.0 + 26 / 60. + 19.1 / 3600., // longitude
            0,                             // elevation
            HUGE_VAL                       // hz...
        }};

        printf("central (lam,phi) = %.16f, %.16f\n", coord.lpzt.lam,
               coord.lpzt.phi);

        PJ_COORD a = proj_trans(p, PJ_FWD, coord);

        EXPECT_NEAR(2600000, a.enu.e, 0.5);
        EXPECT_NEAR(1200000, a.enu.n, 0.5);
    }
    {
        // x_0=2600000 +y_0=1200000
        PJ_COORD xy_origin = proj_coord(2600000, 1200000, 0, 0);
        PJ_COORD a = proj_trans(p, PJ_INV, xy_origin);

        EXPECT_NEAR(46.95240555555556, a.lpzt.lam, 5E-3);
        EXPECT_NEAR(7.439583333333333, a.enu.n, 1E-3);
    }
    proj_context_destroy(ctx);
}

// A `parametrized` test fixture

// clang-format off
class ProjP : public testing::TestWithParam<const char *> 
{};
// clang-format on

TEST_P(ProjP, swiss_wkt) {
    const char *LV95 = GetParam();
    PJ_CONTEXT *ctx = proj_context_create();
    PJ *p = proj_create_crs_to_crs(ctx, "EPSG:4326", LV95, NULL);
    ASSERT_TRUE(p);

    {
        PROJ_STRING_LIST out_warnings{};
        PROJ_STRING_LIST out_grammar_errors{};
        PJ *p_verify = proj_create_from_wkt(ctx, LV95, NULL, &out_warnings,
                                            &out_grammar_errors);
        ASSERT_TRUE(p_verify);
        EXPECT_FALSE(out_warnings);
        EXPECT_FALSE(out_grammar_errors);
        if (out_grammar_errors) {
            for (char **lp = out_grammar_errors; *lp; ++lp)
                CONSOLE(*lp);
        }
    }

    // project same coordinates as in lv03 test
    {
        PJ_COORD coord = {{
            46.9524055555556, // latitude
            7.43958333333333, // longitude
            0,                // elevation
            HUGE_VAL          // hz...
        }};

        printf("central (lam,phi) = %.16f, %.16f\n", coord.lpzt.lam,
               coord.lpzt.phi);

        PJ_COORD a = proj_trans(p, PJ_FWD, coord);

        EXPECT_NEAR(2600072, a.enu.e, 0.5);
        EXPECT_NEAR(1200147, a.enu.n, 0.5);
    }
}

namespace {
const char *lv95_good() {
    static
#include "lv95_good.h"
        return LV95;
}

const char *lv95_globalmapper() {
    static
#include "lv95_globalmapper.h"
        return LV95;
}

const char *lv95_original() {
    static
#include "lv95_original.h"
        return LV95;
}
} // namespace

INSTANTIATE_TEST_SUITE_P(two_wkt_definitions, ProjP,
                         testing::Values(lv95_good(), lv95_globalmapper()));

//

///
TEST_F(ProjF, lv95_GM20040_2) {

    auto LV95 = lv95_original();
    CONSOLE_EVAL(LV95);

    PJ_CONTEXT *ctx = proj_context_create();
    PJ *p = proj_create_crs_to_crs(ctx, "EPSG:4326", LV95, NULL);
    ASSERT_TRUE(p);

    {
        PROJ_STRING_LIST out_warnings{};
        PROJ_STRING_LIST out_grammar_errors{};
        PJ *p_verify = proj_create_from_wkt(ctx, LV95, NULL, &out_warnings,
                                            &out_grammar_errors);
        ASSERT_TRUE(p_verify);
        EXPECT_FALSE(out_warnings);
        EXPECT_FALSE(out_grammar_errors);
        if (out_grammar_errors) {
            for (char **lp = out_grammar_errors; *lp; ++lp)
                CONSOLE(*lp);
        }
    }

    // project same coordinates as in lv03 test
    {
        PJ_COORD coord = {{
            46.9524055555556, // latitude
            7.43958333333333, // longitude
            0,                // elevation
            HUGE_VAL          // hz...
        }};

        printf("central (lam,phi) = %.16f, %.16f\n", coord.lpzt.lam,
               coord.lpzt.phi);

        PJ_COORD a = proj_trans(p, PJ_FWD, coord);

        EXPECT_NEAR(2600000, a.enu.e, 0.1E-8);
        EXPECT_NEAR(1200000, a.enu.n, 0.1E-8);
    }
}
