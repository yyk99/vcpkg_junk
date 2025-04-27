//
//
//

#include <gtest/gtest.h>

#include <iostream>
#include <vector>

#include "../common/DebuggingConsole.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/polygon/voronoi.hpp>

namespace bg = boost::geometry;
typedef bg::model::d2::point_xy<double> point2d_t;
typedef bg::model::polygon<point2d_t> polygon2d_t;

TEST(boost_geometry, t0) {
    // Define a point
    point2d_t point(1.0, 1.0);

    polygon2d_t polygon;
    bg::read_wkt("POLYGON((0 0, 0 2, 2 2, 2 0, 0 0))", polygon);

    // Check if the point is within the polygon
    bool is_within = bg::within(point, polygon);
    std::cout << "Is the point within the polygon? "
              << (is_within ? "yes" : "no") << std::endl;
    ASSERT_TRUE(is_within);

    point2d_t point_outside(3.0, 3.0);
    // Check if the point is within the polygon
    bool is_within_outside = bg::within(point_outside, polygon);
    std::cout << "Is the point outside within the polygon? "
              << (is_within_outside ? "yes" : "no") << std::endl;
    ASSERT_FALSE(is_within_outside);

    // Check if a point on the border is within the polygon
    point2d_t point_on_border(0.0, 0.0);
    bool is_within_border = bg::within(point_on_border, polygon);
    std::cout << "Is the point on the border within the polygon? "
              << (is_within_border ? "yes" : "no") << std::endl;
    ASSERT_FALSE(is_within_border);
}

TEST(boost_geometry, within) {
    polygon2d_t p1;
    bg::read_wkt("POLYGON((0 0, 0 2, 2 2, 2 0, 0 0))", p1);

    polygon2d_t p2;
    bg::read_wkt("POLYGON((0 0, 0 2, 2 2, 2 0, 0 0))", p2);

    EXPECT_TRUE(bg::within(p1, p2));
    EXPECT_TRUE(bg::within(p2, p1));
}

TEST(boost_geometry, intersects) {
    polygon2d_t p1;
    bg::read_wkt("POLYGON((0 0, 0 2, 2 2, 2 0, 0 0))", p1);

    polygon2d_t p2;
    bg::read_wkt("POLYGON((0 0, 0 2, 2 2, 2 0, 0 0))", p2);

    EXPECT_TRUE(bg::intersects(p1, p2));
    EXPECT_TRUE(bg::intersects(p2, p1));
}

TEST(boost_geometry, intersects_within) {
    polygon2d_t p1;
    bg::read_wkt("POLYGON((0.5 0.5, 0.5 1.5, 1.5 1.5, 1.5 0.5, 0.5 0.5))", p1);

    polygon2d_t p2;
    bg::read_wkt("POLYGON((0 0, 0 2, 2 2, 2 0, 0 0))", p2);

    EXPECT_TRUE(bg::within(p1, p2));
    EXPECT_FALSE(bg::within(p2, p1));
    EXPECT_TRUE(bg::intersects(p1, p2));
    EXPECT_TRUE(bg::intersects(p2, p1));
}

/// @brief 
/// @param --gtest_filter=boost_geometry.weird_case_1
/// @param  
TEST(boost_geometry, weird_case_1) {
    polygon2d_t p1;
    bg::read_wkt("POLYGON((0 0, 0 10, 5 5, 10 10, 10 0, 5 5, 0 0))", p1);
    ASSERT_FALSE(bg::is_valid(p1));

    polygon2d_t p2;
    bg::read_wkt("POLYGON((0 0, 0 2, 2 2, 2 0, 0 0))", p2);
    ASSERT_TRUE(bg::is_valid(p2));

    // Split p1 into two "simple" polygons

}

#if 0
std::ostream &operator<<(std::ostream &ss, std::vector<Point> const &v) {
    ss << "{ ";
    if (v.size()) {
        ss << v.front();
        for (int i = 1; i < v.size(); ++i)
            ss << ", " << v[i];
    }
    ss << "}";
    return ss;
}
#endif