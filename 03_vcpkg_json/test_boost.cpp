//
//
//

#include <boost/filesystem.hpp>
#include <gtest/gtest.h>
#include <iostream>

TEST(t1, t1) { std::cout << "Hello...\n"; }

TEST(t1, fs1) {
    boost::filesystem::path actual(__FILE__);

    EXPECT_TRUE(boost::filesystem::exists(actual));
}

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <iostream>

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

TEST(boost_geometry, weird_case_1) {
    polygon2d_t p1;
    bg::read_wkt("POLYGON((0 0, 0 10, 5 5, 10 10, 10 0, 5 5, 0 0))", p1);

    polygon2d_t p2;
    bg::read_wkt("POLYGON((0 0, 0 2, 2 2, 2 0, 0 0))", p2);
}

//
// Voronoi tests
//

#include "boost/polygon/voronoi.hpp"
#include <iostream>
#include <vector>

class voronoy_f : public testing::Test {};

typedef boost::geometry::model::d2::point_xy<double> Point;
typedef boost::geometry::model::polygon<Point> Polygon;

struct Segment {
    Point p0;
    Point p1;
    Segment(double x1, double y1, double x2, double y2)
        : p0(x1, y1), p1(x2, y2) {}
};

namespace boost {
namespace polygon {

template <> struct geometry_concept<Point> {
    typedef point_concept type;
};

template <> struct point_traits<Point> {
    typedef int coordinate_type;

    static inline coordinate_type get(const Point &point,
                                      orientation_2d orient) {
        return (orient == HORIZONTAL) ? point.get<0>() : point.get<1>();
    }
};

template <> struct geometry_concept<Segment> {
    typedef segment_concept type;
};

template <> struct segment_traits<Segment> {
    typedef double coordinate_type;
    typedef Point point_type;

    static inline point_type get(const Segment &segment, direction_1d dir) {
        return dir.to_int() ? segment.p1 : segment.p0;
    }
};
} // namespace polygon
} // namespace boost

std::ostream &operator<<(std::ostream &ss, Point const &v) {
    ss << v.get<0>() << " " << v.get<1>();
    return ss;
}

std::ostream &operator<<(std::ostream &ss, std::vector<Point> const &v) {
    ss << "{ ";
    if (v.size()) {
        ss << v.front();
        for (int i = 1 ; i < v.size() ; ++i)
            ss << ", " << v[i];
    }
    ss << "}";
    return ss;
}

TEST_F(voronoy_f, example) {

    std::vector<Point> vertices;
    // add your input vertices

    // const char polygon_wkt[] = "POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))";
    // const char polygon_wkt[] = "POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))";
    // const char polygon_wkt[] = "POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))";

    // const char polygon_wkt[] = "POLYGON((0 0, 10 0, 5 5, 10 10, 0 10, 0 0))";

    const char polygon_wkt[] = "POLYGON((0 0, 10 0, 5 5, 10 10, 0 10, 5 5, 0 0))";

    Polygon polygon;
    boost::geometry::read_wkt(polygon_wkt, polygon);

    for (auto &p : polygon.outer())
        vertices.push_back(p);

    boost::polygon::voronoi_diagram<double> vd;
    boost::polygon::construct_voronoi(vertices.begin(), vertices.end(), &vd);

    for (const auto &vertex : vd.vertices()) {
        std::vector<Point> triangle;
        auto edge = vertex.incident_edge();
        do {
            auto cell = edge->cell();
            assert(cell->contains_point());

            triangle.push_back(vertices[cell->source_index()]);
            if (triangle.size() == 3) {
                // process output triangles
                std::cout << "Got triangle:" << triangle << std::endl;
                triangle.erase(triangle.begin() + 1);
            }

            edge = edge->rot_next();
        } while (edge != vertex.incident_edge());
    }
}
