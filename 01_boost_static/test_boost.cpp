//
//
//

#include <gtest/gtest.h>

#include <boost/filesystem.hpp>
#include <iostream>
#include <vector>

#include "../common/CONSOLE.h"

TEST(t1, t1) { std::cout << "Hello...\n"; }

TEST(t1, fs1) {
    boost::filesystem::path actual(__FILE__);

    EXPECT_TRUE(boost::filesystem::exists(actual));
}

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

//
// Voronoi tests
//

class voronoy_f : public testing::Test {};

typedef boost::geometry::model::d2::point_xy<double> Point;
typedef boost::geometry::model::polygon<Point> Polygon;
typedef boost::geometry::model::segment<Point> Segment;

namespace boost {
namespace polygon {

template <> struct geometry_concept<Point> {
    typedef point_concept type;
};

template <> struct point_traits<Point> {
    typedef double coordinate_type;

    static inline coordinate_type get(const Point &point,
                                      orientation_2d orient) {
        return (orient == HORIZONTAL) ? point.get<0>() : point.get<1>();
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
        for (int i = 1; i < v.size(); ++i)
            ss << ", " << v[i];
    }
    ss << "}";
    return ss;
}

TEST_F(voronoy_f, example) {

    std::vector<Point> vertices;
    std::vector<Polygon> triangulation;
    // add your input vertices

    // const char polygon_wkt[] = "POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))";
    // const char polygon_wkt[] = "POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))";
    // const char polygon_wkt[] = "POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))";

    // const char polygon_wkt[] = "POLYGON((0 0, 10 0, 5 5, 10 10, 0 10, 0 0))";

    const char polygon_wkt[] =
        "POLYGON((0 0, 10 0, 5 5, 10 10, 0 10, 5 5, 0 0))";

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
                Polygon tr;
                tr.outer().push_back(triangle[0]);
                tr.outer().push_back(triangle[1]);
                tr.outer().push_back(triangle[2]);
                triangulation.push_back(tr);

                CONSOLE("Got triangle:" << triangle);
                triangle.erase(triangle.begin() + 1);
            }

            edge = edge->rot_next();
        } while (edge != vertex.incident_edge());
    }

    CONSOLE("...");
    ASSERT_EQ(4, triangulation.size());
}

// Boost.Polygon library voronoi_basic_tutorial.cpp file

//          Copyright Andrii Sydorchuk 2010-2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

// See http://www.boost.org for updates, documentation, and revision history.

namespace boost {
namespace polygon {

// template <> struct geometry_concept<Point> {
//     typedef point_concept type;
// };
//
// template <> struct point_traits<Point> {
//     typedef double coordinate_type;
//
//     static inline coordinate_type get(const Point &point,
//                                       orientation_2d orient) {
//         return (orient == HORIZONTAL) ? point.get<0>() : point.get<1>();
//     }
// };

template <> struct geometry_concept<Segment> {
    typedef segment_concept type;
};

template <> struct segment_traits<Segment> {
    typedef double coordinate_type;
    typedef Point point_type;

    static inline point_type get(const Segment &segment, direction_1d dir) {
        return dir.to_int() ? segment.second : segment.first;
    }
};
} // namespace polygon
} // namespace boost

using boost::polygon::high;
using boost::polygon::low;
using boost::polygon::voronoi_builder;
using boost::polygon::voronoi_diagram;
using boost::polygon::x;
using boost::polygon::y;

class voronoi_basic_tutorial_f : public testing::Test {
public:
    // Traversing Voronoi edges using edge iterator.
    int iterate_primary_edges1(const voronoi_diagram<double> &vd) {
        int result = 0;
        for (voronoi_diagram<double>::const_edge_iterator it =
                 vd.edges().begin();
             it != vd.edges().end(); ++it) {
            if (it->is_primary())
                ++result;
        }
        return result;
    }

    // Traversing Voronoi edges using cell iterator.
    int iterate_primary_edges2(const voronoi_diagram<double> &vd) {
        int result = 0;
        for (voronoi_diagram<double>::const_cell_iterator it =
                 vd.cells().begin();
             it != vd.cells().end(); ++it) {
            const voronoi_diagram<double>::cell_type &cell = *it;
            const voronoi_diagram<double>::edge_type *edge =
                cell.incident_edge();
            // This is convenient way to iterate edges around Voronoi cell.
            do {
                if (edge->is_primary())
                    ++result;
                edge = edge->next();
            } while (edge != cell.incident_edge());
        }
        return result;
    }

    // Traversing Voronoi edges using vertex iterator.
    // As opposite to the above two functions this one will not iterate through
    // edges without finite endpoints and will iterate only once through edges
    // with single finite endpoint.
    int iterate_primary_edges3(const voronoi_diagram<double> &vd) {
        int result = 0;
        for (voronoi_diagram<double>::const_vertex_iterator it =
                 vd.vertices().begin();
             it != vd.vertices().end(); ++it) {
            const voronoi_diagram<double>::vertex_type &vertex = *it;
            const voronoi_diagram<double>::edge_type *edge =
                vertex.incident_edge();
            // This is convenient way to iterate edges around Voronoi vertex.
            do {
                if (edge->is_primary())
                    ++result;
                edge = edge->rot_next();
            } while (edge != vertex.incident_edge());
        }
        return result;
    }
};

/// @brief Converted voronoi_basic_tutorial example
/// @param --gtest_filter=voronoi_basic_tutorial_f.main
/// @param
TEST_F(voronoi_basic_tutorial_f, main) {
    // Preparing Input Geometries.
    std::vector<Point> points;
    points.push_back(Point(0, 0));
    points.push_back(Point(1, 6));
    std::vector<Segment> segments;
    segments.push_back(Segment({-4, 5}, {5, -1}));
    segments.push_back(Segment({3, -11}, {13, -1}));

    // Construction of the Voronoi Diagram.
    voronoi_diagram<double> vd;
    construct_voronoi(points.begin(), points.end(), segments.begin(),
                      segments.end(), &vd);

    // Traversing Voronoi Graph.
    {
        printf("Traversing Voronoi graph.\n");
        printf("Number of visited primary edges using edge iterator: %d\n",
               iterate_primary_edges1(vd));
        ASSERT_EQ(24, iterate_primary_edges1(vd));
        printf("Number of visited primary edges using cell iterator: %d\n",
               iterate_primary_edges2(vd));
        ASSERT_EQ(24, iterate_primary_edges2(vd));
        printf("Number of visited primary edges using vertex iterator: %d\n",
               iterate_primary_edges3(vd));
        ASSERT_EQ(21, iterate_primary_edges3(vd));
        printf("\n");
    }

    ASSERT_EQ(8, vd.cells().size());

    // Using color member of the Voronoi primitives to store the average number
    // of edges around each cell (including secondary edges).
    {
        printf("Number of edges (including secondary) around the Voronoi "
               "cells:\n");
        for (voronoi_diagram<double>::const_edge_iterator it =
                 vd.edges().begin();
             it != vd.edges().end(); ++it) {
            std::size_t cnt = it->cell()->color();
            it->cell()->color(cnt + 1);
        }
        for (voronoi_diagram<double>::const_cell_iterator it =
                 vd.cells().begin();
             it != vd.cells().end(); ++it) {
            printf("%lu ", (unsigned long)it->color());
        }
        printf("\n");
        printf("\n");
        ASSERT_EQ(4, vd.cells()[0].color());
        ASSERT_EQ(6, vd.cells()[1].color());
        ASSERT_EQ(4, vd.cells()[2].color());
        ASSERT_EQ(3, vd.cells()[3].color());
        ASSERT_EQ(3, vd.cells()[4].color());
        ASSERT_EQ(5, vd.cells()[5].color());
        ASSERT_EQ(3, vd.cells()[6].color());
        ASSERT_EQ(4, vd.cells()[7].color());
    }

    // Linking Voronoi cells with input geometries.
    {
        unsigned int cell_index = 0;
        for (voronoi_diagram<double>::const_cell_iterator it =
                 vd.cells().begin();
             it != vd.cells().end(); ++it) {
            if (it->contains_point()) {
                if (it->source_category() ==
                    boost::polygon::SOURCE_CATEGORY_SINGLE_POINT) {
                    std::size_t index = it->source_index();
                    Point p = points[index];
                    printf("Cell #%u contains a point: (%g, %g).\n", cell_index,
                           x(p), y(p));
                } else if (it->source_category() ==
                           boost::polygon::
                               SOURCE_CATEGORY_SEGMENT_START_POINT) {
                    std::size_t index = it->source_index() - points.size();
                    Point p0 = low(segments[index]);
                    printf("Cell #%u contains segment start point: (%g, %g).\n",
                           cell_index, x(p0), y(p0));
                } else if (it->source_category() ==
                           boost::polygon::SOURCE_CATEGORY_SEGMENT_END_POINT) {
                    std::size_t index = it->source_index() - points.size();
                    Point p1 = high(segments[index]);
                    printf("Cell #%u contains segment end point: (%g, %g).\n",
                           cell_index, x(p1), y(p1));
                }
            } else {
                std::size_t index = it->source_index() - points.size();
                Point p0 = low(segments[index]);
                Point p1 = high(segments[index]);
                printf("Cell #%u contains a segment: ((%g, %g), (%g, %g)). \n",
                       cell_index, x(p0), y(p0), x(p1), y(p1));
            }
            ++cell_index;
        }
    }
}
