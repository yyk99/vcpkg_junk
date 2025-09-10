#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <iostream>
#include <vector>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

// Define 2D point type
typedef bg::model::point<float, 2, bg::cs::cartesian> Point;
typedef bg::model::box<Point> Box;
typedef std::pair<Point, int> Value; // Point with ID

// R* tree with maximum 16 values per node
typedef bgi::rtree<Value, bgi::rstar<16>> RStarTree;

#ifdef _WIN32
#   define WIN32_LEAN_AND_MEAN
#   include "windows.h"
#endif

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
    std::cout << "code page == " << GetACP() << "\n";
#endif
    
    std::cout << "Simple Boost R* Tree Example\n";
    std::cout << "============================\n\n";

    // Create R* tree
    RStarTree rtree;

    // Sample data points
    std::vector<Value> points = {{Point(2, 3), 1}, {Point(5, 4), 2},
                                 {Point(9, 6), 3}, {Point(4, 7), 4},
                                 {Point(8, 1), 5}, {Point(7, 2), 6},
                                 {Point(1, 8), 7}, {Point(6, 5), 8},
                                 {Point(3, 9), 9}, {Point(10, 10), 10}};

    // Insert points into R* tree
    std::cout << "1. Inserting points:\n";
    for (const auto &point : points) {
        rtree.insert(point);
        std::cout << "   Point " << point.second << ": ("
                  << bg::get<0>(point.first) << ", " << bg::get<1>(point.first)
                  << ")\n";
    }
    std::cout << "   Total points inserted: " << rtree.size() << "\n\n";

    // 2. Range Query - Find points in a box
    std::cout << "2. Range Query - Points in box (3,3) to (8,8):\n";
    Box query_box(Point(3, 3), Point(8, 8));

    std::vector<Value> result_within;
    rtree.query(bgi::within(query_box), std::back_inserter(result_within));

    std::cout << "   Found " << result_within.size() << " points:\n";
    for (const auto &val : result_within) {
        std::cout << "   Point " << val.second << ": (" << bg::get<0>(val.first)
                  << ", " << bg::get<1>(val.first) << ")\n";
    }
    std::cout << "\n";

    // 3. K-Nearest Neighbors
    std::cout << "3. K-Nearest Neighbors - 3 closest to (5, 5):\n";
    Point query_point(5, 5);

    std::vector<Value> result_nearest;
    rtree.query(bgi::nearest(query_point, 3),
                std::back_inserter(result_nearest));

    for (const auto &val : result_nearest) {
        float distance = static_cast<float>(bg::distance(query_point, val.first));
        std::cout << "   Point " << val.second << ": (" << bg::get<0>(val.first)
                  << ", " << bg::get<1>(val.first)
                  << ") - distance: " << distance << "\n";
    }
    std::cout << "\n";

    // 4. Intersects Query - Points that intersect with a box
    std::cout
        << "4. Intersects Query - Points intersecting box (0,0) to (4,4):\n";
    Box intersect_box(Point(0, 0), Point(4, 4));

    std::vector<Value> result_intersects;
    rtree.query(bgi::intersects(intersect_box),
                std::back_inserter(result_intersects));

    std::cout << "   Found " << result_intersects.size() << " points:\n";
    for (const auto &val : result_intersects) {
        std::cout << "   Point " << val.second << ": (" << bg::get<0>(val.first)
                  << ", " << bg::get<1>(val.first) << ")\n";
    }
    std::cout << "\n";

    // 5. Custom Query - Points with x > 6
    std::cout << "5. Custom Query - Points with x coordinate > 6:\n";

    std::vector<Value> result_custom;
    rtree.query(bgi::satisfies(
                    [](const Value &v) { return bg::get<0>(v.first) > 6.0f; }),
                std::back_inserter(result_custom));

    std::cout << "   Found " << result_custom.size() << " points:\n";
    for (const auto &val : result_custom) {
        std::cout << "   Point " << val.second << ": (" << bg::get<0>(val.first)
                  << ", " << bg::get<1>(val.first) << ")\n";
    }
    std::cout << "\n";

    // 6. Iterator-based query (early termination)
    std::cout << "6. Iterator Query - First 2 nearest points to (1, 1):\n";
    Point query_point2(1, 1);

    int count = 0;
    for (auto it = rtree.qbegin(bgi::nearest(query_point2, 5));
         it != rtree.qend() && count < 2; ++it, ++count) {

        float distance = static_cast<float>(bg::distance(query_point2, it->first));
        std::cout << "   Point " << it->second << ": (" << bg::get<0>(it->first)
                  << ", " << bg::get<1>(it->first)
                  << ") - distance: " << distance << "\n";
    }
    std::cout << "\n";

    // 7. Remove points
    std::cout << "7. Removal - Remove points with ID 3 and 7:\n";

    // Find and remove point with ID 3
    auto to_remove = std::find_if(points.begin(), points.end(),
                                  [](const Value &v) { return v.second == 3; });
    if (to_remove != points.end()) {
        rtree.remove(*to_remove);
        std::cout << "   Removed point " << to_remove->second << "\n";
    }

    // Find and remove point with ID 7
    to_remove = std::find_if(points.begin(), points.end(),
                             [](const Value &v) { return v.second == 7; });
    if (to_remove != points.end()) {
        rtree.remove(*to_remove);
        std::cout << "   Removed point " << to_remove->second << "\n";
    }

    std::cout << "   Points remaining: " << rtree.size() << "\n\n";

    // 8. Final query to verify removal
    std::cout << "8. Final Query - All remaining points:\n";
    std::vector<Value> all_remaining;
    rtree.query(bgi::satisfies([](const Value &) { return true; }),
                std::back_inserter(all_remaining));

    std::cout << "   Remaining " << all_remaining.size() << " points:\n";
    for (const auto &val : all_remaining) {
        std::cout << "   Point " << val.second << ": (" << bg::get<0>(val.first)
                  << ", " << bg::get<1>(val.first) << ")\n";
    }

    return 0;
}

/*
Compilation:
g++ -std=c++14 boost_rstar_simple.cpp -o boost_rstar_simple

Requirements:
- Boost libraries (header-only for Boost.Geometry)
- C++14 or later

Key R* Features Demonstrated:
1. Optimal splitting algorithm for better query performance
2. Range queries (within/intersects)
3. K-nearest neighbor searches
4. Custom predicates
5. Iterator-based queries with early termination
6. Dynamic insertion and removal
7. Efficient spatial indexing

R* Algorithm Benefits:
- Better query performance than basic R-tree
- Minimizes overlap and coverage
- Good for read-heavy workloads
- Automatic optimization during insertion
*/