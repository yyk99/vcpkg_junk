#include <algorithm>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

// WARNING: These are internal headers and may change between versions
// Use with caution in production code
#include <boost/geometry/index/detail/rtree/utilities/view.hpp>
#include <boost/geometry/index/detail/rtree/node/node.hpp>
#include <boost/geometry/index/detail/rtree/visitors/iterator.hpp>
#include <boost/geometry/index/detail/rtree/utilities/print.hpp>

#include <chrono>
#include <iostream>
#include <memory>
#include <random>
#include <vector>

#include "boost_rtree_zond.h"

#ifndef CONSOLE
#if _DEBUG
#define CONSOLE(x)                                                             \
    do {                                                                       \
        std::cout << __func__ << ":" << x << '\n';                             \
    } while (0)
#else
#define CONSOLE(x)
#endif

#define CONSOLE_EVAL(x) CONSOLE(#x << " : " << (x))

// "test" console output
#if _DEBUG
#define CONSOLE_T(x)                                                           \
    do {                                                                       \
        std::cout << test_case_name() << "." << test_name() << ": " << x       \
                  << '\n';                                                     \
    } while (0)
#else
#define CONSOLE_T(x)
#endif

#define CONSOLE_TE(x) CONSOLE_T(#x << " : " << (x))

#endif // CONSOLE

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

// Define 2D and 3D point types
typedef bg::model::point<float, 2, bg::cs::cartesian> Point2D;
typedef bg::model::point<float, 3, bg::cs::cartesian> Point3D;
typedef bg::model::box<Point2D> Box2D;
typedef bg::model::box<Point3D> Box3D;

// Value types for storing additional data with points
typedef std::pair<Point2D, int> PointValue2D;         // Point + ID
typedef std::pair<Point3D, std::string> PointValue3D; // Point + Name
typedef std::pair<Box2D, int> BoxValue2D;             // Box + ID

// R-tree types with different configurations
typedef bgi::rtree<PointValue2D, bgi::quadratic<16>> RTree2D_Quad;
typedef bgi::rtree<PointValue3D, bgi::rstar<32>> RTree3D_RStar;
typedef bgi::rtree<BoxValue2D, bgi::linear<8>> RTreeBox2D;

typedef bgi::detail::rtree::utilities::view<RTree3D_RStar> RTree3D_RStar_view;


class BoostRTreeExamples {
public:
    // Generate random 2D points
    static std::vector<PointValue2D>
    generateRandom2DPoints(int count = 1000, float minCoord = -100.0f,
                           float maxCoord = 100.0f) {
        std::vector<PointValue2D> points;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<float> dis(minCoord, maxCoord);

        for (int i = 0; i < count; ++i) {
            Point2D point(dis(gen), dis(gen));
            points.emplace_back(point, i);
        }

        return points;
    }

    // Generate random 3D points
    static std::vector<PointValue3D>
    generateRandom3DPoints(int count = 1000, float minCoord = -50.0f,
                           float maxCoord = 50.0f) {
        std::vector<PointValue3D> points;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<float> dis(minCoord, maxCoord);

        for (int i = 0; i < count; ++i) {
            Point3D point(dis(gen), dis(gen), dis(gen));
            points.emplace_back(point, "Point_" + std::to_string(i));
        }

        return points;
    }

    // Generate random boxes
    static std::vector<BoxValue2D> generateRandomBoxes(int count = 500) {
        std::vector<BoxValue2D> boxes;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<float> pos_dis(-80.0f, 80.0f);
        std::uniform_real_distribution<float> size_dis(1.0f, 10.0f);

        for (int i = 0; i < count; ++i) {
            float x = pos_dis(gen);
            float y = pos_dis(gen);
            float w = size_dis(gen);
            float h = size_dis(gen);

            Point2D min_corner(x, y);
            Point2D max_corner(x + w, y + h);
            Box2D box(min_corner, max_corner);

            boxes.emplace_back(box, i);
        }

        return boxes;
    }

    // Example 1: Basic 2D Point R-tree Operations
    static void basic2DPointOperations() {
        std::cout << "\n=== Basic 2D Point R-tree Operations ===" << std::endl;

        // Generate test data
        auto points = generateRandom2DPoints(1000);
        std::cout << "Generated " << points.size() << " random 2D points"
                  << std::endl;

        // Create R-tree with quadratic split strategy, max 16 elements per node
        RTree2D_Quad rtree;

        // Insert points
        auto start = std::chrono::high_resolution_clock::now();
        for (const auto &point : points) {
            rtree.insert(point);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto insert_time =
            std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        std::cout << "Inserted " << rtree.size() << " points in "
                  << insert_time.count() << " μs" << std::endl;

        // Query operations
        Point2D query_point(10.0f, 20.0f);

        // 1. K-nearest neighbors
        std::cout << "\n1. K-Nearest Neighbors (k=5):" << std::endl;
        std::vector<PointValue2D> knn_result;
        rtree.query(bgi::nearest(query_point, 5),
                    std::back_inserter(knn_result));

        for (const auto &result : knn_result) {
            float x = bg::get<0>(result.first);
            float y = bg::get<1>(result.first);
            float distance = float(bg::distance(query_point, result.first));
            std::cout << "  Point " << result.second << ": (" << x << ", " << y
                      << ") distance: " << distance << std::endl;
        }

        // 2. Range query (within box)
        Box2D query_box(Point2D(0.0f, 0.0f), Point2D(30.0f, 30.0f));
        std::cout << "\n2. Range Query (box: 0,0 to 30,30):" << std::endl;

        std::vector<PointValue2D> range_result;
        rtree.query(bgi::within(query_box), std::back_inserter(range_result));
        std::cout << "  Found " << range_result.size() << " points within box"
                  << std::endl;

        for (size_t i = 0; i < std::min(size_t(5), range_result.size()); ++i) {
            float x = bg::get<0>(range_result[i].first);
            float y = bg::get<1>(range_result[i].first);
            std::cout << "    Point " << range_result[i].second << ": (" << x
                      << ", " << y << ")" << std::endl;
        }

        // 3. Intersects query
        std::cout << "\n3. Intersects Query (same box):" << std::endl;
        std::vector<PointValue2D> intersects_result;
        rtree.query(bgi::intersects(query_box),
                    std::back_inserter(intersects_result));
        std::cout << "  Found " << intersects_result.size()
                  << " points intersecting box" << std::endl;

        // 4. Remove some points
        std::cout << "\n4. Removal Operations:" << std::endl;
        size_t initial_size = rtree.size();
        size_t removed = rtree.remove(
            range_result.begin(),
            range_result.begin() + std::min(size_t(10), range_result.size()));
        std::cout << "  Removed " << removed << " points" << std::endl;
        std::cout << "  R-tree size: " << initial_size << " -> " << rtree.size()
                  << std::endl;
    }
    // Example 2: 3D Point R-tree with R* algorithm
    static void rstar3DOperations() {
        std::cout << "\n=== 3D Point R-tree with R* Algorithm ===" << std::endl;

        auto points = generateRandom3DPoints(2000);
        std::cout << "Generated " << points.size() << " random 3D points"
                  << std::endl;

        // R* tree generally provides better query performance but slower
        // insertions
        RTree3D_RStar rtree;

        auto start = std::chrono::high_resolution_clock::now();
        rtree.insert(points.begin(), points.end());
        auto end = std::chrono::high_resolution_clock::now();
        auto insert_time =
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        std::cout << "Bulk inserted " << rtree.size() << " points in "
                  << insert_time.count() << " ms" << std::endl;

        Point3D query_point(5.0f, 10.0f, -5.0f);

        // Nearest neighbor search
        std::cout << "\nNearest neighbors to point (5, 10, -5):" << std::endl;
        std::vector<PointValue3D> nn_result;
        rtree.query(bgi::nearest(query_point, 3),
                    std::back_inserter(nn_result));

        for (const auto &result : nn_result) {
            float x = bg::get<0>(result.first);
            float y = bg::get<1>(result.first);
            float z = bg::get<2>(result.first);
            float distance = float(bg::distance(query_point, result.first));
            std::cout << "  " << result.second << ": (" << x << ", " << y
                      << ", " << z << ") distance: " << distance << std::endl;
        }

        // 3D box query
        Box3D query_box(Point3D(-10.0f, -10.0f, -10.0f),
                        Point3D(10.0f, 10.0f, 10.0f));
        std::cout << "\n3D Box query (-10,-10,-10) to (10,10,10):" << std::endl;

        size_t count = 0;
#if 0
        rtree.query(bgi::within(query_box), [&count](const PointValue3D &v) {
            if (count < 5) {
                float x = bg::get<0>(v.first);
                float y = bg::get<1>(v.first);
                float z = bg::get<2>(v.first);
                std::cout << "  " << v.second << ": (" << x << ", " << y << ", "
                          << z << ")" << std::endl;
            }
            count++;
            return true; // Continue iteration
        });
#else
        // std::back_inserter(result)
        std::vector<PointValue3D> result;
        rtree.query(bgi::within(query_box), std::back_inserter(result));
        count = result.size();
#endif
        std::cout << "  Total points in box: " << count << std::endl;

        if(0) {
            auto traverse_with_iterators = [](RTree3D_RStar const &tree) {
                std::cout << "=== Iterator-based Traversal ===\n";

                // Forward iteration
                std::cout << "Forward iteration:\n";
                for (auto it = tree.begin(); it != tree.end(); ++it) {
                    auto const &pt = it->first;
                    auto data = it->second;
                    std::cout << "Point: (" << bg::get<0>(pt) << ", "
                              << bg::get<1>(pt) << ") Data: " << data << "\n";
                }

                // Range-based for loop
                std::cout << "\nRange-based for loop:\n";
                for (const auto &entry : tree) {
                    auto const &pt = entry.first;
                    auto data = entry.second;
                    std::cout << "Point: (" << bg::get<0>(pt) << ", "
                              << bg::get<1>(pt) << ") Data: " << data << "\n";
                }
            };

            traverse_with_iterators(rtree);
        }
        {
            struct zond_t : public boost::static_visitor<> {
                void operator()(
                    bgi::detail::rtree::variant_leaf<PointValue3D, bgi::rstar<32>, Point3D, bgi::detail::rtree::allocators<boost::container::new_allocator<PointValue3D>, PointValue3D, bgi::rstar<32>, Point3D, bgi::detail::rtree::node_variant_static_tag>, bgi::detail::rtree::node_variant_static_tag> const &hz
                    ) 
                { 
                    std::cout << "here...\n";
                };
            };
            RTree3D_RStar_view backdoor(rtree);

            zond_t zond;
            backdoor.apply_visitor(zond);

            //auto root = rtree.root;
        }
    }

    static void rstar3DOperationsToo() {
        std::cout << "\n=== 3D Point R-tree with R* Algorithm (Too) ===" << std::endl;

        auto points = generateRandom3DPoints(100, 0, 100);
        std::cout << "Generated " << points.size() << " random 3D points"
                  << std::endl;

        // R* tree generally provides better query performance but slower
        // insertions
        RTree3D_RStar rtree;

        rtree.insert(points.begin(), points.end());

        Point3D query_point(5.0f, 10.0f, -5.0f);

        struct zond_t {
            void
            operator()(bgi::detail::rtree::variant_leaf<
                       PointValue3D, bgi::rstar<32>, Point3D,
                       bgi::detail::rtree::allocators<
                           boost::container::new_allocator<PointValue3D>,
                           PointValue3D, bgi::rstar<32>, Point3D,
                           bgi::detail::rtree::node_variant_static_tag>,
                       bgi::detail::rtree::node_variant_static_tag> const &hz) {
                std::cout << "here...\n";
            };
        };
        RTree3D_RStar_view backdoor(rtree);
        zond_t zond;
        backdoor.apply_visitor(zond);

        {
            auto depth = backdoor.depth();
            std::cout << "depth == " << depth << "\n";
        }

        {
            auto tr = backdoor.translator();
            std::cout << "tr == " << sizeof(tr) << "\n";
        }
        if (0) {
            CONSOLE("bgi::detail::rtree::utilities::print...");
            bgi::detail::rtree::utilities::print(std::cout, rtree);
        }
        {
            CONSOLE("bgi::detail::rtree::utilities::zond...");
            bgi::detail::rtree::utilities::zond(std::cout, rtree);
        }
    }

    // Example 3: Box R-tree for Rectangular Objects
    static void boxRTreeOperations() {
        std::cout << "\n=== Box R-tree Operations ===" << std::endl;

        auto boxes = generateRandomBoxes(500);
        std::cout << "Generated " << boxes.size() << " random boxes"
                  << std::endl;

        RTreeBox2D rtree(boxes.begin(), boxes.end());
        std::cout << "Created R-tree with " << rtree.size() << " boxes"
                  << std::endl;

        // Query operations with boxes
        Box2D query_box(Point2D(20.0f, 20.0f), Point2D(40.0f, 40.0f));

        // 1. Find boxes that intersect with query box
        std::cout
            << "\n1. Boxes intersecting with query box (20,20) to (40,40):"
            << std::endl;
        std::vector<BoxValue2D> intersecting_boxes;
        rtree.query(bgi::intersects(query_box),
                    std::back_inserter(intersecting_boxes));

        std::cout << "  Found " << intersecting_boxes.size()
                  << " intersecting boxes" << std::endl;
        for (size_t i = 0; i < std::min(size_t(5), intersecting_boxes.size());
             ++i) {
            const Box2D &box = intersecting_boxes[i].first;
            float min_x = bg::get<bg::min_corner, 0>(box);
            float min_y = bg::get<bg::min_corner, 1>(box);
            float max_x = bg::get<bg::max_corner, 0>(box);
            float max_y = bg::get<bg::max_corner, 1>(box);

            std::cout << "    Box " << intersecting_boxes[i].second << ": ("
                      << min_x << "," << min_y << ") to (" << max_x << ","
                      << max_y << ")" << std::endl;
        }

        // 2. Find boxes completely within query box
        std::cout << "\n2. Boxes completely within query box:" << std::endl;
        std::vector<BoxValue2D> within_boxes;
        rtree.query(bgi::within(query_box), std::back_inserter(within_boxes));
        std::cout << "  Found " << within_boxes.size()
                  << " boxes completely within query area" << std::endl;

        // 3. Find boxes that contain a specific point
        Point2D test_point(25.0f, 25.0f);
        std::cout << "\n3. Boxes containing point (25, 25):" << std::endl;
        std::vector<BoxValue2D> containing_boxes;
        rtree.query(bgi::contains(test_point),
                    std::back_inserter(containing_boxes));

        std::cout << "  Found " << containing_boxes.size()
                  << " boxes containing the point" << std::endl;
        for (const auto &box_value : containing_boxes) {
            std::cout << "    Box " << box_value.second << " contains the point"
                      << std::endl;
        }
    }

    // Example 4: Custom Predicates and Advanced Queries
    static void advancedQueryOperations() {
        std::cout << "\n=== Advanced Query Operations ===" << std::endl;

        auto points = generateRandom2DPoints(1000);
        RTree2D_Quad rtree(points.begin(), points.end());

        Point2D center(0.0f, 0.0f);

        // Custom predicate: find points within certain distance
        auto distance_predicate = [&center](const PointValue2D &v) {
            return bg::distance(center, v.first) <= 20.0f;
        };

        std::cout
            << "\nCustom predicate - Points within distance 20 from origin:"
            << std::endl;
        std::vector<PointValue2D> custom_result;
        rtree.query(bgi::satisfies(distance_predicate),
                    std::back_inserter(custom_result));
        std::cout << "  Found " << custom_result.size() << " points"
                  << std::endl;

        // Combine predicates
        Box2D region(Point2D(-30.0f, -30.0f), Point2D(30.0f, 30.0f));

        std::cout
            << "\nCombined predicates - Points in box AND within distance 15:"
            << std::endl;
        std::vector<PointValue2D> combined_result;
        rtree.query(bgi::within(region) &&
                        bgi::satisfies([&center](const PointValue2D &v) {
                            return bg::distance(center, v.first) <= 15.0f;
                        }),
                    std::back_inserter(combined_result));

        std::cout << "  Found " << combined_result.size()
                  << " points matching both conditions" << std::endl;

        // Iterative query with early termination
        std::cout << "\nIterative query with early termination (find 3 closest "
                     "points):"
                  << std::endl;
        int found = 0;
        for (auto it = rtree.qbegin(bgi::nearest(center, 10));
             it != rtree.qend() && found < 3; ++it, ++found) {
            float x = bg::get<0>(it->first);
            float y = bg::get<1>(it->first);
            float distance = float(bg::distance(center, it->first));
            std::cout << "  Point " << it->second << ": (" << x << ", " << y
                      << ") distance: " << distance << std::endl;
        }
    }

    // Example 5: Performance Comparison
    static void performanceComparison() {
        std::cout << "\n=== Performance Comparison ===" << std::endl;

        const int POINT_COUNT = 10000;
        const int QUERY_COUNT = 1000;

        auto points = generateRandom2DPoints(POINT_COUNT);
        auto query_points = generateRandom2DPoints(QUERY_COUNT, -50.0f, 50.0f);

        // Build R-tree
        auto start = std::chrono::high_resolution_clock::now();
        RTree2D_Quad rtree(points.begin(), points.end());
        auto build_time = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start);

        std::cout << "Built R-tree with " << POINT_COUNT << " points in "
                  << build_time.count() << " ms" << std::endl;

        // R-tree nearest neighbor queries
        start = std::chrono::high_resolution_clock::now();
        for (const auto &query_point : query_points) {
            std::vector<PointValue2D> result;
            rtree.query(bgi::nearest(query_point.first, 5),
                        std::back_inserter(result));
        }
        auto rtree_time = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start);

        // Brute force comparison
        start = std::chrono::high_resolution_clock::now();
        for (const auto &query_point : query_points) {
            std::vector<std::pair<float, int>> distances;
            for (const auto &point : points) {
                float dist = float(bg::distance(query_point.first, point.first));
                distances.emplace_back(dist, point.second);
            }
            std::partial_sort(distances.begin(), distances.begin() + 5,
                              distances.end());
        }
        auto brute_force_time =
            std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - start);

        std::cout << "\nPerformance Results (" << QUERY_COUNT
                  << " queries):" << std::endl;
        std::cout << "  R-tree queries: " << rtree_time.count() << " ms"
                  << std::endl;
        std::cout << "  Brute force: " << brute_force_time.count() << " ms"
                  << std::endl;
        std::cout << "  Speedup: "
                  << float(brute_force_time.count()) / rtree_time.count() << "x"
                  << std::endl;
    }

    // Example 6: Dynamic Operations
    static void dynamicOperations() {
        std::cout << "\n=== Dynamic Operations ===" << std::endl;

        RTree2D_Quad rtree;

        // Progressive insertion
        std::cout << "Progressive insertion:" << std::endl;
        for (int i = 0; i < 100; ++i) {
            Point2D point(i * 0.5f, i * 0.3f);
            rtree.insert(std::make_pair(point, i));

            if (i % 20 == 0) {
                std::cout << "  Inserted " << (i + 1)
                          << " points, tree size: " << rtree.size()
                          << std::endl;
            }
        }

        // Conditional removal
        std::cout << "\nConditional removal (remove points with x > 30):"
                  << std::endl;
        size_t removed = 0;
        auto it = rtree.qbegin(bgi::intersects(
            Box2D(Point2D(30.0f, -1000.0f), Point2D(1000.0f, 1000.0f))));

        std::vector<PointValue2D> to_remove;
        for (; it != rtree.qend(); ++it) {
            to_remove.push_back(*it);
        }

        for (const auto &point : to_remove) {
            rtree.remove(point);
            removed++;
        }

        std::cout << "  Removed " << removed << " points" << std::endl;
        std::cout << "  Final tree size: " << rtree.size() << std::endl;

        // Clear and rebuild
        std::cout << "\nClear and rebuild:" << std::endl;
        rtree.clear();
        std::cout << "  Cleared tree, size: " << rtree.size() << std::endl;

        auto new_points = generateRandom2DPoints(50);
        rtree.insert(new_points.begin(), new_points.end());
        std::cout << "  Rebuilt with " << rtree.size() << " points"
                  << std::endl;
    }

};

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include "windows.h"
#endif

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
    std::cout << "code page == " << GetACP() << "\n";
#endif
    std::cout << "Boost.Geometry R-tree Examples" << std::endl;
    std::cout << "==============================" << std::endl;

    try {
        BoostRTreeExamples::basic2DPointOperations();
        BoostRTreeExamples::rstar3DOperations();
        BoostRTreeExamples::rstar3DOperationsToo();
        BoostRTreeExamples::boxRTreeOperations();
        BoostRTreeExamples::advancedQueryOperations();
        BoostRTreeExamples::performanceComparison();
        BoostRTreeExamples::dynamicOperations();

        std::cout << "\n=== R-tree Configuration Options ===" << std::endl;
        std::cout << "Available splitting algorithms:" << std::endl;
        std::cout << "  - Linear: Fast insertion, moderate query performance"
                  << std::endl;
        std::cout << "  - Quadratic: Balanced insertion/query performance"
                  << std::endl;
        std::cout << "  - R*: Slower insertion, best query performance"
                  << std::endl;
        std::cout << "\nNode capacity affects memory usage and performance"
                  << std::endl;
        std::cout << "Higher capacity = fewer nodes, better cache performance"
                  << std::endl;

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return -1;
    }

    return 0;
}

// Build instructions:
/*
Compilation requires Boost libraries:

g++ -std=c++14 -O3 boost_rtree_example.cpp -o boost_rtree_example

Or with CMake:

cmake_minimum_required(VERSION 3.10)
project(boost_rtree_example)

find_package(Boost REQUIRED)

add_executable(boost_rtree_example boost_rtree_example.cpp)
target_link_libraries(boost_rtree_example ${Boost_LIBRARIES})
target_include_directories(boost_rtree_example PRIVATE ${Boost_INCLUDE_DIRS})

set_property(TARGET boost_rtree_example PROPERTY CXX_STANDARD 14)
*/