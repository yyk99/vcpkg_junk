//
//
//

#include <gtest/gtest.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Point_2.h>
#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef CGAL::Polygon_2<K> Polygon_2;

class CGAL_Test : public testing::Test
{
};

TEST_F(CGAL_Test, t0) {
    // Create a simple polygon (triangle)
    Polygon_2 triangle;
    triangle.push_back(Point_2(0, 0));
    triangle.push_back(Point_2(1, 0));
    triangle.push_back(Point_2(0.5, 1));
    
    std::cout << "Triangle vertices: " << triangle.size() << std::endl;
    std::cout << "Area: " << triangle.area() << std::endl;
    std::cout << "Is simple: " << triangle.is_simple() << std::endl;
    std::cout << "Is convex: " << triangle.is_convex() << std::endl;
    std::cout << "Orientation: " << 
        (triangle.orientation() == CGAL::CLOCKWISE ? "Clockwise" : "Counter-clockwise") 
        << std::endl;
    
    // Point-in-polygon test
    Point_2 test_point(0.3, 0.3);
    CGAL::Bounded_side side = triangle.bounded_side(test_point);
    
    std::cout << "Point (0.3, 0.3) is ";
    if (side == CGAL::ON_BOUNDED_SIDE) std::cout << "inside";
    else if (side == CGAL::ON_BOUNDARY) std::cout << "on boundary";
    else std::cout << "outside";
    std::cout << " the triangle" << std::endl;
}