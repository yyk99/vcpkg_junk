#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes_2;

int main() {
    // Create first polygon (square)
    Polygon_2 poly1;
    poly1.push_back(Point_2(0, 0));
    poly1.push_back(Point_2(2, 0));
    poly1.push_back(Point_2(2, 2));
    poly1.push_back(Point_2(0, 2));

    // Create second polygon (overlapping square)
    Polygon_2 poly2;
    poly2.push_back(Point_2(1, 1));
    poly2.push_back(Point_2(3, 1));
    poly2.push_back(Point_2(3, 3));
    poly2.push_back(Point_2(1, 3));

    // Union - returns bool and takes single output polygon
    Polygon_with_holes_2 union_result;
    if (CGAL::join(poly1, poly2, union_result)) {
        std::cout << "Union succeeded" << std::endl;
        std::cout << "  Union has " << union_result.number_of_holes() << " holes" << std::endl;
    } else {
        std::cout << "Union failed (polygons are disjoint)" << std::endl;
    }

    // Intersection - can produce multiple polygons, uses output iterator
    std::vector<Polygon_with_holes_2> intersection_result;
    CGAL::intersection(poly1, poly2, std::back_inserter(intersection_result));
    std::cout << "Intersection result has " << intersection_result.size() << " polygon(s)" << std::endl;
    if (!intersection_result.empty()) {
        std::cout << "  First polygon has " << intersection_result[0].number_of_holes() << " holes" << std::endl;
    }

    // Difference - can produce multiple polygons, uses output iterator
    std::vector<Polygon_with_holes_2> difference_result;
    CGAL::difference(poly1, poly2, std::back_inserter(difference_result));
    std::cout << "Difference result has " << difference_result.size() << " polygon(s)" << std::endl;
    if (!difference_result.empty()) {
        std::cout << "  First polygon has " << difference_result[0].number_of_holes() << " holes" << std::endl;
    }

    return 0;
}
