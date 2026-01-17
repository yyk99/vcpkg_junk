#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <cmath>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef CGAL::Polygon_2<K> Polygon_2;

// Helper function to compute perimeter
double compute_perimeter(const Polygon_2 &poly) {
    double perimeter = 0.0;
    for (auto eit = poly.edges_begin(); eit != poly.edges_end(); ++eit) {
        perimeter += std::sqrt(CGAL::to_double(eit->squared_length()));
    }
    return perimeter;
}

int main() {
    // Create a self-intersecting polygon
    Polygon_2 poly;
    poly.push_back(Point_2(0, 0));
    poly.push_back(Point_2(2, 2));
    poly.push_back(Point_2(2, 0));
    poly.push_back(Point_2(0, 2)); // This creates self-intersection

    std::cout << "Original polygon:" << std::endl;
    std::cout << "Is simple: " << poly.is_simple() << std::endl;
    std::cout << "Is convex: " << poly.is_convex() << std::endl;
    std::cout << "Area: " << poly.area() << std::endl;

    // Create a simple polygon
    Polygon_2 simple_poly;
    simple_poly.push_back(Point_2(0, 0));
    simple_poly.push_back(Point_2(3, 0));
    simple_poly.push_back(Point_2(3, 2));
    simple_poly.push_back(Point_2(1, 3));
    simple_poly.push_back(Point_2(0, 2));

    std::cout << "\nSimple polygon:" << std::endl;
    std::cout << "Is simple: " << simple_poly.is_simple() << std::endl;
    std::cout << "Is convex: " << simple_poly.is_convex() << std::endl;
    std::cout << "Area: " << simple_poly.area() << std::endl;
    std::cout << "Perimeter: " << compute_perimeter(simple_poly) << std::endl;

    return 0;
}
