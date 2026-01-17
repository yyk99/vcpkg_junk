#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes_2;

int main() {
    // Create outer boundary (square)
    Polygon_2 outer;
    outer.push_back(Point_2(0, 0));
    outer.push_back(Point_2(4, 0));
    outer.push_back(Point_2(4, 4));
    outer.push_back(Point_2(0, 4));
    
    // Create hole (smaller square inside)
    Polygon_2 hole;
    hole.push_back(Point_2(1, 1));
    hole.push_back(Point_2(3, 1));
    hole.push_back(Point_2(3, 3));
    hole.push_back(Point_2(1, 3));
    
    // Create polygon with holes
    std::vector<Polygon_2> holes;
    holes.push_back(hole);
    
    Polygon_with_holes_2 polygon_with_hole(outer, holes.begin(), holes.end());
    
    std::cout << "Outer boundary area: " << outer.area() << std::endl;
    std::cout << "Hole area: " << hole.area() << std::endl;
    std::cout << "Number of holes: " << polygon_with_hole.number_of_holes() << std::endl;
    
    return 0;
}