#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/convex_hull_2.h>
#include <vector>
#include <random>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef CGAL::Polygon_2<K> Polygon_2;

int main() {
    // Generate random points
    std::vector<Point_2> points;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 10.0);
    
    for (int i = 0; i < 20; ++i) {
        points.push_back(Point_2(dis(gen), dis(gen)));
    }
    
    std::cout << "Generated " << points.size() << " random points" << std::endl;
    
    // Compute convex hull
    Polygon_2 hull;
    CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(hull));
    
    std::cout << "Convex hull has " << hull.size() << " vertices" << std::endl;
    std::cout << "Hull area: " << hull.area() << std::endl;
    
    // Print hull vertices
    std::cout << "Hull vertices:" << std::endl;
    for (auto v = hull.vertices_begin(); v != hull.vertices_end(); ++v) {
        std::cout << "(" << v->x() << ", " << v->y() << ")" << std::endl;
    }
    
    return 0;
}