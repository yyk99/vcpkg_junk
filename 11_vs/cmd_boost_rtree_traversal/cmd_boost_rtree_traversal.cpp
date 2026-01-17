#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <iostream>
#include <vector>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

// Define types
typedef bg::model::point<double, 2, bg::cs::cartesian> point_t;
typedef bg::model::box<point_t> box_t;
typedef std::pair<point_t, int> value_t;
typedef bgi::rtree<value_t, bgi::quadratic<4>> rtree_t;

// =============================================================================
// Method 1: Iterator-based Traversal (Public API)
// =============================================================================

void traverse_with_iterators(const rtree_t& tree) {
    std::cout << "=== Iterator-based Traversal ===\n";
    
    // Forward iteration
    std::cout << "Forward iteration:\n";
    for (auto it = tree.begin(); it != tree.end(); ++it) {
        const point_t& pt = it->first;
        int data = it->second;
        std::cout << "Point: (" << bg::get<0>(pt) << ", " << bg::get<1>(pt) 
                  << ") Data: " << data << "\n";
    }
    
    // Range-based for loop
    std::cout << "\nRange-based for loop:\n";
    for (const auto& entry : tree) {
        const point_t& pt = entry.first;
        int data = entry.second;
        std::cout << "Point: (" << bg::get<0>(pt) << ", " << bg::get<1>(pt) 
                  << ") Data: " << data << "\n";
    }
}

// =============================================================================
// Method 2: Query-based Traversal
// =============================================================================

void traverse_with_queries(const rtree_t& tree) {
    std::cout << "\n=== Query-based Traversal ===\n";
    
    // Get all elements using a query that covers everything
    box_t universe(point_t(-1000, -1000), point_t(1000, 1000));
    
    std::vector<value_t> results;
    tree.query(bgi::intersects(universe), std::back_inserter(results));
    
    std::cout << "Query results (all elements):\n";
    for (const auto& entry : results) {
        const point_t& pt = entry.first;
        int data = entry.second;
        std::cout << "Point: (" << bg::get<0>(pt) << ", " << bg::get<1>(pt) 
                  << ") Data: " << data << "\n";
    }
    
    // Spatial range query
    box_t range(point_t(0, 0), point_t(5, 5));
    results.clear();
    tree.query(bgi::intersects(range), std::back_inserter(results));
    
    std::cout << "\nRange query [0,0 to 5,5]:\n";
    for (const auto& entry : results) {
        const point_t& pt = entry.first;
        int data = entry.second;
        std::cout << "Point: (" << bg::get<0>(pt) << ", " << bg::get<1>(pt) 
                  << ") Data: " << data << "\n";
    }
}

// =============================================================================
// Method 3: Visitor Pattern (for custom processing)
// =============================================================================

struct PrintVisitor {
    void operator()(const value_t& entry) const {
        const point_t& pt = entry.first;
        int data = entry.second;
        std::cout << "Visited: (" << bg::get<0>(pt) << ", " << bg::get<1>(pt) 
                  << ") Data: " << data << "\n";
    }
};

void traverse_with_visitor(const rtree_t& tree) {
    std::cout << "\n=== Visitor Pattern Traversal ===\n";
#if 0
    box_t universe(point_t(-1000, -1000), point_t(1000, 1000));
    tree.query(bgi::intersects(universe), PrintVisitor());
#else
    std::cout << "********** not implemented yet **********\n";
#endif
}

// =============================================================================
// Method 4: Internal Node Structure Access (Advanced)
// =============================================================================

// Note: This requires access to internal implementation details
// and may not work with all Boost versions. Use with caution.

template<typename Node>
struct NodeTraverser {
    static void traverse_node(const Node* node, int depth = 0) {
        if (!node) return;
        
        std::string indent(depth * 2, ' ');
        
        if (node->is_leaf()) {
            std::cout << indent << "Leaf Node (depth " << depth << "):\n";
            // Access leaf elements - this is implementation-specific
            // In practice, you'd need to cast to the actual leaf node type
        } else {
            std::cout << indent << "Internal Node (depth " << depth << "):\n";
            // Access internal node children - implementation-specific
            // Would need to cast to internal node type and iterate children
        }
    }
};

// =============================================================================
// Method 5: Custom Traversal with Predicates
// =============================================================================

void traverse_with_predicates(const rtree_t& tree) {
    std::cout << "\n=== Predicate-based Traversal ===\n";
    
    // Find all points within distance 3.0 from origin
    point_t origin(0, 0);
    std::vector<value_t> nearby;
    
    tree.query(bgi::satisfies([&origin](const value_t& entry) {
        return bg::distance(entry.first, origin) <= 3.0;
    }), std::back_inserter(nearby));
    
    std::cout << "Points within distance 3.0 from origin:\n";
    for (const auto& entry : nearby) {
        const point_t& pt = entry.first;
        double dist = bg::distance(pt, origin);
        std::cout << "Point: (" << bg::get<0>(pt) << ", " << bg::get<1>(pt) 
                  << ") Distance: " << dist << " Data: " << entry.second << "\n";
    }
}

// =============================================================================
// Method 6: Nearest Neighbor Traversal
// =============================================================================

void traverse_nearest_neighbors(const rtree_t& tree) {
    std::cout << "\n=== Nearest Neighbor Traversal ===\n";
    
    point_t query_point(2.5, 2.5);
    std::vector<value_t> nearest;
    
    // Find 3 nearest neighbors
    tree.query(bgi::nearest(query_point, 3), std::back_inserter(nearest));
    
    std::cout << "3 nearest neighbors to (2.5, 2.5):\n";
    for (size_t i = 0; i < nearest.size(); ++i) {
        const point_t& pt = nearest[i].first;
        double dist = bg::distance(pt, query_point);
        std::cout << i+1 << ". Point: (" << bg::get<0>(pt) << ", " << bg::get<1>(pt) 
                  << ") Distance: " << dist << " Data: " << nearest[i].second << "\n";
    }
}

// =============================================================================
// Method 7: Depth-First Traversal Simulation
// =============================================================================

void subdivide_and_query(const rtree_t &tree, const box_t &region, int depth);

void simulate_depth_first_traversal(const rtree_t& tree) {
    std::cout << "\n=== Simulated Depth-First Traversal ===\n";
    
    // We can simulate depth-first by using spatial queries
    // This isn't true node traversal but gives similar results
    
    auto bounds = tree.bounds();
    if (bg::is_valid(bounds)) {
        std::cout << "Tree bounds: min(" 
                  << bg::get<bg::min_corner, 0>(bounds) << ", "
                  << bg::get<bg::min_corner, 1>(bounds) << ") max("
                  << bg::get<bg::max_corner, 0>(bounds) << ", "
                  << bg::get<bg::max_corner, 1>(bounds) << ")\n";
        
        // Recursively subdivide space
        subdivide_and_query(tree, bounds, 0);
    }
}

void subdivide_and_query(const rtree_t& tree, const box_t& region, int depth) {
    if (depth > 3) return; // Limit recursion depth
    
    std::string indent(depth * 2, ' ');
    std::vector<value_t> elements;
    tree.query(bgi::intersects(region), std::back_inserter(elements));
    
    std::cout << indent << "Region depth " << depth << " contains " 
              << elements.size() << " elements\n";
    
    if (elements.size() > 1 && depth < 3) {
        // Subdivide into quadrants
        point_t min_pt = region.min_corner();
        point_t max_pt = region.max_corner();
        
        double mid_x = (bg::get<0>(min_pt) + bg::get<0>(max_pt)) / 2.0;
        double mid_y = (bg::get<1>(min_pt) + bg::get<1>(max_pt)) / 2.0;
        
        // Create four quadrants
        box_t quadrants[4] = {
            box_t(min_pt, point_t(mid_x, mid_y)),                    // SW
            box_t(point_t(mid_x, bg::get<1>(min_pt)), point_t(bg::get<0>(max_pt), mid_y)), // SE
            box_t(point_t(bg::get<0>(min_pt), mid_y), point_t(mid_x, bg::get<1>(max_pt))), // NW
            box_t(point_t(mid_x, mid_y), max_pt)                     // NE
        };
        
        for (int i = 0; i < 4; ++i) {
            subdivide_and_query(tree, quadrants[i], depth + 1);
        }
    }
}

// =============================================================================
// Main function demonstrating all traversal methods
// =============================================================================

int main() {
    // Create and populate R-tree
    rtree_t tree;
    
    // Insert sample data
    std::vector<value_t> sample_data = {
        {point_t(1, 1), 10},
        {point_t(2, 3), 20},
        {point_t(4, 2), 30},
        {point_t(0, 4), 40},
        {point_t(3, 0), 50},
        {point_t(5, 5), 60},
        {point_t(-1, 2), 70},
        {point_t(2, -1), 80}
    };
    
    for (const auto& entry : sample_data) {
        tree.insert(entry);
    }
    
    std::cout << "R-tree contains " << tree.size() << " elements\n\n";
    
    // Demonstrate all traversal methods
    traverse_with_iterators(tree);
    traverse_with_queries(tree);
    traverse_with_visitor(tree);
    traverse_with_predicates(tree);
    traverse_nearest_neighbors(tree);
    simulate_depth_first_traversal(tree);
    
    return 0;
}

// =============================================================================
// Additional Utility Functions
// =============================================================================

// Function to print tree statistics
void print_tree_stats(const rtree_t& tree) {
    std::cout << "\n=== Tree Statistics ===\n";
    std::cout << "Size: " << tree.size() << "\n";
    std::cout << "Empty: " << (tree.empty() ? "true" : "false") << "\n";
    
    if (!tree.empty()) {
        auto bounds = tree.bounds();
        std::cout << "Bounds: [" 
                  << bg::get<bg::min_corner, 0>(bounds) << ", "
                  << bg::get<bg::min_corner, 1>(bounds) << "] to ["
                  << bg::get<bg::max_corner, 0>(bounds) << ", "
                  << bg::get<bg::max_corner, 1>(bounds) << "]\n";
    }
}

// Template function for generic traversal
template<typename RTree, typename Visitor>
void generic_traverse(const RTree& tree, Visitor visitor) {
    for (const auto& entry : tree) {
        visitor(entry);
    }
}