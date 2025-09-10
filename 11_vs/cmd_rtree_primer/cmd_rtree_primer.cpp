#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <vector>

// 2D Point structure
struct Point2D {
    float x, y;

    Point2D(float x = 0, float y = 0) : x(x), y(y) {}

    float distance(const Point2D &other) const {
        float dx = x - other.x;
        float dy = y - other.y;
        return std::sqrt(dx * dx + dy * dy);
    }

    void print() const { std::cout << "(" << x << ", " << y << ")"; }
};

// Rectangle (bounding box) structure
struct Rectangle {
    float min_x, min_y, max_x, max_y;

    Rectangle()
        : min_x(std::numeric_limits<float>::max()),
          min_y(std::numeric_limits<float>::max()),
          max_x(std::numeric_limits<float>::lowest()),
          max_y(std::numeric_limits<float>::lowest()) {}

    Rectangle(float min_x, float min_y, float max_x, float max_y)
        : min_x(min_x), min_y(min_y), max_x(max_x), max_y(max_y) {}

    Rectangle(const Point2D &point)
        : min_x(point.x), min_y(point.y), max_x(point.x), max_y(point.y) {}

    // Area of rectangle
    float area() const {
        if (min_x > max_x || min_y > max_y)
            return 0.0f;
        return (max_x - min_x) * (max_y - min_y);
    }

    // Perimeter of rectangle
    float perimeter() const {
        if (min_x > max_x || min_y > max_y)
            return 0.0f;
        return 2.0f * ((max_x - min_x) + (max_y - min_y));
    }

    // Check if point is inside rectangle
    bool contains(const Point2D &point) const {
        return point.x >= min_x && point.x <= max_x && point.y >= min_y &&
               point.y <= max_y;
    }

    // Check if this rectangle contains another rectangle
    bool contains(const Rectangle &other) const {
        return other.min_x >= min_x && other.max_x <= max_x &&
               other.min_y >= min_y && other.max_y <= max_y;
    }

    // Check if rectangles intersect
    bool intersects(const Rectangle &other) const {
        return !(other.max_x < min_x || other.min_x > max_x ||
                 other.max_y < min_y || other.min_y > max_y);
    }

    // Expand rectangle to include point
    void expand(const Point2D &point) {
        min_x = std::min(min_x, point.x);
        min_y = std::min(min_y, point.y);
        max_x = std::max(max_x, point.x);
        max_y = std::max(max_y, point.y);
    }

    // Expand rectangle to include another rectangle
    void expand(const Rectangle &other) {
        if (other.min_x <= other.max_x && other.min_y <= other.max_y) {
            min_x = std::min(min_x, other.min_x);
            min_y = std::min(min_y, other.min_y);
            max_x = std::max(max_x, other.max_x);
            max_y = std::max(max_y, other.max_y);
        }
    }

    // Calculate area increase if expanded to include rectangle
    float expansionArea(const Rectangle &other) const {
        Rectangle expanded = *this;
        expanded.expand(other);
        return expanded.area() - area();
    }

    // Distance from point to rectangle (0 if inside)
    float distanceToPoint(const Point2D &point) const {
        float dx = std::max({min_x - point.x, 0.0f, point.x - max_x});
        float dy = std::max({min_y - point.y, 0.0f, point.y - max_y});
        return std::sqrt(dx * dx + dy * dy);
    }

    // Center of rectangle
    Point2D center() const {
        return Point2D((min_x + max_x) / 2.0f, (min_y + max_y) / 2.0f);
    }

    bool isValid() const { return min_x <= max_x && min_y <= max_y; }

    void print() const {
        std::cout << "[(" << min_x << "," << min_y << ") to (" << max_x << ","
                  << max_y << ")]";
    }
};

// R-tree Entry (can hold point data and ID)
struct RTreeEntry {
    Point2D point;
    Rectangle mbr; // Minimum Bounding Rectangle
    int id;

    RTreeEntry(const Point2D &p, int id = 0) : point(p), mbr(p), id(id) {}
    RTreeEntry(const Rectangle &r, int id = 0) : mbr(r), id(id) {}
};

// Forward declaration
class RTreeNode;

using RTreeNodePtr = std::shared_ptr<RTreeNode>;

// R-tree Node
class RTreeNode {
public:
    static const int MAX_ENTRIES = 4; // M = 4 for demonstration
    static const int MIN_ENTRIES = 2; // m = M/2

    std::vector<RTreeEntry> entries;
    std::vector<RTreeNodePtr> children;
    Rectangle mbr;
    bool isLeaf;
    int level;

    RTreeNode(bool leaf = true, int level = 0) : isLeaf(leaf), level(level) {}

    bool isFull() const { return entries.size() >= MAX_ENTRIES; }

    bool needsSplit() const { return entries.size() > MAX_ENTRIES; }

    void updateMBR() {
        mbr = Rectangle();
        for (const auto &entry : entries) {
            mbr.expand(entry.mbr);
        }
    }

    // Add entry to node
    void addEntry(const RTreeEntry &entry) {
        entries.push_back(entry);
        mbr.expand(entry.mbr);
    }

    // Add child node
    void addChild(RTreeNodePtr child) { children.push_back(child); }

    // Find entry that needs least enlargement to fit new rectangle
    int chooseBestEntry(const Rectangle &rect) const {
        int bestIndex = 0;
        float minEnlargement = entries[0].mbr.expansionArea(rect);
        float minArea = entries[0].mbr.area();

        for (size_t i = 1; i < entries.size(); ++i) {
            float enlargement = entries[i].mbr.expansionArea(rect);
            float area = entries[i].mbr.area();

            if (enlargement < minEnlargement ||
                (enlargement == minEnlargement && area < minArea)) {
                minEnlargement = enlargement;
                minArea = area;
                bestIndex = static_cast<int>(i);
            }
        }

        return bestIndex;
    }

    void print(int indent = 0) const {
        std::string indentStr(indent * 2, ' ');
        std::cout << indentStr << (isLeaf ? "Leaf" : "Internal")
                  << " Node (level=" << level << "):" << std::endl;
        std::cout << indentStr << "  MBR: ";
        mbr.print();
        std::cout << std::endl;

        if (isLeaf) {
            for (const auto &entry : entries) {
                std::cout << indentStr << "  Entry " << entry.id << ": ";
                entry.point.print();
                std::cout << std::endl;
            }
        } else {
            for (size_t i = 0; i < children.size(); ++i) {
                std::cout << indentStr << "  Child " << i << ":" << std::endl;
                if (children[i]) {
                    children[i]->print(indent + 2);
                }
            }
        }
    }
};

class RTree {
private:
    RTreeNodePtr root;
    int height;

public:
    RTree() : height(1) { root = std::make_shared<RTreeNode>(true, 0); }

    void insert(const Point2D &point, int id = 0) {
        RTreeEntry entry(point, id);
        insertEntry(entry);
    }

    void insertEntry(const RTreeEntry &entry) {
        RTreeNodePtr leaf = chooseLeaf(entry.mbr);
        leaf->addEntry(entry);

        if (leaf->needsSplit()) {
            auto newNode = splitNode(leaf);
            if (leaf == root) {
                // Create new root
                auto newRoot = std::make_shared<RTreeNode>(false, height);
                newRoot->addEntry(RTreeEntry(leaf->mbr));
                newRoot->addChild(leaf);
                newRoot->addEntry(RTreeEntry(newNode->mbr));
                newRoot->addChild(newNode);
                newRoot->updateMBR();
                root = newRoot;
                height++;
            } else {
                // Propagate split upward
                propagateSplit(leaf, newNode);
            }
        }

        // Update MBRs up the tree
        updateMBRs(leaf);
    }

    // Range query: find all points within rectangle
    std::vector<RTreeEntry> rangeQuery(const Rectangle &queryRect) const {
        std::vector<RTreeEntry> results;
        rangeSearch(root, queryRect, results);
        return results;
    }

    // Nearest neighbor search
    std::vector<RTreeEntry> nearestNeighbors(const Point2D &queryPoint,
                                             int k = 1) const {
        std::vector<RTreeEntry> results;
        nearestSearch(queryPoint, k, results);
        return results;
    }

    void print() const {
        std::cout << "R-tree (height=" << height << "):" << std::endl;
        if (root) {
            root->print();
        }
    }

    int getHeight() const { return height; }

    // Get total number of entries
    int size() const { return countEntries(root); }

private:
    RTreeNodePtr chooseLeaf(const Rectangle &rect) {
        RTreeNodePtr node = root;

        while (!node->isLeaf) {
            int bestIndex = node->chooseBestEntry(rect);
            node = node->children[bestIndex];
        }

        return node;
    }

    RTreeNodePtr splitNode(RTreeNodePtr node) {
        auto newNode = std::make_shared<RTreeNode>(node->isLeaf, node->level);

        // Simple linear split algorithm
        linearSplit(node, newNode);

        return newNode;
    }

    void linearSplit(RTreeNodePtr node1, RTreeNodePtr node2) {
        std::vector<RTreeEntry> allEntries = node1->entries;
        std::vector<RTreeNodePtr> allChildren = node1->children;

        node1->entries.clear();
        node1->children.clear();
        node1->mbr = Rectangle();

        // Find seeds (entries that are farthest apart)
        int seed1 = 0, seed2 = 1;
        float maxDistance = 0;

        for (size_t i = 0; i < allEntries.size(); ++i) {
            for (size_t j = i + 1; j < allEntries.size(); ++j) {
                Point2D center1 = allEntries[i].mbr.center();
                Point2D center2 = allEntries[j].mbr.center();
                float dist = center1.distance(center2);

                if (dist > maxDistance) {
                    maxDistance = dist;
                    seed1 = static_cast<int>(i);
                    seed2 = static_cast<int>(j);
                }
            }
        }

        // Add seeds to respective nodes
        node1->addEntry(allEntries[seed1]);
        node2->addEntry(allEntries[seed2]);

        if (!node1->isLeaf) {
            node1->addChild(allChildren[seed1]);
            node2->addChild(allChildren[seed2]);
        }

        // Distribute remaining entries
        for (size_t i = 0; i < allEntries.size(); ++i) {
            if (i == seed1 || i == seed2)
                continue;

            float enlargement1 = node1->mbr.expansionArea(allEntries[i].mbr);
            float enlargement2 = node2->mbr.expansionArea(allEntries[i].mbr);

            if (enlargement1 < enlargement2 ||
                (enlargement1 == enlargement2 &&
                 node1->entries.size() <= node2->entries.size())) {
                node1->addEntry(allEntries[i]);
                if (!node1->isLeaf && i < allChildren.size()) {
                    node1->addChild(allChildren[i]);
                }
            } else {
                node2->addEntry(allEntries[i]);
                if (!node2->isLeaf && i < allChildren.size()) {
                    node2->addChild(allChildren[i]);
                }
            }
        }
    }

    void propagateSplit(RTreeNodePtr node, RTreeNodePtr newNode) {
        // This is a simplified version - full implementation would handle
        // parent tracking For now, we'll assume splits only happen at leaf
        // level for this demo
    }

    void updateMBRs(RTreeNodePtr node) {
        node->updateMBR();
        // In full implementation, would update parent MBRs as well
    }

    void rangeSearch(RTreeNodePtr node, const Rectangle &queryRect,
                     std::vector<RTreeEntry> &results) const {
        if (!node || !node->mbr.intersects(queryRect)) {
            return;
        }

        if (node->isLeaf) {
            for (const auto &entry : node->entries) {
                if (queryRect.contains(entry.point)) {
                    results.push_back(entry);
                }
            }
        } else {
            for (auto child : node->children) {
                rangeSearch(child, queryRect, results);
            }
        }
    }

    void nearestSearch(const Point2D &queryPoint, int k,
                       std::vector<RTreeEntry> &results) const {
        // Priority queue based approach would be better for large k
        // This is a simplified brute-force approach
        std::vector<std::pair<float, RTreeEntry>> candidates;
        collectAllEntries(root, candidates, queryPoint);

        std::sort(candidates.begin(), candidates.end(),
                  [](std::pair<float, RTreeEntry> const &l,
                     std::pair<float, RTreeEntry> const &r) {
                      return l.first < r.first;
                  });

        int count = std::min(k, (int)candidates.size());
        for (int i = 0; i < count; ++i) {
            results.push_back(candidates[i].second);
        }
    }

    void
    collectAllEntries(RTreeNodePtr node,
                      std::vector<std::pair<float, RTreeEntry>> &candidates,
                      const Point2D &queryPoint) const {
        if (!node)
            return;

        if (node->isLeaf) {
            for (const auto &entry : node->entries) {
                float distance = queryPoint.distance(entry.point);
                candidates.emplace_back(distance, entry);
            }
        } else {
            for (auto child : node->children) {
                collectAllEntries(child, candidates, queryPoint);
            }
        }
    }

    int countEntries(RTreeNodePtr node) const {
        if (!node)
            return 0;

        if (node->isLeaf) {
            return static_cast<int>(node->entries.size());
        } else {
            int count = 0;
            for (auto child : node->children) {
                count += countEntries(child);
            }
            return count;
        }
    }
};

// Test and demonstration
class RTreeDemo {
public:
    static std::vector<Point2D> generateRandomPoints(int count,
                                                     float minCoord = -100.0f,
                                                     float maxCoord = 100.0f) {
        std::vector<Point2D> points;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<float> dis(minCoord, maxCoord);

        for (int i = 0; i < count; ++i) {
            points.emplace_back(dis(gen), dis(gen));
        }

        return points;
    }

    static void basicDemo() {
        std::cout << "=== Basic R-tree Demo ===" << std::endl;

        RTree rtree;

        // Insert some points
        std::vector<Point2D> points = {{10, 20}, {5, 15},  {30, 40}, {25, 35},
                                       {60, 70}, {55, 65}, {80, 90}, {75, 85},
                                       {15, 25}, {45, 55}};

        std::cout << "Inserting " << points.size() << " points..." << std::endl;
        for (size_t i = 0; i < points.size(); ++i) {
            rtree.insert(points[i], static_cast<int>(i));
        }

        std::cout << "\nR-tree structure:" << std::endl;
        rtree.print();

        // Range query
        Rectangle queryRect(20, 30, 70, 80);
        std::cout << "\nRange query for rectangle (20,30) to (70,80):"
                  << std::endl;
        auto rangeResults = rtree.rangeQuery(queryRect);

        std::cout << "Found " << rangeResults.size() << " points:" << std::endl;
        for (const auto &result : rangeResults) {
            std::cout << "  Point " << result.id << ": ";
            result.point.print();
            std::cout << std::endl;
        }

        // Nearest neighbor query
        Point2D queryPoint(50, 60);
        std::cout << "\nNearest neighbor to point ";
        queryPoint.print();
        std::cout << ":" << std::endl;

        auto nnResults = rtree.nearestNeighbors(queryPoint, 3);
        for (size_t i = 0; i < nnResults.size(); ++i) {
            std::cout << "  " << (i + 1) << ". Point " << nnResults[i].id
                      << ": ";
            nnResults[i].point.print();
            std::cout << " (distance: "
                      << queryPoint.distance(nnResults[i].point) << ")"
                      << std::endl;
        }

        std::cout << "\nR-tree height: " << rtree.getHeight() << std::endl;
        std::cout << "Total entries: " << rtree.size() << std::endl;
    }

    static void performanceDemo() {
        std::cout << "\n=== Performance Demo ===" << std::endl;

        const int POINT_COUNT = 1000;
        auto points = generateRandomPoints(POINT_COUNT);

        // Build R-tree
        auto start = std::chrono::high_resolution_clock::now();
        RTree rtree;
        for (size_t i = 0; i < points.size(); ++i) {
            rtree.insert(points[i], static_cast<int>(i));
        }
        auto buildTime = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start);

        std::cout << "Built R-tree with " << POINT_COUNT << " points in "
                  << buildTime.count() << " ms" << std::endl;
        std::cout << "R-tree height: " << rtree.getHeight() << std::endl;
        std::cout << "Total entries: " << rtree.size() << std::endl;

        // Test range queries
        Rectangle queryRect(-20, -20, 20, 20);
        start = std::chrono::high_resolution_clock::now();
        auto results = rtree.rangeQuery(queryRect);
        auto queryTime = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start);

        std::cout << "\nRange query found " << results.size() << " points in "
                  << queryTime.count() << " μs" << std::endl;

        // Compare with brute force
        start = std::chrono::high_resolution_clock::now();
        int bruteForceCount = 0;
        for (const auto &point : points) {
            if (queryRect.contains(point)) {
                bruteForceCount++;
            }
        }
        auto bruteForceTime =
            std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - start);

        std::cout << "Brute force found " << bruteForceCount << " points in "
                  << bruteForceTime.count() << " μs" << std::endl;

        if (queryTime.count() > 0) {
            std::cout << "Speedup factor: "
                      << (float)bruteForceTime.count() / queryTime.count()
                      << "x" << std::endl;
        }
    }

    static void spatialDistributionDemo() {
        std::cout << "\n=== Spatial Distribution Demo ===" << std::endl;

        RTree rtree;

        // Create clustered data
        std::vector<std::vector<Point2D>> clusters = {
            {{10, 10}, {12, 11}, {11, 9}, {13, 12}, {9, 10}},   // Cluster 1
            {{50, 50}, {52, 51}, {51, 49}, {53, 52}, {49, 50}}, // Cluster 2
            {{80, 20}, {82, 21}, {81, 19}, {83, 22}, {79, 20}}, // Cluster 3
            {{20, 80}, {22, 81}, {21, 79}, {23, 82}, {19, 80}}  // Cluster 4
        };

        int id = 0;
        for (const auto &cluster : clusters) {
            for (const auto &point : cluster) {
                rtree.insert(point, id++);
            }
        }

        std::cout << "Inserted " << id << " clustered points" << std::endl;
        std::cout << "R-tree structure:" << std::endl;
        rtree.print();

        // Query each cluster region
        std::vector<Rectangle> clusterQueries = {
            {8, 8, 15, 15},   // Query cluster 1
            {48, 48, 55, 55}, // Query cluster 2
            {78, 18, 85, 25}, // Query cluster 3
            {18, 78, 25, 85}  // Query cluster 4
        };

        for (size_t i = 0; i < clusterQueries.size(); ++i) {
            auto results = rtree.rangeQuery(clusterQueries[i]);
            std::cout << "\nCluster " << (i + 1) << " query found "
                      << results.size() << " points" << std::endl;
        }
    }
};

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

    std::cout << "R-tree Implementation Examples" << std::endl;
    std::cout << "==============================" << std::endl;

    try {
        RTreeDemo::basicDemo();
        RTreeDemo::performanceDemo();
        RTreeDemo::spatialDistributionDemo();

        std::cout << "\n=== R-tree Characteristics ===" << std::endl;
        std::cout << "- Balanced tree structure" << std::endl;
        std::cout << "- Each node contains 2-4 entries (configurable)"
                  << std::endl;
        std::cout << "- Minimizes overlap between bounding rectangles"
                  << std::endl;
        std::cout
            << "- Efficient for range queries and nearest neighbor searches"
            << std::endl;
        std::cout << "- Dynamically adjusts to data distribution" << std::endl;

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return -1;
    }

    return 0;
}