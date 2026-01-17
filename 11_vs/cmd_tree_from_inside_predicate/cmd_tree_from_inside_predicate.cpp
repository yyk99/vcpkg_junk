//
// build a tree of objects from a list of the objects using an "is_inside"
// predicate
//

#include <algorithm>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

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
#endif /* CONSOLE */

// Example object class - could be rectangles, polygons, etc.
struct Rectangle {
    int id;
    std::string name;
    double x1, y1, x2, y2; // bottom-left and top-right corners

    Rectangle(int id, const std::string &name, double x1, double y1, double x2,
              double y2)
        : id(id), name(name), x1(x1), y1(y1), x2(x2), y2(y2) {}

    double area() const { return (x2 - x1) * (y2 - y1); }

    // Check if this rectangle completely contains another rectangle
    bool contains(const Rectangle &other) const {
        return x1 <= other.x1 && y1 <= other.y1 && x2 >= other.x2 &&
               y2 >= other.y2;
    }

    std::string toString() const {
        return name + "[" + std::to_string(id) + "]" + "(" +
               std::to_string(x1) + "," + std::to_string(y1) + " to " +
               std::to_string(x2) + "," + std::to_string(y2) + ")";
    }
};

// Tree node structure
template <typename T> struct TreeNode {
    T data;
    std::vector<std::unique_ptr<TreeNode<T>>> children;
    TreeNode<T> *parent;

    TreeNode(const T &item) : data(item), parent(nullptr) {}

    void addChild(std::unique_ptr<TreeNode<T>> child) {
        child->parent = this;
        children.push_back(std::move(child));
    }

    // Print tree structure
    void print(int depth = 0) const {
        std::string indent(depth * 2, ' ');
        std::cout << indent << data.toString() << " (area: " << data.area()
                  << ")\n";
        for (const auto &child : children) {
            child->print(depth + 1);
        }
    }

    // Get all nodes at this level and below
    std::vector<const TreeNode<T> *> getAllNodes() const {
        std::vector<const TreeNode<T> *> nodes;
        nodes.push_back(this);
        for (const auto &child : children) {
            auto childNodes = child->getAllNodes();
            nodes.insert(nodes.end(), childNodes.begin(), childNodes.end());
        }
        return nodes;
    }
};

// Tree builder class
template <typename T> class TreeBuilder {
public:
    using IsInsidePredicate = std::function<bool(const T &, const T &)>;

private:
    IsInsidePredicate is_inside;

    // Find the most specific parent for a given object
    TreeNode<T> *findMostSpecificParent(
        const T &item, const std::vector<std::unique_ptr<TreeNode<T>>> &roots) {
        TreeNode<T> *bestParent = nullptr;
        double smallestArea = std::numeric_limits<double>::max();

        // Search all existing nodes
        for (const auto &root : roots) {
            auto candidate = findMostSpecificParentInSubtree(item, root.get());
            if (candidate && candidate->data.area() < smallestArea) {
                bestParent = candidate;
                smallestArea = candidate->data.area();
            }
        }

        return bestParent;
    }

    TreeNode<T> *findMostSpecificParentInSubtree(const T &item,
                                                 TreeNode<T> *node) {
        TreeNode<T> *bestParent = nullptr;

        // Check if current node can contain the item
        if (is_inside(item, node->data)) {
            bestParent = node;

            // Check children for more specific parent
            for (const auto &child : node->children) {
                auto candidate =
                    findMostSpecificParentInSubtree(item, child.get());
                if (candidate) {
                    bestParent = candidate; // Child is more specific
                }
            }
        }

        return bestParent;
    }

    // Check if any existing node should become a child of the new item
    void
    reparentExistingNodes(TreeNode<T> *newNode,
                          std::vector<std::unique_ptr<TreeNode<T>>> &roots) {
        auto it = roots.begin();
        while (it != roots.end()) {
            if (is_inside((*it)->data, newNode->data)) {
                // This root should become a child of newNode
                auto child = std::move(*it);
                it = roots.erase(it);
                newNode->addChild(std::move(child));
            } else {
                // Check if any children of this root should be reparented
                reparentChildrenRecursive(newNode, it->get());
                ++it;
            }
        }
    }

    void reparentChildrenRecursive(TreeNode<T> *newNode,
                                   TreeNode<T> *currentNode) {
        auto it = currentNode->children.begin();
        while (it != currentNode->children.end()) {
            if (is_inside((*it)->data, newNode->data)) {
                // This child should move to newNode
                auto child = std::move(*it);
                it = currentNode->children.erase(it);
                newNode->addChild(std::move(child));
            } else {
                // Recursively check this child's children
                reparentChildrenRecursive(newNode, it->get());
                ++it;
            }
        }
    }

public:
    TreeBuilder(IsInsidePredicate predicate) : is_inside(predicate) {}

    std::vector<std::unique_ptr<TreeNode<T>>>
    buildTree(const std::vector<T> &objects) {
        std::vector<std::unique_ptr<TreeNode<T>>> roots;

        // Sort objects by area (largest first) for better tree construction
        std::vector<T> sortedObjects = objects;
        std::sort(sortedObjects.begin(), sortedObjects.end(),
                  [](const T &a, const T &b) { return a.area() > b.area(); });

        for (const auto &obj : sortedObjects) {
            auto newNode = std::make_unique<TreeNode<T>>(obj);

            // Find the most specific parent for this object
            TreeNode<T> *parent = findMostSpecificParent(obj, roots);

            if (parent) {
                // Check if any of parent's children should become children of
                // newNode
                reparentChildrenRecursive(newNode.get(), parent);
                // Add newNode as child of parent
                parent->addChild(std::move(newNode));
            } else {
                // Check if any existing roots should become children of newNode
                reparentExistingNodes(newNode.get(), roots);
                // Add as new root
                roots.push_back(std::move(newNode));
            }
        }

        return roots;
    }
};

// Alternative approach: Build tree with validation
template <typename T> class ValidatingTreeBuilder {
private:
    using IsInsidePredicate = std::function<bool(const T &, const T &)>;
    IsInsidePredicate is_inside;

public:
    ValidatingTreeBuilder(IsInsidePredicate predicate) : is_inside(predicate) {}

    std::vector<std::unique_ptr<TreeNode<T>>>
    buildTree(const std::vector<T> &objects) {
        // Create all nodes first
        std::vector<std::unique_ptr<TreeNode<T>>> allNodes;
        for (const auto &obj : objects) {
            allNodes.push_back(std::make_unique<TreeNode<T>>(obj));
        }

        // Build containment relationships
        for (size_t i = 0; i < allNodes.size(); ++i) {
            TreeNode<T> *currentNode = allNodes[i].get();
            TreeNode<T> *bestParent = nullptr;
            double smallestParentArea = std::numeric_limits<double>::max();

            // Find the smallest container for this node
            for (size_t j = 0; j < allNodes.size(); ++j) {
                if (i == j)
                    continue;

                TreeNode<T> *candidateParent = allNodes[j].get();
                if (is_inside(currentNode->data, candidateParent->data) &&
                    candidateParent->data.area() < smallestParentArea) {
                    bestParent = candidateParent;
                    smallestParentArea = candidateParent->data.area();
                }
            }

            if (bestParent) {
                currentNode->parent = bestParent;
            }
        }

        // Extract roots and build final tree structure
        std::vector<std::unique_ptr<TreeNode<T>>> roots;
        for (auto &node : allNodes) {
            if (node && node->parent == nullptr) {
                // This is a root node
                buildChildrenRecursive(node.get(), allNodes);
                roots.push_back(std::move(node));
            }
        }
        return roots;
    }

private:
    void buildChildrenRecursive(
        TreeNode<T> *parent,
        std::vector<std::unique_ptr<TreeNode<T>>> &allNodes) {
        for (auto &node : allNodes) {
            if (node && node->parent == parent) {
                buildChildrenRecursive(node.get(), allNodes);
                parent->children.push_back(std::move(node));
            }
        }
    }
};

// Example usage and test
TEST(main, Tree_built_using_TreeBuilder) {
    std::cout
        << "=== Building Tree from Objects using is_inside Predicate ===\n\n";

    // Create test rectangles
    std::vector<Rectangle> rectangles = {
        {1, "Outer", 0, 0, 20, 20},   // Large outer rectangle
        {2, "Middle1", 2, 2, 18, 18}, // Medium rectangle inside outer
#if 0
        {3, "Inner1", 4, 4, 8, 8},     // Small rectangle inside middle1
        {4, "Inner2", 10, 10, 16, 16}, // Another small rectangle inside middle1
        {5, "Separate", 25, 25, 35, 35}, // Separate rectangle (another root)
        {6, "InSep", 27, 27, 33, 33},    // Inside the separate rectangle
        {7, "Tiny", 5, 5, 6, 6}          // Very small, inside Inner1
#endif
    };

    std::cout << "Input rectangles:\n";
    for (const auto &rect : rectangles) {
        std::cout << "  " << rect.toString() << " (area: " << rect.area()
                  << ")\n";
    }
    std::cout << "\n";

    // Define the is_inside predicate
    auto is_inside = [](const Rectangle &inner, const Rectangle &outer) {
        return outer.contains(inner) && inner.id != outer.id;
    };

    // Build tree using first approach
    std::cout << "=== Tree built using TreeBuilder ===\n";
    TreeBuilder<Rectangle> builder(is_inside);
    auto tree1 = builder.buildTree(rectangles);

    std::cout << "Tree structure:\n";
    for (const auto &root : tree1) {
        root->print();
        std::cout << "\n";
    }

    // Demonstrate tree traversal
    std::cout << "=== Tree Traversal Example ===\n";
    if (!tree1.empty()) {
        std::cout << "All nodes in first tree:\n";
        auto allNodes = tree1[0]->getAllNodes();
        for (const auto *node : allNodes) {
            std::cout << "  " << node->data.name << " (depth: ";
            int depth = 0;
            const TreeNode<Rectangle> *current = node;
            while (current->parent) {
                depth++;
                current = current->parent;
            }
            std::cout << depth << ")\n";
        }
    }
}

TEST(main, Tree_built_using_ValidatingTreeBuilder) {
    std::cout
        << "=== Building Tree from Objects using is_inside Predicate ===\n\n";

    // Create test rectangles
    std::vector<Rectangle> rectangles = {
        {1, "Outer", 0, 0, 20, 20},   // Large outer rectangle
        {2, "Middle1", 2, 2, 18, 18}, // Medium rectangle inside outer
        {3, "Inner1", 4, 4, 8, 8},     // Small rectangle inside middle1
        {4, "Inner2", 10, 10, 16, 16}, // Another small rectangle inside middle1
        {5, "Separate", 25, 25, 35, 35}, // Separate rectangle (another root)
        {6, "InSep", 27, 27, 33, 33},    // Inside the separate rectangle
        {7, "Tiny", 5, 5, 6, 6}          // Very small, inside Inner1
    };

    std::cout << "Input rectangles:\n";
    for (const auto &rect : rectangles) {
        std::cout << "  " << rect.toString() << " (area: " << rect.area()
                  << ")\n";
    }
    std::cout << "\n";

    // Define the is_inside predicate
    auto is_inside = [](const Rectangle &inner, const Rectangle &outer) {
        return outer.contains(inner) && inner.id != outer.id;
    };

    // Build tree using validation approach
    std::cout << "=== Tree built using ValidatingTreeBuilder ===\n";
    ValidatingTreeBuilder<Rectangle> validatingBuilder(is_inside);
    auto tree2 = validatingBuilder.buildTree(rectangles);

    std::cout << "Tree structure:\n";
    for (const auto &root : tree2) {
        root->print();
        std::cout << "\n";
    }

    // Demonstrate tree traversal
    std::cout << "=== Tree Traversal Example ===\n";
    if (!tree2.empty()) {
        std::cout << "All nodes in first tree:\n";
        auto allNodes = tree2[0]->getAllNodes();
        for (const auto *node : allNodes) {
            std::cout << "  " << node->data.name << " (depth: ";
            int depth = 0;
            const TreeNode<Rectangle> *current = node;
            while (current->parent) {
                depth++;
                current = current->parent;
            }
            std::cout << depth << ")\n";
        }
    }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
