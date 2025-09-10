//
//
//

// Define implementation macros once per project.
#define TINYGLTF_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <tiny_gltf.h>


#include <iostream>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>

#ifndef M_PI
/** PI definition */
#define M_PI 3.14159265358979323846
/* 3.1415926535897932384626433832795 */
#endif

// 3D Vector helper
struct Vec3 {
    float x, y, z;
    
    Vec3(float x = 0, float y = 0, float z = 0) : x(x), y(y), z(z) {}
    
    Vec3 operator+(const Vec3& other) const {
        return Vec3(x + other.x, y + other.y, z + other.z);
    }
    
    Vec3 operator-(const Vec3& other) const {
        return Vec3(x - other.x, y - other.y, z - other.z);
    }
    
    Vec3 operator*(float scalar) const {
        return Vec3(x * scalar, y * scalar, z * scalar);
    }
    
    float length() const {
        return std::sqrt(x * x + y * y + z * z);
    }
    
    float dot(const Vec3& other) const {
        return x * other.x + y * other.y + z * other.z;
    }
};

// Axis-Aligned Bounding Box
struct AABB {
    Vec3 min, max;
    
    AABB() : min(Vec3(std::numeric_limits<float>::max())),
             max(Vec3(std::numeric_limits<float>::lowest())) {}
    
    AABB(const Vec3& min, const Vec3& max) : min(min), max(max) {}
    
    void expand(const Vec3& point) {
        min.x = std::min(min.x, point.x);
        min.y = std::min(min.y, point.y);
        min.z = std::min(min.z, point.z);
        max.x = std::max(max.x, point.x);
        max.y = std::max(max.y, point.y);
        max.z = std::max(max.z, point.z);
    }
    
    void expand(const AABB& other) {
        expand(other.min);
        expand(other.max);
    }
    
    Vec3 center() const {
        return Vec3((min.x + max.x) * 0.5f, (min.y + max.y) * 0.5f, (min.z + max.z) * 0.5f);
    }
    
    Vec3 size() const {
        return max - min;
    }
    
    float diagonal() const {
        return size().length();
    }
    
    float volume() const {
        Vec3 s = size();
        return s.x * s.y * s.z;
    }
    
    bool isValid() const {
        return min.x <= max.x && min.y <= max.y && min.z <= max.z;
    }
    
    void print() const {
        std::cout << "AABB: min(" << min.x << ", " << min.y << ", " << min.z 
                  << ") max(" << max.x << ", " << max.y << ", " << max.z << ")\n";
        std::cout << "  Center: (" << center().x << ", " << center().y << ", " << center().z << ")\n";
        std::cout << "  Size: (" << size().x << ", " << size().y << ", " << size().z << ")\n";
        std::cout << "  Diagonal: " << diagonal() << "\n";
        std::cout << "  Volume: " << volume() << "\n";
    }
};

// Bounding Sphere
struct BoundingSphere {
    Vec3 center;
    float radius;
    
    BoundingSphere() : center(), radius(0) {}
    BoundingSphere(const Vec3& center, float radius) : center(center), radius(radius) {}
    
    void print() const {
        std::cout << "Bounding Sphere: center(" << center.x << ", " << center.y 
                  << ", " << center.z << ") radius(" << radius << ")\n";
        std::cout << "  Volume: " << (4.0f / 3.0f * M_PI * radius * radius * radius) << "\n";
    }
};

// Oriented Bounding Box (simplified)
struct OBB {
    Vec3 center;
    Vec3 extents;  // Half-sizes along each axis
    // Note: Full OBB would include orientation matrix
    
    OBB() : center(), extents() {}
    OBB(const Vec3& center, const Vec3& extents) : center(center), extents(extents) {}
    
    void print() const {
        std::cout << "OBB: center(" << center.x << ", " << center.y << ", " << center.z 
                  << ") extents(" << extents.x << ", " << extents.y << ", " << extents.z << ")\n";
    }
};

class GLTFBoundingVolumeCalculator {
private:
    const tinygltf::Model& model;
    
public:
    explicit GLTFBoundingVolumeCalculator(const tinygltf::Model& model) : model(model) {}
    
    // Get vertices from a primitive
    std::vector<Vec3> extractVertices(const tinygltf::Primitive& primitive) const {
        std::vector<Vec3> vertices;
        
        auto pos_it = primitive.attributes.find("POSITION");
        if (pos_it == primitive.attributes.end()) {
            return vertices;
        }
        
        int pos_accessor_idx = pos_it->second;
        const tinygltf::Accessor& pos_accessor = model.accessors[pos_accessor_idx];
        const tinygltf::BufferView& pos_buffer_view = model.bufferViews[pos_accessor.bufferView];
        const tinygltf::Buffer& pos_buffer = model.buffers[pos_buffer_view.buffer];
        
        // Get pointer to position data
        const float* positions = reinterpret_cast<const float*>(
            &pos_buffer.data[pos_buffer_view.byteOffset + pos_accessor.byteOffset]);
        
        size_t stride = pos_buffer_view.byteStride;
        if (stride == 0) {
            stride = sizeof(float) * 3; // Default stride for Vec3
        }
        
        vertices.reserve(pos_accessor.count);
        
        for (size_t i = 0; i < pos_accessor.count; ++i) {
            const float* vertex_data = reinterpret_cast<const float*>(
                reinterpret_cast<const uint8_t*>(positions) + i * stride);
            vertices.emplace_back(vertex_data[0], vertex_data[1], vertex_data[2]);
        }
        
        return vertices;
    }
    
    // Calculate AABB for a single primitive
    AABB calculatePrimitiveAABB(const tinygltf::Primitive& primitive) const {
        AABB aabb;
        
        // Check if accessor already has min/max (common optimization in GLTF)
        auto pos_it = primitive.attributes.find("POSITION");
        if (pos_it != primitive.attributes.end()) {
            const tinygltf::Accessor& pos_accessor = model.accessors[pos_it->second];
            
            if (pos_accessor.minValues.size() == 3 && pos_accessor.maxValues.size() == 3) {
                aabb.min = Vec3(static_cast<float>(pos_accessor.minValues[0]),
                               static_cast<float>(pos_accessor.minValues[1]),
                               static_cast<float>(pos_accessor.minValues[2]));
                aabb.max = Vec3(static_cast<float>(pos_accessor.maxValues[0]),
                               static_cast<float>(pos_accessor.maxValues[1]),
                               static_cast<float>(pos_accessor.maxValues[2]));
                return aabb;
            }
        }
        
        // Fall back to computing from vertices
        std::vector<Vec3> vertices = extractVertices(primitive);
        for (const Vec3& vertex : vertices) {
            aabb.expand(vertex);
        }
        
        return aabb;
    }
    
    // Calculate AABB for a mesh
    AABB calculateMeshAABB(const tinygltf::Mesh& mesh) const {
        AABB mesh_aabb;
        
        for (const auto& primitive : mesh.primitives) {
            AABB primitive_aabb = calculatePrimitiveAABB(primitive);
            if (primitive_aabb.isValid()) {
                mesh_aabb.expand(primitive_aabb);
            }
        }
        
        return mesh_aabb;
    }
    
    // Calculate AABB for entire model
    AABB calculateModelAABB() const {
        AABB model_aabb;
        
        for (const auto& mesh : model.meshes) {
            AABB mesh_aabb = calculateMeshAABB(mesh);
            if (mesh_aabb.isValid()) {
                model_aabb.expand(mesh_aabb);
            }
        }
        
        return model_aabb;
    }
    
    // Calculate bounding sphere from AABB (simple method)
    BoundingSphere calculateBoundingSphereFromAABB(const AABB& aabb) const {
        Vec3 center = aabb.center();
        float radius = aabb.diagonal() * 0.5f;
        return BoundingSphere(center, radius);
    }
    
    // Calculate optimal bounding sphere using Ritter's algorithm
    BoundingSphere calculateOptimalBoundingSphere(const std::vector<Vec3>& points) const {
        if (points.empty()) {
            return BoundingSphere();
        }
        
        if (points.size() == 1) {
            return BoundingSphere(points[0], 0);
        }
        
        // Find two points that are farthest apart
        Vec3 p1 = points[0];
        Vec3 p2 = points[0];
        float max_dist_sq = 0;
        
        for (size_t i = 0; i < points.size(); ++i) {
            for (size_t j = i + 1; j < points.size(); ++j) {
                Vec3 diff = points[i] - points[j];
                float dist_sq = diff.dot(diff);
                if (dist_sq > max_dist_sq) {
                    max_dist_sq = dist_sq;
                    p1 = points[i];
                    p2 = points[j];
                }
            }
        }
        
        // Initial sphere
        Vec3 center = (p1 + p2) * 0.5f;
        float radius = (p2 - p1).length() * 0.5f;
        
        // Expand sphere to include all points
        for (const Vec3& point : points) {
            Vec3 diff = point - center;
            float dist = diff.length();
            
            if (dist > radius) {
                float new_radius = (radius + dist) * 0.5f;
                float factor = new_radius - radius;
                center = center + diff * (factor / dist);
                radius = new_radius;
            }
        }
        
        return BoundingSphere(center, radius);
    }
    
    // Calculate model bounding sphere
    BoundingSphere calculateModelBoundingSphere() const {
        std::vector<Vec3> all_vertices;
        
        // Collect all vertices from all meshes
        for (const auto& mesh : model.meshes) {
            for (const auto& primitive : mesh.primitives) {
                std::vector<Vec3> vertices = extractVertices(primitive);
                all_vertices.insert(all_vertices.end(), vertices.begin(), vertices.end());
            }
        }
        
        return calculateOptimalBoundingSphere(all_vertices);
    }
    
    // Calculate OBB (simplified - uses AABB as approximation)
    OBB calculateModelOBB() const {
        AABB aabb = calculateModelAABB();
        Vec3 center = aabb.center();
        Vec3 extents = aabb.size() * 0.5f;
        return OBB(center, extents);
    }
    
    // Print all bounding volume statistics
    void printAllBoundingVolumes() const {
        std::cout << "\n=== MODEL BOUNDING VOLUMES ===\n";
        
        // AABB
        AABB model_aabb = calculateModelAABB();
        std::cout << "\n--- Axis-Aligned Bounding Box ---\n";
        model_aabb.print();
        
        // Bounding Sphere (from AABB)
        BoundingSphere sphere_from_aabb = calculateBoundingSphereFromAABB(model_aabb);
        std::cout << "\n--- Bounding Sphere (from AABB) ---\n";
        sphere_from_aabb.print();
        
        // Optimal Bounding Sphere
        BoundingSphere optimal_sphere = calculateModelBoundingSphere();
        std::cout << "\n--- Optimal Bounding Sphere ---\n";
        optimal_sphere.print();
        
        // OBB (simplified)
        OBB model_obb = calculateModelOBB();
        std::cout << "\n--- Oriented Bounding Box (simplified) ---\n";
        model_obb.print();
        
        std::cout << "\n=== PER-MESH ANALYSIS ===\n";
        for (size_t i = 0; i < model.meshes.size(); ++i) {
            const auto& mesh = model.meshes[i];
            AABB mesh_aabb = calculateMeshAABB(mesh);
            
            std::cout << "\nMesh " << i;
            if (!mesh.name.empty()) {
                std::cout << " ('" << mesh.name << "')";
            }
            std::cout << ":\n";
            mesh_aabb.print();
        }
    }
};

// Example usage
int main() {

    namespace fs = std::filesystem;
    auto test_data = [](const char *relative_path) {
        auto test_data_dir = fs::absolute(__FILE__).parent_path() / "test_data";
        return test_data_dir / relative_path;
    };

    auto filename = test_data("tileset2/bounding_boxes.glb").string();

    filename = R"(C:\home\work\GM-19017\out\model.glb)";

    std::cout << "=== TINYGLTF BOUNDING VOLUME CALCULATOR ===\n";
    
    // Load GLTF model
    tinygltf::TinyGLTF loader;
    tinygltf::Model model;
    std::string err, warn;
    
    // std::string filename = "model.gltf"; // Replace with your GLTF file
    
    bool success = loader.LoadBinaryFromFile(&model, &err, &warn, filename);
#if 0
    if (!success) {
        std::cerr << "Failed to load GLTF: " << err << std::endl;
        
        // Create a simple test case if file doesn't exist
        std::cout << "\nCreating test data since file couldn't be loaded...\n";
        
        // This is just for demonstration - you'd normally load a real file
        model.meshes.resize(1);
        model.accessors.resize(1);
        model.bufferViews.resize(1);
        model.buffers.resize(1);
        
        // Add some test vertex data
        std::vector<float> test_positions = {
            -1.0f, -1.0f, -1.0f,  // Vertex 0
             1.0f, -1.0f, -1.0f,  // Vertex 1
             1.0f,  1.0f, -1.0f,  // Vertex 2
            -1.0f,  1.0f,  1.0f   // Vertex 3
        };
        
        model.buffers[0].data.resize(test_positions.size() * sizeof(float));
        std::memcpy(model.buffers[0].data.data(), test_positions.data(), 
                   model.buffers[0].data.size());
        
        model.bufferViews[0].buffer = 0;
        model.bufferViews[0].byteLength = model.buffers[0].data.size();
        model.bufferViews[0].byteOffset = 0;
        
        model.accessors[0].bufferView = 0;
        model.accessors[0].count = 4;
        model.accessors[0].type = TINYGLTF_TYPE_VEC3;
        model.accessors[0].componentType = TINYGLTF_COMPONENT_TYPE_FLOAT;
        model.accessors[0].minValues = {-1.0, -1.0, -1.0};
        model.accessors[0].maxValues = {1.0, 1.0, 1.0};
        
        model.meshes[0].primitives.resize(1);
        model.meshes[0].primitives[0].attributes["POSITION"] = 0;
    }
#else
    if (!success) {
        std::cout << "Cannot load " << filename << ". Error: " << err << "\n";
        return 1;
    }
#endif
    if (!warn.empty()) {
        std::cout << "Warning: " << warn << std::endl;
    }
    
    std::cout << "Model loaded successfully!\n";
    std::cout << "Meshes: " << model.meshes.size() << "\n";
    std::cout << "Nodes: " << model.nodes.size() << "\n";
    std::cout << "Accessors: " << model.accessors.size() << "\n";
    std::cout << std::setprecision(15);
    
    // Calculate bounding volumes
    GLTFBoundingVolumeCalculator calculator(model);
    calculator.printAllBoundingVolumes();
    
    // Example: Get just the model AABB
    std::cout << "\n=== QUICK USAGE EXAMPLES ===\n";
    
    AABB model_bounds = calculator.calculateModelAABB();
    std::cout << "Model center: (" << model_bounds.center().x << ", " 
              << model_bounds.center().y << ", " << model_bounds.center().z << ")\n";
    std::cout << "Model size: " << model_bounds.diagonal() << "\n";
    
    BoundingSphere sphere = calculator.calculateModelBoundingSphere();
    std::cout << "Bounding sphere radius: " << sphere.radius << "\n";
    
    return 0;
}