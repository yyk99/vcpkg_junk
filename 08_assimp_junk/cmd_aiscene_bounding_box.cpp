#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <glm/glm.hpp>
//#include <glm/gtc/matrix_transform.hpp>
#include <limits>
//#include <algorithm>


#include <assimp/version.h>

// Method 1: Using version macros (compile-time)
void printCompileTimeVersion() {
    printf("=== COMPILE-TIME VERSION INFO ===\n");
    printf("Assimp Major Version: %d\n", aiGetVersionMajor());
    printf("Assimp Minor Version: %d\n", aiGetVersionMinor());
    printf("Assimp Patch Version: %d\n", aiGetVersionPatch());
    printf("Assimp Revision: %d\n", aiGetVersionRevision());
    printf("Full Version: %d.%d.%d.%d\n", aiGetVersionMajor(),
           aiGetVersionMinor(), aiGetVersionPatch(), aiGetVersionRevision());
    printf("\n");
}


struct BoundingBox {
    glm::vec3 min{std::numeric_limits<float>::max()};
    glm::vec3 max{std::numeric_limits<float>::lowest()};
    
    BoundingBox() = default;
    
    void expand(const glm::vec3& point) {
        min = glm::min(min, point);
        max = glm::max(max, point);
    }
    
    void merge(const BoundingBox& other) {
        if (other.min.x != std::numeric_limits<float>::max()) {
            expand(other.min);
            expand(other.max);
        }
    }
    
    bool isValid() const {
        return min.x != std::numeric_limits<float>::max();
    }
    
    glm::vec3 getCenter() const {
        return (min + max) * 0.5f;
    }
    
    glm::vec3 getSize() const {
        return max - min;
    }
};

// Convert aiMatrix4x4 to glm::mat4
glm::mat4 aiMatrixToGlm(const aiMatrix4x4& aiMat) {
    return glm::mat4(
        aiMat.a1, aiMat.b1, aiMat.c1, aiMat.d1,
        aiMat.a2, aiMat.b2, aiMat.c2, aiMat.d2,
        aiMat.a3, aiMat.b3, aiMat.c3, aiMat.d3,
        aiMat.a4, aiMat.b4, aiMat.c4, aiMat.d4
    );
}

// Compute bounding box for a single mesh with given transform
BoundingBox computeMeshBoundingBox(const aiMesh* mesh, const glm::mat4& transform) {
    BoundingBox bbox;
    
    for (unsigned int i = 0; i < mesh->mNumVertices; ++i) {
        // Transform vertex to world space
        glm::vec4 vertex(mesh->mVertices[i].x, mesh->mVertices[i].y, mesh->mVertices[i].z, 1.0f);
        glm::vec4 transformedVertex = transform * vertex;
        
        bbox.expand(glm::vec3(transformedVertex));
    }
    
    return bbox;
}

// Recursively traverse scene nodes and compute bounding box
void computeNodeBoundingBox(const aiScene* scene, const aiNode* node, 
                           const glm::mat4& parentTransform, BoundingBox& sceneBBox) {
    
    // Combine parent transform with current node's transform
    glm::mat4 nodeTransform = aiMatrixToGlm(node->mTransformation);
    glm::mat4 currentTransform = parentTransform * nodeTransform;
    
    // Process all meshes attached to this node
    for (unsigned int i = 0; i < node->mNumMeshes; ++i) {
        unsigned int meshIndex = node->mMeshes[i];
        const aiMesh* mesh = scene->mMeshes[meshIndex];
        
        BoundingBox meshBBox = computeMeshBoundingBox(mesh, currentTransform);
        sceneBBox.merge(meshBBox);
    }
    
    // Recursively process child nodes
    for (unsigned int i = 0; i < node->mNumChildren; ++i) {
        computeNodeBoundingBox(scene, node->mChildren[i], currentTransform, sceneBBox);
    }
}

// Main function to compute scene bounding box
BoundingBox computeSceneBoundingBox(const aiScene* scene) {
    if (!scene || !scene->mRootNode) {
        return BoundingBox(); // Return invalid bounding box
    }
    
    BoundingBox sceneBBox;
    glm::mat4 identity(1.0f); // Identity matrix as initial transform
    
    computeNodeBoundingBox(scene, scene->mRootNode, identity, sceneBBox);
    
    return sceneBBox;
}

// Alternative: Compute bounding box with custom root transform
BoundingBox computeSceneBoundingBox(const aiScene* scene, const glm::mat4& rootTransform) {
    if (!scene || !scene->mRootNode) {
        return BoundingBox();
    }
    
    BoundingBox sceneBBox;
    computeNodeBoundingBox(scene, scene->mRootNode, rootTransform, sceneBBox);
    
    return sceneBBox;
}

// Utility function to print bounding box info
void printBoundingBox(const BoundingBox& bbox) {
    if (!bbox.isValid()) {
        printf("Invalid bounding box\n");
        return;
    }
    
    glm::vec3 center = bbox.getCenter();
    glm::vec3 size = bbox.getSize();
    
    printf("Bounding Box:\n");
    printf("  Min: (%.3f, %.3f, %.3f)\n", bbox.min.x, bbox.min.y, bbox.min.z);
    printf("  Max: (%.3f, %.3f, %.3f)\n", bbox.max.x, bbox.max.y, bbox.max.z);
    printf("  Center: (%.3f, %.3f, %.3f)\n", center.x, center.y, center.z);
    printf("  Size: (%.3f, %.3f, %.3f)\n", size.x, size.y, size.z);
}

// Example usage

#include <filesystem>
#include <assimp/Importer.hpp>

int main() {

    printCompileTimeVersion();

    namespace fs = std::filesystem;
    auto test_data = [](const char *relative_path) {
        auto test_data_dir = fs::absolute(__FILE__).parent_path() / "test_data";
        return test_data_dir / relative_path;
    };

    auto bounding_boxes_glb = test_data("tileset2/bounding_boxes.glb").string();

    bounding_boxes_glb = R"(C:\home\work\GM-19017\out\model.glb)";

    // Assuming you have loaded an aiScene
    Assimp::Importer importer;
    const aiScene *scene = importer.ReadFile(bounding_boxes_glb, 0);

    if (!scene) {
        printf("Cannot open %s. Error: %s\n", bounding_boxes_glb.c_str(),
               importer.GetErrorString());
        return 1;
    }
    
    // Compute bounding box
    BoundingBox bbox = computeSceneBoundingBox(scene);
    
    // Print results
    printBoundingBox(bbox);
    
    return 0;
}
