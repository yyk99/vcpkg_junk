#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <queue>
#include <sstream>
#include <string>
#include <vector>

#include <cassert>
#include <filesystem>

#ifndef M_PI
#define M_PI 3.14159265358979323846 // pi
#endif

// Point structure for 3D point cloud data
struct Point3D {
    double x, y, z;
    uint8_t r, g, b; // RGB color
    double intensity;

    Point3D(double x = 0, double y = 0, double z = 0, uint8_t r = 255,
            uint8_t g = 255, uint8_t b = 255, double intensity = 0)
        : x(x), y(y), z(z), r(r), g(g), b(b), intensity(intensity) {}
};

// Bounding box for spatial partitioning
struct BoundingBox {
    double min_x, min_y, min_z;
    double max_x, max_y, max_z;

    BoundingBox()
        : min_x(1e10), min_y(1e10), min_z(1e10), max_x(-1e10), max_y(-1e10),
          max_z(-1e10) {}

    BoundingBox(double min_x, double min_y, double min_z, double max_x,
                double max_y, double max_z)
        : min_x(min_x), min_y(min_y), min_z(min_z), max_x(max_x), max_y(max_y),
          max_z(max_z) {}

    void expand(const Point3D &point) {
        min_x = std::min(min_x, point.x);
        min_y = std::min(min_y, point.y);
        min_z = std::min(min_z, point.z);
        max_x = std::max(max_x, point.x);
        max_y = std::max(max_y, point.y);
        max_z = std::max(max_z, point.z);
    }

    bool contains(const Point3D &point) const {
        return point.x >= min_x && point.x <= max_x && point.y >= min_y &&
               point.y <= max_y && point.z >= min_z && point.z <= max_z;
    }

    double diagonal() const {
        double dx = max_x - min_x;
        double dy = max_y - min_y;
        double dz = max_z - min_z;
        return sqrt(dx * dx + dy * dy + dz * dz);
    }

    double volume() const {
        return (max_x - min_x) * (max_y - min_y) * (max_z - min_z);
    }

    Point3D center() const {
        return Point3D((min_x + max_x) / 2, (min_y + max_y) / 2,
                       (min_z + max_z) / 2);
    }
};

// Cesium 3D Tile node
struct TileNode {
    BoundingBox bounds;
    std::vector<Point3D> points;
    std::vector<std::shared_ptr<TileNode>> children;
    int level;
    std::string tile_id;
    double geometric_error;

    TileNode(int level = 0) : level(level), geometric_error(0) {}

    bool isLeaf() const { return children.empty(); }

    size_t totalPoints() const {
        size_t count = points.size();
        for (const auto &child : children) {
            count += child->totalPoints();
        }
        return count;
    }
};

// Point Cloud to 3D Tiles Converter
class PointCloudTo3DTiles {
private:
    std::vector<Point3D> original_points;
    std::shared_ptr<TileNode> root_tile;
    size_t max_points_per_tile;
    int max_depth;
    std::string output_directory;
    Point3D rtc_center;

public:
    PointCloudTo3DTiles(size_t max_points = 50000, int max_depth = 10,
                        const std::string &output_dir = "./3dtiles/")
        : max_points_per_tile(max_points), max_depth(max_depth),
          output_directory(output_dir) {}

    // Load point cloud from various formats
    bool loadFromPLY(const std::string &filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Cannot open PLY file: " << filename << std::endl;
            return false;
        }

        std::string line;
        size_t vertex_count = 0;
        bool in_header = true;

        // Parse PLY header
        while (std::getline(file, line) && in_header) {
            if (line.find("element vertex") != std::string::npos) {
                std::istringstream iss(line);
                std::string word1, word2;
                iss >> word1 >> word2 >> vertex_count;
            }
            if (line == "end_header") {
                in_header = false;
            }
        }

        // Read vertex data
        original_points.reserve(vertex_count);
        for (size_t i = 0; i < vertex_count && std::getline(file, line); ++i) {
            std::istringstream iss(line);
            double x, y, z;
            int r, g, b;

            if (iss >> x >> y >> z >> r >> g >> b) {
                original_points.emplace_back(x, y, z, static_cast<uint8_t>(r),
                                             static_cast<uint8_t>(g),
                                             static_cast<uint8_t>(b));
            }
        }

        std::cout << "Loaded " << original_points.size()
                  << " points from PLY file\n";
        return true;
    }

    bool loadFromLAS(const std::string &filename) {
        // Simplified LAS reader (real implementation would use a LAS library)
        std::cout << "LAS loading not implemented in this example\n";
        std::cout << "Consider using libraries like: libLAS, PDAL, or LASlib\n";
        return false;
    }

    bool loadFromCSV(const std::string &filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Cannot open CSV file: " << filename << std::endl;
            return false;
        }

        std::string line;
        // Skip header if present
        if (std::getline(file, line)) {
            // Check if first line is numeric (no header) or text (has header)
            std::istringstream test(line);
            double test_val;
            if (!(test >> test_val)) {
                // First line is not numeric, assume it's a header
            } else {
                // First line is numeric, process it
                std::istringstream iss(line);
                double x, y, z;
                int r = 255, g = 255, b = 255;
                char comma;

                if (iss >> x >> comma >> y >> comma >> z) {
                    iss >> comma >> r >> comma >> g >> comma >> b;
                    original_points.emplace_back(
                        x, y, z, static_cast<uint8_t>(r),
                        static_cast<uint8_t>(g), static_cast<uint8_t>(b));
                }
            }
        }

        // Read remaining lines
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            double x, y, z;
            int r = 255, g = 255, b = 255;
            char comma;

            if (iss >> x >> comma >> y >> comma >> z) {
                iss >> comma >> r >> comma >> g >> comma >> b;
                original_points.emplace_back(x, y, z, static_cast<uint8_t>(r),
                                             static_cast<uint8_t>(g),
                                             static_cast<uint8_t>(b));
            }
        }

        std::cout << "Loaded " << original_points.size()
                  << " points from CSV file\n";
        return true;
    }

    // Generate sample point cloud for testing
    void generateSamplePointCloud(size_t count = 100000) {
        original_points.clear();
        original_points.reserve(count);

        std::cout << "Generating " << count << " sample points...\n";

        for (size_t i = 0; i < count; ++i) {
            // Create a more interesting distribution
            double angle = 2 * M_PI * i / count;
            double radius = 100 + 50 * sin(8 * angle);
            double height = 10 * sin(4 * angle);

            double x = radius * cos(angle) + (rand() % 20 - 10);
            double y = radius * sin(angle) + (rand() % 20 - 10);
            double z = height + (rand() % 10 - 5);

            // Color based on height
            uint8_t r = static_cast<uint8_t>(128 + 127 * sin(height / 10));
            uint8_t g = static_cast<uint8_t>(128 + 127 * cos(height / 10));
            uint8_t b = static_cast<uint8_t>(200);

            original_points.emplace_back(x, y, z, r, g, b);
        }
    }

    // Build octree structure for tiling
    void buildTileHierarchy() {
        if (original_points.empty()) {
            std::cerr << "No points loaded!\n";
            return;
        }

        std::cout << "Building tile hierarchy...\n";

        // Calculate overall bounding box
        BoundingBox overall_bounds;
        for (const auto &point : original_points) {
            overall_bounds.expand(point);
        }

        // Create root tile
        root_tile = std::make_shared<TileNode>(0);
        root_tile->bounds = overall_bounds;
        root_tile->points = original_points;
        root_tile->tile_id = "root";
        root_tile->geometric_error = overall_bounds.diagonal();

        // Recursively subdivide
        subdivide(root_tile);

        std::cout << "Tile hierarchy built. Total tiles: "
                  << countTiles(root_tile) << "\n";
    }

    // Generate 3D Tiles output
    void generate3DTiles() {
        if (!root_tile) {
            std::cerr << "No tile hierarchy built!\n";
            return;
        }

        std::cout << "Generating 3D Tiles...\n";

        // Create output directory
#if 0
        std::string mkdir_cmd = "mkdir -p " + output_directory;
        auto rc = system(mkdir_cmd.c_str());
        assert(rc == 0);
#else
        std::error_code rc;
        std::filesystem::create_directories(output_directory, rc);
        assert(rc.value() == 0);
#endif

        // Generate tileset.json
        generateTilesetJson();

        // Generate individual tile files
        generateTileFiles(root_tile);

        std::cout << "3D Tiles generation complete!\n";
        std::cout << "Output directory: " << output_directory << "\n";
        std::cout << "Main file: " << output_directory << "tileset.json\n";
    }

    void printStatistics() const {
        if (!root_tile) {
            std::cout << "No tiles generated yet.\n";
            return;
        }

        std::cout << "\n=== 3D TILES STATISTICS ===\n";
        std::cout << "Total points: " << original_points.size() << "\n";
        std::cout << "Total tiles: " << countTiles(root_tile) << "\n";
        std::cout << "Max depth: " << getMaxDepth(root_tile) << "\n";
        std::cout << "Root bounds: (" << root_tile->bounds.min_x << ", "
                  << root_tile->bounds.min_y << ", " << root_tile->bounds.min_z
                  << ") to (" << root_tile->bounds.max_x << ", "
                  << root_tile->bounds.max_y << ", " << root_tile->bounds.max_z
                  << ")\n";
        std::cout << "Root geometric error: " << root_tile->geometric_error
                  << "\n";
    }

private:
    void subdivide(std::shared_ptr<TileNode> node) {
        if (node->level >= max_depth ||
            node->points.size() <= max_points_per_tile) {
            return;
        }

        // Create 8 child nodes (octree subdivision)
        Point3D center = node->bounds.center();

        for (int i = 0; i < 8; ++i) {
            auto child = std::make_shared<TileNode>(node->level + 1);

            // Determine child bounds
            double min_x = (i & 1) ? center.x : node->bounds.min_x;
            double max_x = (i & 1) ? node->bounds.max_x : center.x;
            double min_y = (i & 2) ? center.y : node->bounds.min_y;
            double max_y = (i & 2) ? node->bounds.max_y : center.y;
            double min_z = (i & 4) ? center.z : node->bounds.min_z;
            double max_z = (i & 4) ? node->bounds.max_z : center.z;

            child->bounds =
                BoundingBox(min_x, min_y, min_z, max_x, max_y, max_z);
            child->tile_id = node->tile_id + "_" + std::to_string(i);
            child->geometric_error = node->geometric_error / 2.0;

            // Distribute points to children
            for (const auto &point : node->points) {
                if (child->bounds.contains(point)) {
                    child->points.push_back(point);
                }
            }

            if (!child->points.empty()) {
                node->children.push_back(child);
                subdivide(child);
            }
        }

        // Clear points from internal nodes (keep only in leaves for this
        // example)
        if (!node->children.empty()) {
            node->points.clear();
        }
    }

    void generateTilesetJson() const {
        std::ofstream file(output_directory + "tileset.json");

        assert(file.good());

        file << "{\n";
        file << "  \"asset\": {\n";
        file << "    \"version\": \"1.0\",\n";
        file << "    \"generator\": \"Custom Point Cloud to 3D Tiles "
                "Converter\"\n";
        file << "  },\n";
        file << "  \"geometricError\": " << root_tile->geometric_error << ",\n";
        file << "  \"root\": {\n";

        writeTileJson(file, root_tile, 2);

        file << "  }\n";
        file << "}\n";

        file.close();
    }

    void writeTileJson(std::ofstream &file, std::shared_ptr<TileNode> tile,
                       int indent) const {
        std::string ind(indent, ' ');

        file << ind << "\"boundingVolume\": {\n";
        file << ind << "  \"box\": [\n";
        file << ind << "    " << (tile->bounds.min_x + tile->bounds.max_x) / 2
             << ",\n";
        file << ind << "    " << (tile->bounds.min_y + tile->bounds.max_y) / 2
             << ",\n";
        file << ind << "    " << (tile->bounds.min_z + tile->bounds.max_z) / 2
             << ",\n";
        file << ind << "    " << (tile->bounds.max_x - tile->bounds.min_x) / 2
             << ", 0, 0,\n";
        file << ind << "    0, "
             << (tile->bounds.max_y - tile->bounds.min_y) / 2 << ", 0,\n";
        file << ind << "    0, 0, "
             << (tile->bounds.max_z - tile->bounds.min_z) / 2 << "\n";
        file << ind << "  ]\n";
        file << ind << "},\n";
        file << ind << "\"geometricError\": " << tile->geometric_error << ",\n";

        if (tile->isLeaf() && !tile->points.empty()) {
            file << ind << "\"content\": {\n";
            file << ind << "  \"uri\": \"" << tile->tile_id << ".pnts\"\n";
            file << ind << "}";
        } else {
            file << ind << "\"refine\": \"REPLACE\"";
        }

        if (!tile->children.empty()) {
            file << ",\n" << ind << "\"children\": [\n";
            for (size_t i = 0; i < tile->children.size(); ++i) {
                file << ind << "  {\n";
                writeTileJson(file, tile->children[i], indent + 4);
                file << ind << "  }";
                if (i < tile->children.size() - 1)
                    file << ",";
                file << "\n";
            }
            file << ind << "]\n";
        } else {
            file << "\n";
        }
    }

    void generateTileFiles(std::shared_ptr<TileNode> tile) const {
        if (tile->isLeaf() && !tile->points.empty()) {
            generatePntsFile(tile);
        }

        for (const auto &child : tile->children) {
            generateTileFiles(child);
        }
    }

    void generatePntsFile(std::shared_ptr<TileNode> tile) const {
        std::string filename = output_directory + tile->tile_id + ".pnts";
        std::ofstream file(filename, std::ios::binary);

        if (!file.is_open()) {
            std::cerr << "Cannot create PNTS file: " << filename << std::endl;
            return;
        }

        // PNTS Header
        // struct header_t
        //{
        //    char magic[4];
        //    uint32_t version;
        //    uint32_t byteLength;
        //    uint32_t featureTableJSONByteLength;
        //    uint32_t featureTableBinaryByteLength;
        //    uint32_t batchTableJSONByteLength;
        //    uint32_t batchTableBinaryByteLength;
        //};

        uint32_t version = 1;
        uint32_t point_count = static_cast<uint32_t>(tile->points.size());
        uint32_t positions_size = point_count * 12; // 3 floats per point
        uint32_t colors_size = point_count * 3;     // 3 bytes per point
        uint32_t feature_table_json_size = 256;     // Fixed size for simplicity
        uint32_t feature_table_binary_size = positions_size + colors_size;
        uint32_t total_size =
            28 + feature_table_json_size + feature_table_binary_size;

        // Write PNTS header
        file.write("pnts", 4);                                   // Magic: 0
        file.write(reinterpret_cast<const char *>(&version), 4); // Version: 4
        file.write(reinterpret_cast<const char *>(&total_size),
                   4); // byteLength: 8
        file.write(reinterpret_cast<const char *>(&feature_table_json_size),
                   4); // featureTableJSONByteLength:12
        file.write(reinterpret_cast<const char *>(&feature_table_binary_size),
                   4); // featureTableBinaryByteLength:16
        uint32_t batch_table_json_size = 0;
        uint32_t batch_table_binary_size = 0;
        file.write(reinterpret_cast<const char *>(&batch_table_json_size),
                   4); // batchTableJSONByteLength:20
        file.write(reinterpret_cast<const char *>(&batch_table_binary_size),
                   4); // batchTableBinaryByteLength:24

        // Feature table JSON
        std::string feature_json =
            "{"
            "\"POINTS_LENGTH\":" +
            std::to_string(point_count) +
            ","
            "\"POSITION\":{\"byteOffset\":0},"
            "\"RGB\":{\"byteOffset\":" +
            std::to_string(positions_size) + "}" + "," + "\"RTC_CENTER\":[" +
            std::to_string(rtc_center.x) + "," + std::to_string(rtc_center.y) +
            "," + std::to_string(rtc_center.y) +
            "]"
            "}";

        feature_json.resize(feature_table_json_size, ' ');
        file.write(feature_json.c_str(), feature_table_json_size);

        // Feature table binary - positions
        for (const auto &point : tile->points) {
            float pos[3] = {static_cast<float>(point.x),
                            static_cast<float>(point.y),
                            static_cast<float>(point.z)};
            file.write(reinterpret_cast<const char *>(pos), 12);
        }

        // Feature table binary - colors
        for (const auto &point : tile->points) {
            uint8_t color[3] = {point.r, point.g, point.b};
            file.write(reinterpret_cast<const char *>(color), 3);
        }

        file.close();
    }

    size_t countTiles(std::shared_ptr<TileNode> tile) const {
        size_t count = 1;
        for (const auto &child : tile->children) {
            count += countTiles(child);
        }
        return count;
    }

    int getMaxDepth(std::shared_ptr<TileNode> tile) const {
        int max_depth = tile->level;
        for (const auto &child : tile->children) {
            max_depth = std::max(max_depth, getMaxDepth(child));
        }
        return max_depth;
    }
};

// Example usage and demonstration
int main() {
    std::cout << "=== CESIUM 3D TILES POINT CLOUD CONVERTER ===\n\n";

    // Create converter
    PointCloudTo3DTiles converter(10000, 8, "./output_3dtiles/");

    std::cout << "1. LOADING POINT CLOUD DATA\n";

    // Option 1: Generate sample data
    converter.generateSamplePointCloud(50000);

    // Option 2: Load from CSV (uncomment to use)
    // if (!converter.loadFromCSV("pointcloud.csv")) {
    //     std::cerr << "Failed to load CSV file\n";
    //     return 1;
    // }

    // Option 3: Load from PLY (uncomment to use)
    // if (!converter.loadFromPLY("pointcloud.ply")) {
    //     std::cerr << "Failed to load PLY file\n";
    //     return 1;
    // }

    std::cout << "\n2. BUILDING TILE HIERARCHY\n";
    converter.buildTileHierarchy();

    std::cout << "\n3. GENERATING 3D TILES\n";
    converter.generate3DTiles();

    std::cout << "\n4. STATISTICS\n";
    converter.printStatistics();

    std::cout << "\n=== USAGE INSTRUCTIONS ===\n";
    std::cout << "1. Copy the generated 3D tiles to a web server\n";
    std::cout << "2. Use CesiumJS to load the tileset:\n\n";

    std::cout << "JavaScript example:\n";
    std::cout << "const tileset = new Cesium.Cesium3DTileset({\n";
    std::cout << "  url: 'path/to/your/tileset.json'\n";
    std::cout << "});\n";
    std::cout << "viewer.scene.primitives.add(tileset);\n\n";

    std::cout << "HTML example:\n";
    std::cout << "<!DOCTYPE html>\n";
    std::cout << "<html>\n";
    std::cout << "<head>\n";
    std::cout << "  <script "
                 "src=\"https://cesium.com/downloads/cesiumjs/releases/1.95/"
                 "Build/Cesium/Cesium.js\"></script>\n";
    std::cout << "  <link "
                 "href=\"https://cesium.com/downloads/cesiumjs/releases/1.95/"
                 "Build/Cesium/Widgets/widgets.css\" rel=\"stylesheet\">\n";
    std::cout << "</head>\n";
    std::cout << "<body>\n";
    std::cout << "  <div id=\"cesiumContainer\" style=\"width: 100%; height: "
                 "500px;\"></div>\n";
    std::cout << "  <script>\n";
    std::cout << "    const viewer = new Cesium.Viewer('cesiumContainer');\n";
    std::cout << "    const tileset = new Cesium.Cesium3DTileset({\n";
    std::cout << "      url: 'tileset.json'\n";
    std::cout << "    });\n";
    std::cout << "    viewer.scene.primitives.add(tileset);\n";
    std::cout << "    viewer.zoomTo(tileset);\n";
    std::cout << "  </script>\n";
    std::cout << "</body>\n";
    std::cout << "</html>\n\n";

    std::cout << "=== CONVERSION COMPLETE ===\n";

    return 0;
}
