/*

This comprehensive example shows how to decode Draco-compressed Cesium PNTS
(point cloud tile) files. Here are the key components:

Key Features:

1. PNTS File Format Handling

Parses the Cesium PNTS header (magic number, version, byte lengths)
Extracts Feature Table JSON and binary data
Automatically detects Draco compression

2. Draco Decompression

Uses the Draco decoder library to decompress point cloud data
Extracts standard attributes: POSITION, COLOR, NORMAL
Handles generic attributes (like BATCH_ID)
Supports both RGB and RGBA color formats

3. Attribute Extraction

Positions: 3D coordinates (x, y, z)
Colors: RGB or RGBA per point
Normals: Optional surface normals
Batch IDs: For feature identification

*/

#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdint.h>
#include <vector>

// Include Draco decoder headers
#include "draco/compression/decode.h"
#include "draco/point_cloud/point_cloud.h"

// Cesium PNTS file format structures
#pragma pack(push, 1)
struct PntsHeader {
    char magic[4];       // "pnts"
    uint32_t version;    // Version number (usually 1)
    uint32_t byteLength; // Total byte length of the file
    uint32_t featureTableJSONByteLength;
    uint32_t featureTableBinaryByteLength;
    uint32_t batchTableJSONByteLength;
    uint32_t batchTableBinaryByteLength;
};
#pragma pack(pop)

struct PointCloudData {
    std::vector<float> positions;   // x, y, z coordinates
    std::vector<uint8_t> colors;    // RGB or RGBA colors
    std::vector<float> normals;     // Optional normals
    std::vector<uint32_t> batchIds; // Optional batch IDs
    uint32_t pointCount;

    PointCloudData() : pointCount(0) {}
};

class CesiumPntsDecoder {
public:
    // Decode a PNTS file from memory buffer
    static bool decodePnts(const std::vector<uint8_t> &fileData,
                           PointCloudData &output) {
        if (fileData.size() < sizeof(PntsHeader)) {
            std::cerr << "File too small to contain PNTS header\n";
            return false;
        }

        // Read header
        PntsHeader header;
        std::memcpy(&header, fileData.data(), sizeof(PntsHeader));

        // Validate magic number
        if (std::strncmp(header.magic, "pnts", 4) != 0) {
            std::cerr << "Invalid PNTS magic number\n";
            return false;
        }

        std::cout << "PNTS Header:\n";
        std::cout << "  Version: " << header.version << "\n";
        std::cout << "  Byte Length: " << header.byteLength << "\n";
        std::cout << "  Feature Table JSON Length: "
                  << header.featureTableJSONByteLength << "\n";
        std::cout << "  Feature Table Binary Length: "
                  << header.featureTableBinaryByteLength << "\n";

        // Calculate offsets
        size_t offset = sizeof(PntsHeader);
        size_t jsonStart = offset;
        size_t jsonEnd = jsonStart + header.featureTableJSONByteLength;
        size_t binaryStart = jsonEnd;
        size_t binaryEnd = binaryStart + header.featureTableBinaryByteLength;

        // Extract feature table JSON (optional - parse to get metadata)
        if (header.featureTableJSONByteLength > 0) {
            std::string featureTableJSON(
                reinterpret_cast<const char *>(fileData.data() + jsonStart),
                header.featureTableJSONByteLength);
            std::cout << "\nFeature Table JSON:\n"
                      << featureTableJSON << "\n\n";
        }

        // Check if data is Draco compressed
        // Look for Draco magic bytes in binary section
        if (binaryEnd > fileData.size()) {
            std::cerr << "Invalid binary data range\n";
            return false;
        }

        const uint8_t *binaryData = fileData.data() + binaryStart;
        size_t binarySize = header.featureTableBinaryByteLength;

        // Check for Draco magic number "DRACO" at the start of binary data
        bool isDraco = (binarySize > 5 && binaryData[0] == 'D' &&
                        binaryData[1] == 'R' && binaryData[2] == 'A' &&
                        binaryData[3] == 'C' && binaryData[4] == 'O');

        if (isDraco) {
            std::cout << "Detected Draco compression\n";
            return decodeDracoPointCloud(binaryData, binarySize, output);
        } else {
            std::cout
                << "No Draco compression detected, reading raw binary data\n";
            return decodeRawPointCloud(binaryData, binarySize, output);
        }
    }

    // Decode Draco compressed point cloud
    static bool decodeDracoPointCloud(const uint8_t *data, size_t size,
                                      PointCloudData &output) {
        // Create Draco decoder buffer
        draco::DecoderBuffer buffer;
        buffer.Init(reinterpret_cast<const char *>(data), size);

        // Decode point cloud
        draco::Decoder decoder;
        auto statusor = decoder.DecodePointCloudFromBuffer(&buffer);

        if (!statusor.ok()) {
            std::cerr << "Failed to decode Draco point cloud: "
                      << statusor.status().error_msg_string() << "\n";
            return false;
        }

        std::unique_ptr<draco::PointCloud> pointCloud =
            std::move(statusor).value();
        output.pointCount = pointCloud->num_points();

        std::cout << "Successfully decoded " << output.pointCount
                  << " points\n";

        // Extract position attribute
        const draco::PointAttribute *posAttr =
            pointCloud->GetNamedAttribute(draco::GeometryAttribute::POSITION);

        if (posAttr) {
            std::cout << "Found POSITION attribute\n";
            output.positions.reserve(output.pointCount * 3);

            for (uint32_t i = 0; i < output.pointCount; ++i) {
                float pos[3];
                posAttr->GetMappedValue(draco::PointIndex(i), pos);
                output.positions.push_back(pos[0]);
                output.positions.push_back(pos[1]);
                output.positions.push_back(pos[2]);
            }
        }

        // Extract color attribute
        const draco::PointAttribute *colorAttr =
            pointCloud->GetNamedAttribute(draco::GeometryAttribute::COLOR);

        if (colorAttr) {
            std::cout << "Found COLOR attribute\n";
            int numComponents = colorAttr->num_components();
            output.colors.reserve(output.pointCount * numComponents);

            for (uint32_t i = 0; i < output.pointCount; ++i) {
                if (numComponents == 3) {
                    uint8_t rgb[3];
                    colorAttr->GetMappedValue(draco::PointIndex(i), rgb);
                    output.colors.push_back(rgb[0]);
                    output.colors.push_back(rgb[1]);
                    output.colors.push_back(rgb[2]);
                } else if (numComponents == 4) {
                    uint8_t rgba[4];
                    colorAttr->GetMappedValue(draco::PointIndex(i), rgba);
                    output.colors.push_back(rgba[0]);
                    output.colors.push_back(rgba[1]);
                    output.colors.push_back(rgba[2]);
                    output.colors.push_back(rgba[3]);
                }
            }
        }

        // Extract normal attribute (if present)
        const draco::PointAttribute *normalAttr =
            pointCloud->GetNamedAttribute(draco::GeometryAttribute::NORMAL);

        if (normalAttr) {
            std::cout << "Found NORMAL attribute\n";
            output.normals.reserve(output.pointCount * 3);

            for (uint32_t i = 0; i < output.pointCount; ++i) {
                float normal[3];
                normalAttr->GetMappedValue(draco::PointIndex(i), normal);
                output.normals.push_back(normal[0]);
                output.normals.push_back(normal[1]);
                output.normals.push_back(normal[2]);
            }
        }

        // Extract generic attributes (like BATCH_ID)
        for (int32_t i = 0; i < pointCloud->num_attributes(); ++i) {
            const draco::PointAttribute *attr = pointCloud->attribute(i);

            // Check for custom attributes by metadata or unique_id
            if (attr->attribute_type() == draco::GeometryAttribute::GENERIC) {
                std::cout << "Found GENERIC attribute at index " << i << "\n";

                // Try to extract as batch IDs (assuming uint32_t)
                if (attr->num_components() == 1) {
                    output.batchIds.reserve(output.pointCount);
                    for (uint32_t j = 0; j < output.pointCount; ++j) {
                        uint32_t batchId;
                        attr->GetMappedValue(draco::PointIndex(j), &batchId);
                        output.batchIds.push_back(batchId);
                    }
                }
            }
        }

        return true;
    }

    // Decode raw (uncompressed) point cloud data
    static bool decodeRawPointCloud(const uint8_t *data, size_t size,
                                    PointCloudData &output) {
        // This is a simplified example - actual implementation depends on
        // the feature table JSON which specifies the data layout

        // Assuming positions are stored as float32[3] per point
        if (size % (3 * sizeof(float)) != 0) {
            std::cerr << "Invalid raw point data size\n";
            return false;
        }

        output.pointCount = size / (3 * sizeof(float));
        output.positions.resize(size / sizeof(float));

        std::memcpy(output.positions.data(), data, size);

        std::cout << "Decoded " << output.pointCount << " raw points\n";
        return true;
    }

    // Load PNTS file from disk
    static bool loadPntsFile(const std::string &filename,
                             PointCloudData &output) {
        std::ifstream file(filename, std::ios::binary | std::ios::ate);
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << filename << "\n";
            return false;
        }

        std::streamsize size = file.tellg();
        file.seekg(0, std::ios::beg);

        std::vector<uint8_t> buffer(size);
        if (!file.read(reinterpret_cast<char *>(buffer.data()), size)) {
            std::cerr << "Failed to read file\n";
            return false;
        }

        return decodePnts(buffer, output);
    }

    // Helper function to print point cloud statistics
    static void printStatistics(const PointCloudData &data) {
        std::cout << "\n=== Point Cloud Statistics ===\n";
        std::cout << "Point Count: " << data.pointCount << "\n";
        std::cout << "Has Positions: "
                  << (!data.positions.empty() ? "Yes" : "No") << "\n";
        std::cout << "Has Colors: " << (!data.colors.empty() ? "Yes" : "No")
                  << "\n";
        std::cout << "Has Normals: " << (!data.normals.empty() ? "Yes" : "No")
                  << "\n";
        std::cout << "Has Batch IDs: "
                  << (!data.batchIds.empty() ? "Yes" : "No") << "\n";

        // Print sample data
        if (!data.positions.empty() && data.pointCount > 0) {
            std::cout << "\nFirst point position: (" << data.positions[0]
                      << ", " << data.positions[1] << ", " << data.positions[2]
                      << ")\n";
        }

        if (!data.colors.empty() && data.pointCount > 0) {
            std::cout << "First point color: (" << (int)data.colors[0] << ", "
                      << (int)data.colors[1] << ", " << (int)data.colors[2]
                      << ")\n";
        }
    }
};

// Example usage
int main(int argc, char *argv[]) {
    std::cout << "=== Cesium PNTS Draco Decoder ===\n\n";

    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <pnts_file>\n";
        std::cout << "\nExample with simulated data:\n\n";

        // For demonstration, you would normally load an actual .pnts file
        // std::string filename = "example.pnts";
        // PointCloudData pointCloud;
        // if (CesiumPntsDecoder::loadPntsFile(filename, pointCloud)) {
        //     CesiumPntsDecoder::printStatistics(pointCloud);
        // }

        std::cout << "To use this decoder:\n";
        std::cout << "1. Link against Draco library\n";
        std::cout << "2. Provide a Cesium PNTS file path as argument\n";
        std::cout << "3. The decoder will detect and decompress Draco data "
                     "automatically\n\n";

        std::cout << "Build command example:\n";
        std::cout << "g++ -o decoder decoder.cpp -I/path/to/draco/src "
                     "-L/path/to/draco/lib -ldraco\n";

        return 1;
    }

    std::string filename = argv[1];
    PointCloudData pointCloud;

    if (CesiumPntsDecoder::loadPntsFile(filename, pointCloud)) {
        CesiumPntsDecoder::printStatistics(pointCloud);

        // Now you can use the decoded point cloud data
        // For example: render it, process it, export to another format, etc.

        std::cout << "\nDecoding successful!\n";
        return 0;
    } else {
        std::cerr << "Failed to decode PNTS file\n";
        return 1;
    }
}

