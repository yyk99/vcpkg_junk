//
//
//

#include <filesystem>
#include <fstream>
#include <string>

namespace fs = std::filesystem;

#include "test_draco.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdint.h>
#include <vector>

// Include Draco decoder headers
#include "draco/compression/decode.h"
#include "draco/point_cloud/point_cloud.h"

#include <json-c/json.h>
#include <json-c/printbuf.h>

#include "../common/CONSOLE.h"

#define VERBOSE 0

namespace cesium_pnts {
#if 0
    }
#endif

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

struct PointCloudBatchData {
    std::vector<uint8_t> batch_data;
};

struct PointCloudData {
    std::vector<float> positions;    // x, y, z coordinates
    std::vector<uint8_t> colors;     // RGB or RGBA colors
    std::vector<float> normals;      // Optional normals
    std::vector<uint8_t> batch_ids;  // Optional batch IDs. If empty, then the
                                     // id mapping is trivial.
    std::vector<PointCloudBatchData> attributes;
    uint32_t pointCount;

    PointCloudData() : pointCount(0) {}
};

struct batch_array_t {
    std::string attr;              // The attribute name, e.g. "name"
    std::vector<std::string> vals; // The attribute values in a array (json text)
};

struct batch_reference_t {
    std::string attr;
    uint32_t byte_offset;
    draco::DataType component_type; /* BYTE, SHORT, ...*/
    int32_t number_of_components;   /* SCALAR == 1, VEC2 == 2, VEC3 == 3, VEC4 == 4 */
};

struct batch_header_t {
    std::vector<std::variant<batch_array_t, batch_reference_t>> attr;
};

class CesiumPntsDecoder {
public:
    static draco::DataType batch_reference_component_type(const char *ct) {
        // clang-format off
#define R(s, t) if (!strcmp(ct, #s)) return draco::DT_##t
        // clang-format on

        R(BYTE, INT8);
        R(UNSIGNED_BYTE, UINT8);
        R(SHORT, INT16);
        R(UNSIGNED_SHORT, UINT16);
        R(INT, INT32);
        R(UNSIGNED_INT, UINT32);
        R(FLOAT, FLOAT32);
        R(DOUBLE, FLOAT64);
        return draco::DT_INVALID;
#undef R
    }

    PntsHeader header;
    batch_header_t batch_header;
    std::string last_error;

public:
    // Decode a PNTS file from memory buffer
    bool decodePnts(const std::vector<uint8_t> &fileData, PointCloudData &output) {
        if (fileData.size() < sizeof(PntsHeader)) {
            last_error = "File too small to contain PNTS header";
            return false;
        }

        std::memcpy(&header, fileData.data(), sizeof(PntsHeader));

        // Validate magic number
        if (std::strncmp(header.magic, "pnts", 4) != 0) {
            last_error = "Invalid PNTS magic number";
            return false;
        }
#if VERBOSE
        CONSOLE("PNTS Header:");
        CONSOLE_EVAL(header.version);
        CONSOLE_EVAL(header.byteLength);
        CONSOLE_EVAL(header.featureTableJSONByteLength);
        CONSOLE_EVAL(header.featureTableBinaryByteLength);
        CONSOLE_EVAL(header.batchTableJSONByteLength);
        CONSOLE_EVAL(header.batchTableBinaryByteLength);
#endif
        // Calculate offsets
        size_t offset = sizeof(PntsHeader);
        size_t jsonStart = offset;
        size_t jsonEnd = jsonStart + header.featureTableJSONByteLength;
        size_t binaryStart = jsonEnd;
        size_t binaryEnd = binaryStart + header.featureTableBinaryByteLength;

        // Extract feature table JSON (optional - parse to get metadata)
        std::string featureTableJSON;
        if (header.featureTableJSONByteLength > 0) {
            featureTableJSON = std::string(reinterpret_cast<const char *>(fileData.data() + jsonStart),
                                           header.featureTableJSONByteLength);
#if VERBOSE
            CONSOLE_EVAL(featureTableJSON);
#endif
        }

        // Check if data is Draco compressed
        // Look for Draco magic bytes in binary section
        if (binaryEnd > fileData.size()) {
            last_error = "Invalid binary data range";
            return false;
        }

        const uint8_t *binaryData = fileData.data() + binaryStart;
        size_t binarySize = header.featureTableBinaryByteLength;

        /* parse batch part */
        std::string batchTableJSON;
        if (header.batchTableJSONByteLength) {
            const uint8_t *batchData = fileData.data() + binaryStart + binarySize;
            batchTableJSON.assign((const char *)batchData, header.batchTableJSONByteLength);
#if VERBOSE
            CONSOLE_EVAL(batchTableJSON);
#endif
        }

        if (!batchTableJSON.empty()) {
            auto pos = binaryData + binarySize + header.batchTableJSONByteLength;
            if (parse_batch_table_json(batchTableJSON.c_str(), batch_header)) {
                // ...
            }
        }

        // Check for Draco magic number "DRACO" at the start of binary data
        bool isDraco = (binarySize > 5 && binaryData[0] == 'D' && binaryData[1] == 'R' && binaryData[2] == 'A' &&
                        binaryData[3] == 'C' && binaryData[4] == 'O');

        if (isDraco)
            return decodeDracoPointCloud(binaryData, binarySize, output);
        /* non-draco */
        return decodeRawPointCloud(featureTableJSON, batchTableJSON, binaryData, binarySize, output);
    }

    // Decode Draco compressed point cloud
    bool decodeDracoPointCloud(const uint8_t *data, size_t size, PointCloudData &output) {
        // Create Draco decoder buffer
        draco::DecoderBuffer buffer;
        buffer.Init(reinterpret_cast<const char *>(data), size);

        // Decode point cloud
        draco::Decoder decoder;
        auto statusor = decoder.DecodePointCloudFromBuffer(&buffer);

        if (!statusor.ok()) {
            std::cerr << "Failed to decode Draco point cloud: " << statusor.status().error_msg_string() << "\n";
            return false;
        }

        std::unique_ptr<draco::PointCloud> pointCloud = std::move(statusor).value();
        output.pointCount = pointCloud->num_points();

        // Extract position attribute
        const draco::PointAttribute *posAttr = pointCloud->GetNamedAttribute(draco::GeometryAttribute::POSITION);

        if (posAttr) {
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
        const draco::PointAttribute *colorAttr = pointCloud->GetNamedAttribute(draco::GeometryAttribute::COLOR);

        if (colorAttr) {
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
        const draco::PointAttribute *normalAttr = pointCloud->GetNamedAttribute(draco::GeometryAttribute::NORMAL);

        if (normalAttr) {
            output.normals.reserve(output.pointCount * 3);

            for (uint32_t i = 0; i < output.pointCount; ++i) {
                float normal[3];
                normalAttr->GetMappedValue(draco::PointIndex(i), normal);
                output.normals.push_back(normal[0]);
                output.normals.push_back(normal[1]);
                output.normals.push_back(normal[2]);
            }
        }

        // Extract generic attributes
        std::vector<int32_t> generic_attributes_ids;
        for (int32_t i = 0; i < pointCloud->num_attributes(); ++i) {
            const draco::PointAttribute *attr = pointCloud->attribute(i);
            if (attr->attribute_type() == draco::GeometryAttribute::GENERIC) {
                generic_attributes_ids.push_back(i);
            }
        }
        output.attributes.reserve(generic_attributes_ids.size());
        // Extract the entire attribute data
        for (auto i : generic_attributes_ids) {
            const draco::PointAttribute *attr = pointCloud->attribute(i);
            output.attributes.push_back({});
            if (attr->attribute_type() == draco::GeometryAttribute::GENERIC) {
                size_t attribute_size = attr->num_components() * draco::DataTypeLength(attr->data_type());
                PointCloudBatchData &attr_batch = output.attributes.back();
                attr_batch.batch_data.resize(attribute_size * output.pointCount);
                for (uint32_t j = 0; j < output.pointCount; ++j) {
                    attr->GetMappedValue(draco::PointIndex(j), &attr_batch.batch_data[j * attribute_size]);
                }
            }
        }

        return true;
    }

    bool false_because(std::string s) {
        last_error = s;
        return false;
    }

    struct feature_table_props_t {
        std::optional<unsigned> POINTS_LENGTH;
        std::optional<unsigned> POSITION_byteOffset;
        std::optional<unsigned> NORMAL_byteOffset;
        std::optional<unsigned> RGBA_byteOffset;
        std::optional<unsigned> POSITION_QUANTIZED_byteOffset;
        std::optional<unsigned> RGB_byteOffset;
        std::optional<unsigned> BATCH_ID_byteOffset;
        std::optional<draco::DataType> BATCH_ID_componentType;
        std::vector<float> RTC_CENTER;
        std::map<std::string, std::string> extensions;
    };

    static int batch_reference_number_of_components(const char *t) {
        if (!strcmp(t, "SCALAR"))
            return 1;
        if (!strcmp(t, "VEC2"))
            return 2;
        if (!strcmp(t, "VEC3"))
            return 3;
        if (!strcmp(t, "VEC4"))
            return 4;
        return 0; // invalid!
    }

    bool parse_feature_table_json(const char *json_text, feature_table_props_t &props) {
        auto root = json_tokener_parse(json_text);
        if (!root)
            return false_because("Cannot parse featureTableJSON");

        // destroy `root` object on exit of the scope
        std::unique_ptr<json_object, std::function<void(json_object *)>> on_return(
            root, [](json_object *r) { json_object_put(r); });

        {
            auto o = json_object_object_get(root, "POINTS_LENGTH");
            if (!o)
                return false_because("POINTS_LENGTH is not defined");
            props.POINTS_LENGTH = json_object_get_int(o);
        }

        if (auto eo = json_object_object_get(root, "extensions")) {
            if (auto e_draco = json_object_object_get(eo, "3DTILES_draco_point_compression")) {
                // TODO: parse the draco extensions parameters
                props.extensions.insert({"3DTILES_draco_point_compression", "1"});
            }
        }

        if (auto po = json_object_object_get(root, "POSITION")) {
            unsigned byteOffset = 0;
            if (auto bo = json_object_object_get(root, "byteOffset"))
                byteOffset = json_object_get_int(bo);
            props.POSITION_byteOffset = byteOffset;
        }
        if (auto o = json_object_object_get(root, "RTC_CENTER")) {
            if (json_type_array != json_object_get_type(o))
                return false_because("RTC_CENTER must be array");
            array_list *arr = json_object_get_array(o);
            if (arr->length != 3)
                return false_because("RTC_CENTER must be array of 3");
            props.RTC_CENTER.resize(3);
            for (int i = 0; i < arr->length; ++i)
                props.RTC_CENTER[i] = static_cast<float>(json_object_get_double((json_object *)arr->array[i]));
        }
        if (auto po = json_object_object_get(root, "RGB")) {
            unsigned byteOffset = 0;
            if (auto bo = json_object_object_get(po, "byteOffset"))
                byteOffset = json_object_get_int(bo);
            props.RGB_byteOffset = byteOffset;
        }
        if (auto po = json_object_object_get(root, "RGBA")) {
            unsigned byteOffset = 0;
            if (auto bo = json_object_object_get(po, "byteOffset"))
                byteOffset = json_object_get_int(bo);
            props.RGBA_byteOffset = byteOffset;
        }
        if (auto po = json_object_object_get(root, "BATCH_ID")) {
            unsigned byteOffset = 0;
            if (auto bo = json_object_object_get(po, "byteOffset"))
                byteOffset = json_object_get_int(bo);
            props.BATCH_ID_byteOffset = byteOffset;
            if (auto ct = json_object_object_get(po, "componentType")) {
                const char *component_type = json_object_get_string(ct);
                props.BATCH_ID_componentType = batch_reference_component_type(component_type);
            }
        }
        if (auto po = json_object_object_get(root, "NORMAL")) {
            unsigned byteOffset = 0;
            if (auto bo = json_object_object_get(po, "byteOffset"))
                byteOffset = json_object_get_int(bo);
            props.NORMAL_byteOffset = byteOffset;
        }

        if (!props.POSITION_byteOffset.has_value() && !props.POSITION_QUANTIZED_byteOffset.has_value())
            return false_because("Neither POSITION nor POSITION_QUANTIZED defined");

        return true;
    }

    // Decode raw (uncompressed) point cloud data
    bool decodeRawPointCloud(std::string const &featureTableJSON, std::string const &batchTableJSON,
                             const uint8_t *data, size_t size, PointCloudData &output) {

        feature_table_props_t features;

        if (!parse_feature_table_json(featureTableJSON.c_str(), features))
            return false;

        output.pointCount = features.POINTS_LENGTH.value();
        if (features.POSITION_byteOffset.has_value()) {
            auto pos = reinterpret_cast<float const *>(data + features.POSITION_byteOffset.value());
            output.positions.assign(pos, pos + 3 * output.pointCount);
        } else {
            // TODO: implement POSITION_QUANTIZED support
            return false_because("only POSITION currently supported");
        }

        if (features.RGB_byteOffset.has_value()) {
            auto pos = reinterpret_cast<uint8_t const *>(data + features.RGB_byteOffset.value());
            output.colors.assign(pos, pos + 3 * output.pointCount);
        } else if (features.RGBA_byteOffset.has_value()) {
            auto pos = reinterpret_cast<uint8_t const *>(data + features.RGBA_byteOffset.value());
            output.colors.assign(pos, pos + 4 * output.pointCount);
        }

        if (features.BATCH_ID_byteOffset.has_value()) {
            auto pos = reinterpret_cast<uint8_t const *>(data + features.BATCH_ID_byteOffset.value());
            auto type_length = 4;
            if (features.BATCH_ID_componentType.has_value())
                type_length = draco::DataTypeLength(features.BATCH_ID_componentType.value());
            output.batch_ids.assign(pos, pos + type_length * output.pointCount);
        }

        if (features.NORMAL_byteOffset.has_value()) {
            auto pos = reinterpret_cast<float const *>(data + features.NORMAL_byteOffset.value());
            output.normals.assign(pos, pos + 3 * output.pointCount);
        }

        output.attributes.reserve(batch_header.attr.size());
        // Extract the entire attribute data
        for (auto &ap : batch_header.attr) {
            output.attributes.push_back({});
        }

        return true;
    }

    /// @brief
    /// @param batchTableJSON
    /// @param batch_header
    /// @return true if success
    bool parse_batch_table_json(const char *batchTableJSON, batch_header_t &batch_header) {
        auto root = json_tokener_parse(batchTableJSON);
        if (!root)
            return false_because("Cannot parse batchTableJSON");

        // destroy `root` object on exit of the scope
        std::unique_ptr<json_object, std::function<void(json_object *)>> on_return(
            root, [](json_object *r) { json_object_put(r); });

        auto it = json_object_iter_begin(root);
        auto itEnd = json_object_iter_end(root);
        while (!json_object_iter_equal(&it, &itEnd)) {
            const char *key = json_object_iter_peek_name(&it);
            auto val = json_object_iter_peek_value(&it);

            if (!strcmp("extensions", key)) {
                // extensions will be parsed separately
                ; // TODO: skip for now
            } else if (json_object_get_type(val) == json_type_object) {
                /* we are parsing a reference */
                /* e.g.
                 * {"byteOffset":152337,"componentType":"UNSIGNED_SHORT","type":"SCALAR"}
                 */

                batch_reference_t br{};
                br.attr = key;
                if (auto bo = json_object_object_get(val, "byteOffset")) {
                    br.byte_offset = json_object_get_int(bo);
                }
                if (auto ct = json_object_object_get(val, "componentType")) {
                    br.component_type = batch_reference_component_type(json_object_get_string(ct));
                }
                if (auto t = json_object_object_get(val, "type")) {
                    br.number_of_components = batch_reference_number_of_components(json_object_get_string(t));
                }

                batch_header.attr.push_back(br);
            } else if (json_object_get_type(val) == json_type_array) {
                /* we are parsing an array */
                batch_array_t batch_arr{};
                batch_arr.attr = key;
                /* split the val into items and save as individual json strings */
                array_list *arr = json_object_get_array(val);
                for (int i = 0; i < arr->length; ++i) {
                    const char *json_string = json_object_to_json_string((json_object *)arr->array[i]);
                    batch_arr.vals.push_back(std::string(json_string));
                }
                batch_header.attr.push_back(batch_arr);
            } else {
                CONSOLE("Not expected type: " << json_object_get_type(val));
            }

            json_object_iter_next(&it);
        }
        return true;
    }

    // Load PNTS file from disk
    bool loadPntsFile(const std::string &filename, PointCloudData &output) {
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
    void printStatistics(const PointCloudData &data) {
        if (batch_header.attr.size()) {
            CONSOLE("=== Batch REFERENCE definitions: ===");
            for (auto &ap : batch_header.attr) {
                if (std::holds_alternative<batch_reference_t>(ap)) {
                    auto &rp = std::get<batch_reference_t>(ap);
                    CONSOLE(rp.attr << " { off:" << rp.byte_offset << ", type:" << rp.component_type
                                    << ", number_of_components:" << rp.number_of_components << " }");
                } else if (std::holds_alternative<batch_array_t>(ap)) {
                    auto &rp = std::get<batch_array_t>(ap);
                    std::string ss = "[";
                    for (auto &s : rp.vals)
                        ss += " " + s;
                    ss += "]";
                    CONSOLE(rp.attr << " " << ss);
                }
            }
        }

        CONSOLE("=== Point Cloud Statistics ===");
        CONSOLE("Point Count: " << data.pointCount);
        CONSOLE("Has Positions: " << (!data.positions.empty() ? "Yes" : "No"));
        CONSOLE("Has Colors: " << (!data.colors.empty() ? "Yes" : "No"));
        CONSOLE("Has Normals: " << (!data.normals.empty() ? "Yes" : "No"));
        CONSOLE("Has Batch IDs: " << (!data.batch_ids.empty() ? "Yes" : "No"));

        // Print sample data
        if (!data.positions.empty() && data.pointCount > 0) {
            CONSOLE("First point position: (" << data.positions[0] << ", " << data.positions[1] << ", "
                                              << data.positions[2] << ")");
        }

        if (!data.colors.empty() && data.pointCount > 0) {
            CONSOLE("First point color: (" << (int)data.colors[0] << ", " << (int)data.colors[1] << ", "
                                           << (int)data.colors[2] << ")");
        }
    }
};

} // namespace cesium_pnts

class DracoF : public testing::Test {
protected:
    std::string test_name() const { return ::testing::UnitTest::GetInstance()->current_test_info()->name(); }

    std::string test_case_name() const {
        return ::testing::UnitTest::GetInstance()->current_test_info()->test_case_name();
    }

    fs::path create_ws() {
        auto ws = fs::absolute("out") / test_case_name() / test_name();
        fs::remove_all(ws);
        fs::create_directories(ws);
        return ws;
    }

    void save_as(std::string const &text, fs::path const &location) {
        std::ofstream os(location);
        if (!os.good())
            throw std::runtime_error("Cannot create file: " + location.string());
        os << text;
        os.close();
    }

    auto test_data(const char *relative_path) {
        auto test_data_dir = fs::absolute(__FILE__).parent_path() / "test_data";
        return test_data_dir / relative_path;
    }

    bool has_failure() const { return ::testing::Test::HasFailure(); }

protected:
    void SetUp() override {}
};

/// @brief
/// @param --gtest_filter=DracoF.t0
/// @param
TEST_F(DracoF, t0) {
    auto zero_pnts = test_data("cesium/pnts/0-draco.pnts");
    ASSERT_TRUE(fs::is_regular_file(zero_pnts)) << "File: " << zero_pnts.string();
}

/// @brief
/// @param --gtest_filter=DracoF.t0_no_draco
/// @param
TEST_F(DracoF, t0_no_draco) {
    auto zero_pnts = test_data("cesium/pnts/0.pnts");
    ASSERT_TRUE(fs::is_regular_file(zero_pnts)) << "File: " << zero_pnts.string();
}

/// @brief
/// @param --gtest_filter=DracoF.decodePnts
/// @param
TEST_F(DracoF, decodePnts) {
    using namespace cesium_pnts;

    auto zero_pnts = test_data("cesium/pnts/0-draco.pnts");
    ASSERT_TRUE(fs::is_regular_file(zero_pnts)) << "File: " << zero_pnts.string();

    CesiumPntsDecoder sot;
    PointCloudData actual;
    ASSERT_TRUE(sot.loadPntsFile(zero_pnts.string(), actual));
#if VERBOSE
    sot.printStatistics(actual);
#endif
    EXPECT_EQ(50779, actual.pointCount);
    EXPECT_EQ(0, actual.normals.size());
    EXPECT_EQ(50779 * 3, actual.colors.size());

    ASSERT_TRUE(actual.batch_ids.empty());
    ASSERT_EQ(3, actual.attributes.size());
    ASSERT_EQ(3, sot.batch_header.attr.size());
    EXPECT_EQ("Intensity", std::get<1>(sot.batch_header.attr[0]).attr);
    // EXPECT_EQ(0, sot.batch_header.references[0].byte_offset);
    // EXPECT_EQ(draco::DT_UINT16, sot.batch_header.references[0].component_type);
    // EXPECT_EQ(1, sot.batch_header.references[0].number_of_components);
    // ASSERT_EQ(50779 * 2 * 1, actual.attributes[0].batch_data.size());
#if VERBOSE
    {
        uint16_t *uint16_ptr = reinterpret_cast<uint16_t *>(actual.attributes[0].batch_data.data());
        CONSOLE(std::hex);
        CONSOLE_EVAL(uint16_ptr[0]);
        CONSOLE_EVAL(uint16_ptr[1]);
        CONSOLE_EVAL(uint16_ptr[2]);
        CONSOLE_EVAL(uint16_ptr[50779 - 1]);
    }
#endif
#if 0 // TODO:
    EXPECT_EQ("NumberOfReturns", sot.batch_header.references[1].attr);
    EXPECT_EQ(0, sot.batch_header.references[1].byte_offset);
    EXPECT_EQ(draco::DT_UINT8, sot.batch_header.references[1].component_type);
    EXPECT_EQ(1, sot.batch_header.references[1].number_of_components);
    ASSERT_EQ(50779 * 1 * 1, actual.attributes[1].batch_data.size());

    EXPECT_EQ("PointSourceID", sot.batch_header.references[2].attr);
    EXPECT_EQ(0, sot.batch_header.references[2].byte_offset);
    EXPECT_EQ(draco::DT_UINT16, sot.batch_header.references[2].component_type);
    EXPECT_EQ(1, sot.batch_header.references[2].number_of_components);
    ASSERT_EQ(50779 * 2 * 1, actual.attributes[2].batch_data.size());
#endif
    EXPECT_EQ(50779 * 3, actual.positions.size());
}

/// @brief
/// @param --gtest_filter=DracoF.decodePnts_no_draco
/// @param
TEST_F(DracoF, decodePnts_no_draco) {
    using namespace cesium_pnts;

    auto zero_pnts = test_data("cesium/pnts/0.pnts");
    ASSERT_TRUE(fs::is_regular_file(zero_pnts)) << "File: " << zero_pnts.string();

    CesiumPntsDecoder sot;
    PointCloudData actual;
    ASSERT_TRUE(sot.loadPntsFile(zero_pnts.string(), actual));

    sot.printStatistics(actual);

    EXPECT_EQ(1016024, sot.header.byteLength);
    EXPECT_EQ(156, sot.header.featureTableJSONByteLength);
    EXPECT_EQ(761688, sot.header.featureTableBinaryByteLength);
    EXPECT_EQ(256, sot.header.batchTableJSONByteLength);
    EXPECT_EQ(253896, sot.header.batchTableBinaryByteLength);

    EXPECT_EQ(50779, actual.pointCount);
    EXPECT_EQ(0, actual.normals.size());
    EXPECT_EQ(50779 * 3, actual.colors.size());
    EXPECT_EQ(50779 * 3, actual.positions.size());

    EXPECT_EQ(0, actual.batch_ids.size());
}

/// @brief Mixed case of batch data (both array & references)
/// @param --gtest_filter=DracoF.decodePnts_no_draco_2
/// @param
TEST_F(DracoF, decodePnts_no_draco_2) {
    using namespace cesium_pnts;

    auto zero_pnts = test_data("cesium/pnts/PointCloudBatched/pointCloudBatched.pnts");
    ASSERT_TRUE(fs::is_regular_file(zero_pnts)) << "File: " << zero_pnts.string();

    CesiumPntsDecoder sot;
    PointCloudData actual;
    ASSERT_TRUE(sot.loadPntsFile(zero_pnts.string(), actual));

    sot.printStatistics(actual);

    EXPECT_EQ(25632, sot.header.byteLength);
    EXPECT_EQ(236, sot.header.featureTableJSONByteLength);
    EXPECT_EQ(25000, sot.header.featureTableBinaryByteLength);
    EXPECT_EQ(240, sot.header.batchTableJSONByteLength);
    EXPECT_EQ(128, sot.header.batchTableBinaryByteLength);

    EXPECT_EQ(3, sot.batch_header.attr.size());

    EXPECT_EQ(1000, actual.pointCount);
    EXPECT_EQ(1000 * 3, actual.normals.size());
    EXPECT_EQ(0, actual.colors.size());
    EXPECT_EQ(1000 * 3, actual.positions.size());
    EXPECT_EQ(1000, actual.batch_ids.size());
    EXPECT_EQ(3, actual.attributes.size());
}
