//
//
//

#include <gtest/gtest.h>
#include <iostream>
#include <string>

#include <filesystem>
namespace fs = std::filesystem;

#include "vtkNamedColors.h"
#include "vtkSmartPointer.h"
#include "vtk_jsoncpp.h"

#include "../common/DebuggingConsole.h"

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

#define CONSOLE_TE(x) CONSOLE_T(#x << " : " << (x))

/// @brief
class VTK_F : public testing::Test {
public:
    const char *test_case_name() const {
        // https://google.github.io/googletest/advanced.html#getting-the-current-tests-name
        // Gets information about the currently running test.
        // Do NOT delete the returned object - it's managed by the UnitTest
        // class.
        const testing::TestInfo *const test_info =
            testing::UnitTest::GetInstance()->current_test_info();

        return test_info->test_case_name();
    }

    const char *test_name() const {
        // https://google.github.io/googletest/advanced.html#getting-the-current-tests-name
        // Gets information about the currently running test.
        // Do NOT delete the returned object - it's managed by the UnitTest
        // class.
        const testing::TestInfo *const test_info =
            testing::UnitTest::GetInstance()->current_test_info();

        return test_info->name();
    }

    fs::path create_workspace() {
        auto ws = fs::path("out") / test_case_name() / test_name();
        if (fs::is_directory(ws))
            fs::remove_all(ws);
        std::error_code err;
        fs::create_directories(ws, err);

        CONSOLE("ws = " << fs::absolute(ws).string());

        return ws;
    }

    //------------------------------------------------------------------------------
    std::vector<std::string> ParseColorNames(const std::string &colorNames) {
        // The delimiter for a color.
        const std::string colorDelimiter = "\n";
        std::vector<std::string> cn;
        size_t start = 0;
        size_t end = colorNames.find(colorDelimiter);
        while (end != std::string::npos) {
            cn.emplace_back(colorNames.substr(start, end - start));
            start = end + 1;
            end = colorNames.find(colorDelimiter, start);
        }
        // Get the last color.
        if (!colorNames.empty()) {
            cn.emplace_back(
                colorNames.substr(start, colorNames.size() - start));
        }
        return cn;
    }
};

/// @brief
/// @param
/// @param
TEST_F(VTK_F, vtkNamedColors_GetColorNames) {
    vtkSmartPointer<vtkNamedColors> sot =
        vtkSmartPointer<vtkNamedColors>::New();
    std::string name;
    {
        auto actual_colors = sot->GetColorNames();
        CONSOLE_EVAL(actual_colors);
        ASSERT_FALSE(actual_colors.empty());
        ASSERT_NE(0, actual_colors.size());
    }
    {
        auto actual_colors = ParseColorNames(sot->GetColorNames());
        CONSOLE_EVAL(actual_colors.size());
        ASSERT_NE(0, actual_colors.size());
    }
}

/// @brief
/// @param
/// @param
TEST_F(VTK_F, json_parse) {
    Json::CharReaderBuilder sot_builder;

    auto sot = std::unique_ptr<Json::CharReader>(sot_builder.newCharReader());
    ASSERT_TRUE(sot);

    {
        const char text_json[] = R"({"a" = 10, "b" = "hello"})";
        Json::Value actual_root;
        Json::String actual_err;
        bool ok = sot->parse(text_json, text_json + strlen(text_json),
                             &actual_root, &actual_err);
        CONSOLE_EVAL(actual_err);
        ASSERT_FALSE(ok);
    }
    {
        const char text_json[] = R"({"a" : 10, "b" : "hello"})";
        Json::Value actual_root;
        Json::String actual_err;
        bool ok = sot->parse(text_json, text_json + strlen(text_json),
                             &actual_root, &actual_err);
        CONSOLE_EVAL(actual_err);
        ASSERT_TRUE(ok);
        ASSERT_TRUE(actual_err.empty());

        EXPECT_TRUE(actual_root.isObject());
        auto actual_a = actual_root.get("a", Json::Value(0));
        EXPECT_TRUE(actual_a.isInt());
        EXPECT_EQ(10, actual_a.asInt());

        auto actual_b = actual_root.get("b", Json::Value(""));
        EXPECT_TRUE(actual_b.isString());
        EXPECT_EQ(std::string("hello"), actual_b.asString());
    }
}

#include <vtkAppendPolyData.h>
#include <vtkCesium3DTilesWriter.h>
#include <vtkFieldData.h>
#include <vtkLogger.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkOBJReader.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>

//------------------------------------------------------------------------------
std::array<double, 3> ReadOBJOffset(const char *comment) {
    std::array<double, 3> translation = {0, 0, 0};
    if (comment) {
        std::istringstream istr(comment);
        std::array<std::string, 3> axesNames = {"x", "y", "z"};
        for (int i = 0; i < 3; ++i) {
            std::string axis;
            std::string s;
            istr >> axis >> s >> translation[i];
            if (istr.fail()) {
                vtkLog(WARNING,
                       "Cannot read axis " << axesNames[i] << " from comment.");
            }
            if (axis != axesNames[i]) {
                vtkLog(WARNING,
                       "Invalid axis " << axesNames[i] << ": " << axis);
            }
        }
    } else {
        vtkLog(WARNING, "nullptr comment.");
    }
    return translation;
}

vtkSmartPointer<vtkPolyData> ReadOBJMesh(int numberOfBuildings,
                                         int vtkNotUsed(lod),
                                         const std::vector<std::string> &files,
                                         std::array<double, 3> &fileOffset) {
    vtkNew<vtkAppendPolyData> append;
    for (size_t i = 0;
         i < files.size() && i < static_cast<size_t>(numberOfBuildings); ++i) {
        vtkNew<vtkOBJReader> reader;
        reader->SetFileName(files[i].c_str());
        reader->Update();
        if (i == 0) {
            fileOffset = ReadOBJOffset(reader->GetComment());
        }
        auto polyData = reader->GetOutput();
        append->AddInputDataObject(polyData);
    }
    append->Update();
    return append->GetOutput();
}

//------------------------------------------------------------------------------
void SetField(vtkDataObject *obj, const char *name, const char *value) {
    vtkFieldData *fd = obj->GetFieldData();
    if (!fd) {
        vtkNew<vtkFieldData> newfd;
        obj->SetFieldData(newfd);
        fd = newfd;
    }
    vtkNew<vtkStringArray> sa;
    sa->SetNumberOfTuples(1);
    sa->SetValue(0, value);
    sa->SetName(name);
    fd->AddArray(sa);
}

//------------------------------------------------------------------------------
std::string GetOBJTextureFileName(const std::string &file) {
    std::string fileNoExt = fs::path(file).stem().string();
    std::string textureFileName = fileNoExt + ".png";
    return fs::is_regular_file(textureFileName) ? textureFileName
                                                : std::string();
}

vtkSmartPointer<vtkMultiBlockDataSet>
ReadOBJBuildings(int numberOfBuildings, int vtkNotUsed(lod),
                 std::vector<std::string> const &files,
                 std::array<double, 3> &fileOffset) {
    auto root = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    for (size_t i = 0;
         i < files.size() && i < static_cast<size_t>(numberOfBuildings); ++i) {
        vtkNew<vtkOBJReader> reader;
        reader->SetFileName(files[i].c_str());
        reader->Update();
        if (i == 0) {
            fileOffset = ReadOBJOffset(reader->GetComment());
        }
        auto polyData = reader->GetOutput();
        std::string textureFileName = GetOBJTextureFileName(files[i]);
        if (!textureFileName.empty()) {
            SetField(polyData, "texture_uri", textureFileName.c_str());
        }
        auto building = vtkSmartPointer<vtkMultiBlockDataSet>::New();
        building->SetBlock(0, polyData);
        root->SetBlock(root->GetNumberOfBlocks(), building);
    }
    return root;
}

TEST_F(VTK_F, vtk_t00) {
    auto ws = create_workspace();

    using namespace vtksys;

    std::array<double, 3> fileOffset;
    vtkSmartPointer<vtkMultiBlockDataSet> actual = ReadOBJBuildings(
        1, 42, {R"(V:\Asset Importer Test Models\models\OBJ\box.obj)"},
        fileOffset);
    ASSERT_TRUE(actual);
    CONSOLE_TE(actual->GetNumberOfBlocks());
}

/// @brief
/// @see
/// https://gitlab.kitware.com/vtk/vtk/-/blob/v9.3.0/IO/Cesium3DTiles/Testing/Cxx/TestCesium3DTilesWriter.cxx
/// @param
/// @param
TEST_F(VTK_F, Cesium3DTilesWriter) {
    auto ws = create_workspace();

    using namespace vtksys;

    std::string textureBaseDirectory = ws.string(); // TODO: set
    std::array<double, 3> fileOffset;

    const std::vector<std::string> input;
    int inputType = vtkCesium3DTilesWriter::Mesh;
    bool addColor{};
    std::string output = ws.string();
    bool contentGLTF = true;
    int numberOfBuildings{};
    int buildingsPerTile{};
    int lod{};
    std::vector<double> inputOffset;
    bool saveTiles = true;
    bool saveTextures = true;
    std::string crs;
    int utmZone = 19; // Maine
    char utmHemisphere = 'N';

    vtkSmartPointer<vtkMultiBlockDataSet> mbData = ReadOBJBuildings(
        1, 42, {R"(V:\Asset Importer Test Models\models\OBJ\box.obj)"},
        fileOffset);
    ASSERT_TRUE(mbData);
    EXPECT_EQ(1, mbData->GetNumberOfBlocks());

    vtkNew<vtkCesium3DTilesWriter> writer;

    writer->SetInputDataObject(mbData);

    writer->SetContentGLTF(contentGLTF);
    writer->ContentGLTFSaveGLBOff();
    writer->SetInputType(inputType);
    writer->SetDirectoryName(output.c_str());
    writer->SetTextureBaseDirectory(textureBaseDirectory.c_str());
    writer->SetOffset(fileOffset.data());
    writer->SetSaveTextures(saveTextures);
    writer->SetNumberOfFeaturesPerTile(buildingsPerTile);
    writer->SetSaveTiles(saveTiles);
    if (crs.empty()) {
        std::ostringstream ostr;
        ostr << "+proj=utm +zone=" << utmZone
             << (utmHemisphere == 'S' ? "+south" : "");
        crs = ostr.str();
    }
    writer->SetCRS(crs.c_str());
    int rc = writer->Write();

    CONSOLE_T("rc = " << rc);
}