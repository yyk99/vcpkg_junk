//
//
//

#include <gtest/gtest.h>
#include <iostream>
#include <string>

#include "vtkNamedColors.h"
#include "vtkSmartPointer.h"
#include "vtk_jsoncpp.h"

#include "../common/DebuggingConsole.h"

/// @brief 
class VTK_F : public testing::Test {
public:
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