//
//
//

#include <gtest/gtest.h>
#include <iostream>
#include <string>

#include "vtkNamedColors.h"
#include "vtkSmartPointer.h"

#define _DEBUG 1
#include "../common/DebuggingConsole.h"

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

///
TEST_F(VTK_F, t0) {
    CONSOLE("Hello...");

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