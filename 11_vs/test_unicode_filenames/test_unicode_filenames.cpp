// test_unicode_filenames.cpp : Use UNICODE filenames in non-unicode environment
//

#define WIN32_LEAN_AND_MEAN
#define _CRT_SECURE_NO_WARNINGS
#define _SILENCE_CXX17_CODECVT_HEADER_DEPRECATION_WARNING

#include <windows.h>
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
#endif

#include <locale>
#include <codecvt> // Required for std::wstring_convert
#include <fstream>
#include <filesystem>
namespace fs = std::filesystem;

TEST(t1, t1) {
    CONSOLE("Hello");
    EXPECT_EQ(1252, GetACP());

    char utf8_filepath[] = u8"out/Добрый День.txt";
    SetConsoleOutputCP(65001);
    CONSOLE(utf8_filepath);

    if (fs::is_directory("out"))
        fs::remove_all("out");
    fs::create_directories("out");

    {
        FILE *actual = fopen(utf8_filepath, "w");
        ASSERT_TRUE(actual);
        fclose(actual);
    }
    {
        std::wstring utf16_path = fs::path(utf8_filepath).wstring();
        std::ofstream actual{utf16_path};

        ASSERT_TRUE(actual);
        actual << utf8_filepath << "\n";
    }

    auto utf8_to_wchar = [](const std::string &utf8_str) -> std::wstring {
        // std::wstring_convert is deprecated in C++17, but often used for
        // convenience. Ensure proper error handling and consider
        // alternatives for new projects.
        std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;
        return converter.from_bytes(utf8_str);
    };

    {
        std::wstring utf16_path = utf8_to_wchar(utf8_filepath);
        std::ofstream actual{utf16_path};

        ASSERT_TRUE(actual);
        actual << utf8_filepath << "\n";
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
