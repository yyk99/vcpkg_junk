//
//
//

#include <chrono>
#include <map>
#include <vector>

#include <fmt/core.h>

//#include <fmt/chrono.h>
//#include <fmt/color.h>
//#include <fmt/format.h>
//#include <fmt/ranges.h>


#if 0
struct Point {
    int x, y;
};

// Custom formatter example (add this outside main)
template <> struct fmt::formatter<Point> {
    constexpr auto parse(format_parse_context &ctx) -> decltype(ctx.begin()) {
        return ctx.end();
    }

    template <typename FormatContext>
    auto format(const Point &p, FormatContext &ctx) -> decltype(ctx.out()) {
        return fmt::format_to(ctx.out(), "Point({}, {})", p.x, p.y);
    }
};
#endif

int main() {
    // 1. Basic formatting (replaces printf/iostream)
    fmt::print("Hello, {}!\n", "World");
    fmt::print("Number: {}, Float: {:.2f}\n", 42, 3.14159);

    // 2. String formatting
    std::string msg = fmt::format("User {} has {} points", "Alice", 1500);
    fmt::print("{}\n", msg);

    // 3. Positional arguments
    fmt::print("{1} comes before {0}\n", "second", "first");
    fmt::print("{0} {0} {1}\n", "hello", "world");

    // 4. Named arguments
    fmt::print("Point: ({x}, {y})\n", fmt::arg("x", 10), fmt::arg("y", 20));

    // 5. Number formatting
    fmt::print("Binary: {:b}\n", 42);     // 101010
    fmt::print("Hex: {:x}\n", 255);       // ff
    fmt::print("Hex upper: {:X}\n", 255); // FF
    fmt::print("Octal: {:o}\n", 64);      // 100
    fmt::print("Padded: {:05d}\n", 42);   // 00042
    fmt::print("Signed: {:+d}\n", 42);    // +42

#if 0
    // 6. Float formatting
    fmt::print("Default: {}\n", 3.14159);
    fmt::print("Fixed: {:.2f}\n", 3.14159);     // 3.14
    fmt::print("Scientific: {:.2e}\n", 1234.5); // 1.23e+03
    fmt::print("General: {:.3g}\n", 1234.5);    // 1.23e+03
    fmt::print("Percentage: {:.1%}\n", 0.75);   // 75.0%

    // 7. Alignment and width
    fmt::print("Left: '{:<10}'\n", "text");   // 'text      '
    fmt::print("Right: '{:>10}'\n", "text");  // '      text'
    fmt::print("Center: '{:^10}'\n", "text"); // '   text   '
    fmt::print("Fill: '{:*^10}'\n", "text");  // '***text***'

    // 8. Containers (requires fmt/ranges.h)
    std::vector<int> vec = {1, 2, 3, 4, 5};
    fmt::print("Vector: {}\n", vec);      // [1, 2, 3, 4, 5]
    fmt::print("Custom: {::02d}\n", vec); // [01, 02, 03, 04, 05]

    std::map<std::string, int> scores = {{"Alice", 100}, {"Bob", 85}};
    fmt::print("Map: {}\n", scores);

    // 9. Date/time formatting (requires fmt/chrono.h)
    auto now = std::chrono::system_clock::now();
    fmt::print("Current time: {}\n", now);
    fmt::print("Custom time: {:%Y-%m-%d %H:%M:%S}\n", now);

    // 10. Colored output (requires fmt/color.h)
    fmt::print(fg(fmt::color::red), "Red text\n");
    fmt::print(fg(fmt::color::green) | bg(fmt::color::yellow),
               "Green on yellow\n");
    fmt::print(fmt::emphasis::bold | fg(fmt::color::blue), "Bold blue\n");
    // 11. Error handling with exceptions
    try {
        // This will throw because there's no argument for {}
        // fmt::format("Missing argument: {}");
    } catch (const fmt::format_error &e) {
        fmt::print("Format error: {}\n", e.what());
    }

    // 12. Compile-time format string checking (C++20)
    // fmt::print(FMT_STRING("Checked: {} {}\n"), 42, "hello");

    // 13. Custom types formatting
    // You can define custom formatter (shown separately below)

    // 14. Locale-aware formatting
    fmt::print("Locale number: {:L}\n", 1234567);

    // 15. Memory buffer (no heap allocation for small strings)
    fmt::memory_buffer buf;
    fmt::format_to(std::back_inserter(buf), "Buffer: {} {}", 42, "test");
    fmt::print("From buffer: {}\n", fmt::to_string(buf));

    // 16. String views and efficient formatting
    std::string_view sv = "view";
    fmt::print("String view: {}\n", sv);

    // 17. Conditional formatting
    bool condition = true;
    fmt::print("Status: {}\n", condition ? "OK" : "FAIL");

    // 18. Multiple format specifiers
    fmt::print("Complex: {:#08x} {:.2f} {:>10}\n", 255, 3.14, "right");
#endif
    return 0;
}

// Example usage of custom formatter:
// Point p{10, 20};
// fmt::print("Custom point: {}\n", p);  // Output: Custom point: Point(10, 20)