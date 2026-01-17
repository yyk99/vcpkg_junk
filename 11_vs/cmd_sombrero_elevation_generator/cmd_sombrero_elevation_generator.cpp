//
// Sombrero image generator
//

#if _WIN32
#define _CRT_SECURE_NO_WARNINGS
#endif

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <fstream>
#include <iomanip>

class SombreroElevationGrid {
private:
    int width;
    int height;
    double xMin, xMax, yMin, yMax;
    std::vector<std::vector<double>> elevationData;

public:
    SombreroElevationGrid(int w, int h, double xmin = -10.0, double xmax = 10.0,
                          double ymin = -10.0, double ymax = 10.0)
        : width(w), height(h), xMin(xmin), xMax(xmax), yMin(ymin), yMax(ymax) {
        elevationData.resize(height, std::vector<double>(width));
    }

    // Sombrero function: z = sin(sqrt(x² + y²)) / sqrt(x² + y²)
    double sombrero(double x, double y) {
        double r = sqrt(x * x + y * y);
        if (r < 1e-10) {
            return 1.0; // Limit as r approaches 0
        }
        return sin(r) / r;
    }

    // Generate the elevation grid
    void generateGrid() {
        std::cout << "Generating Sombrero elevation grid (" << width << "x"
                  << height << ")..." << std::endl;

        for (int row = 0; row < height; ++row) {
            for (int col = 0; col < width; ++col) {
                // Convert pixel coordinates to world coordinates
                double x = xMin + (xMax - xMin) * col / (width - 1);
                double y = yMin + (yMax - yMin) * row / (height - 1);

                elevationData[row][col] = sombrero(x, y);
            }
        }

        std::cout << "Grid generation completed." << std::endl;
    }

    // Save as grayscale PNG (elevation as intensity)
    bool saveAsGrayscalePNG(const std::string &filename) {
        std::vector<unsigned char> imageData(width * height);

        // Find min and max elevations for normalization
        double minElev = elevationData[0][0];
        double maxElev = elevationData[0][0];

        for (int row = 0; row < height; ++row) {
            for (int col = 0; col < width; ++col) {
                minElev = std::min(minElev, elevationData[row][col]);
                maxElev = std::max(maxElev, elevationData[row][col]);
            }
        }

        std::cout << "Elevation range: " << minElev << " to " << maxElev
                  << std::endl;

        // Normalize and convert to 8-bit grayscale
        for (int row = 0; row < height; ++row) {
            for (int col = 0; col < width; ++col) {
                double normalized =
                    (elevationData[row][col] - minElev) / (maxElev - minElev);
                imageData[row * width + col] =
                    static_cast<unsigned char>(normalized * 255);
            }
        }

        int result = stbi_write_png(filename.c_str(), width, height, 1,
                                    imageData.data(), width);

        if (result) {
            std::cout << "✓ Grayscale PNG saved: " << filename << std::endl;
        } else {
            std::cout << "✗ Failed to save grayscale PNG: " << filename
                      << std::endl;
        }

        return result != 0;
    }

    // Save as color-mapped PNG (using a color scheme for elevation)
    bool saveAsColorPNG(std::string const &filename) {
        std::vector<unsigned char> imageData(width * height * 3); // RGB

        // Find min and max elevations for normalization
        double minElev = elevationData[0][0];
        double maxElev = elevationData[0][0];

        for (int row = 0; row < height; ++row) {
            for (int col = 0; col < width; ++col) {
                minElev = std::min(minElev, elevationData[row][col]);
                maxElev = std::max(maxElev, elevationData[row][col]);
            }
        }

        // Convert to RGB using a color map
        for (int row = 0; row < height; ++row) {
            for (int col = 0; col < width; ++col) {
                double normalized =
                    (elevationData[row][col] - minElev) / (maxElev - minElev);

                // Color mapping: blue (low) -> green -> yellow -> red (high)
                unsigned char r, g, b;
                elevationToColor(normalized, r, g, b);

                int pixelIndex = (row * width + col) * 3;
                imageData[pixelIndex] = r;     // Red
                imageData[pixelIndex + 1] = g; // Green
                imageData[pixelIndex + 2] = b; // Blue
            }
        }

        int result = stbi_write_png(filename.c_str(), width, height, 3,
                                    imageData.data(), width * 3);

        if (result) {
            std::cout << "✓ Color PNG saved: " << filename << std::endl;
        } else {
            std::cout << "✗ Failed to save color PNG: " << filename
                      << std::endl;
        }

        return result != 0;
    }

    // Save as heightmap PNG (16-bit grayscale for higher precision)
    bool saveAsHeightmapPNG(std::string const &filename) {
        std::vector<unsigned short> imageData(width * height);

        // Find min and max elevations
        double minElev = elevationData[0][0];
        double maxElev = elevationData[0][0];

        for (int row = 0; row < height; ++row) {
            for (int col = 0; col < width; ++col) {
                minElev = std::min(minElev, elevationData[row][col]);
                maxElev = std::max(maxElev, elevationData[row][col]);
            }
        }

        // Normalize to 16-bit range
        for (int row = 0; row < height; ++row) {
            for (int col = 0; col < width; ++col) {
                double normalized =
                    (elevationData[row][col] - minElev) / (maxElev - minElev);
                imageData[row * width + col] =
                    static_cast<unsigned short>(normalized * 65535);
            }
        }

        // Convert to 8-bit for PNG (stb_image_write doesn't support 16-bit PNG
        // directly)
        std::vector<unsigned char> imageData8bit(width * height);
        for (size_t i = 0; i < imageData.size(); ++i) {
            imageData8bit[i] =
                static_cast<unsigned char>(imageData[i] >> 8); // Take high byte
        }

        int result = stbi_write_png(filename.c_str(), width, height, 1,
                                    imageData8bit.data(), width);

        if (result) {
            std::cout << "✓ Heightmap PNG saved: " << filename << std::endl;
        } else {
            std::cout << "✗ Failed to save heightmap PNG: " << filename
                      << std::endl;
        }

        return result != 0;
    }

    /// @brief Exports elevation data in standard CSV format with X,Y,Elevation columns
    /// @param filename 
    /// @return true if success
    bool saveAsCSV(std::string const &filename) {
        std::ofstream csvFile(filename);

        if (!csvFile.is_open()) {
            std::cout << "✗ Failed to open CSV file: " << filename << std::endl;
            return false;
        }

        std::cout << "Saving CSV file: " << filename << "..." << std::endl;

        // Write header with coordinate information as comments
        csvFile << "# Sombrero Elevation Grid CSV Export\n";
        csvFile << "# Grid size: " << width << "x" << height << "\n";
        csvFile << "# Domain: x ∈ [" << xMin << ", " << xMax << "], y ∈ ["
                << yMin << ", " << yMax << "]\n";
        csvFile << "# Format: X,Y,Elevation\n";

        // Write data points
        for (int row = 0; row < height; ++row) {
            for (int col = 0; col < width; ++col) {
                // Convert pixel coordinates to world coordinates
                double x = xMin + (xMax - xMin) * col / (width - 1);
                double y = yMin + (yMax - yMin) * row / (height - 1);
                double elevation = elevationData[row][col];

                csvFile << std::fixed << std::setprecision(6) << x << "," << y
                        << "," << elevation << "\n";
            }
        }

        csvFile.close();

        if (csvFile.good()) {
            std::cout << "✓ CSV saved: " << filename << std::endl;
            return true;
        } else {
            std::cout << "✗ Error occurred while writing CSV: " << filename
                      << std::endl;
            return false;
        }
    }

private:
    // Convert elevation value (0.0 to 1.0) to RGB color
    void elevationToColor(double value, unsigned char &r, unsigned char &g,
                          unsigned char &b) {
        // Clamp value to [0, 1]
        value = std::max(0.0, std::min(1.0, value));

        if (value < 0.25) {
            // Blue to Cyan
            double t = value / 0.25;
            r = 0;
            g = static_cast<unsigned char>(t * 255);
            b = 255;
        } else if (value < 0.5) {
            // Cyan to Green
            double t = (value - 0.25) / 0.25;
            r = 0;
            g = 255;
            b = static_cast<unsigned char>((1.0 - t) * 255);
        } else if (value < 0.75) {
            // Green to Yellow
            double t = (value - 0.5) / 0.25;
            r = static_cast<unsigned char>(t * 255);
            g = 255;
            b = 0;
        } else {
            // Yellow to Red
            double t = (value - 0.75) / 0.25;
            r = 255;
            g = static_cast<unsigned char>((1.0 - t) * 255);
            b = 0;
        }
    }


public:
    // Print some statistics about the elevation data
    void printStatistics() {
        double minElev = elevationData[0][0];
        double maxElev = elevationData[0][0];
        double sum = 0.0;

        for (int row = 0; row < height; ++row) {
            for (int col = 0; col < width; ++col) {
                double elev = elevationData[row][col];
                minElev = std::min(minElev, elev);
                maxElev = std::max(maxElev, elev);
                sum += elev;
            }
        }

        double mean = sum / (width * height);

        std::cout << "\n=== Elevation Statistics ===" << std::endl;
        std::cout << "Grid size: " << width << " x " << height << std::endl;
        std::cout << "Domain: x ∈ [" << xMin << ", " << xMax << "], y ∈ ["
                  << yMin << ", " << yMax << "]" << std::endl;
        std::cout << "Min elevation: " << minElev << std::endl;
        std::cout << "Max elevation: " << maxElev << std::endl;
        std::cout << "Mean elevation: " << mean << std::endl;
        std::cout << "Range: " << (maxElev - minElev) << std::endl;
    }
};

#if _WIN32
#include <windows.h>
#endif

int main() {
#if _WIN32
    SetConsoleOutputCP(65001);
#endif
    std::cout << "=== Sombrero Function Elevation Grid Generator ==="
              << std::endl
              << std::endl;

    // Create different resolutions and domains for testing
    std::vector<std::pair<int, int>> configurations = {
        {512, 512},
        {1024, 1024},
        {256, 256},
    };

    for (const auto &config : configurations) {
        int width = std::get<0>(config);
        int height = std::get<1>(config);
        std::string baseName =
            "sombrero_" + std::to_string(width) + "x" + std::to_string(height);

        std::cout << "\n--- Generating " << baseName << " ---" << std::endl;

        // For the close-up version, use a smaller domain
        double domain =
            (baseName.find("closeup") != std::string::npos) ? 5.0 : 10.0;

        SombreroElevationGrid grid(width, height, -domain, domain, -domain,
                                   domain);

        // Generate the elevation data
        grid.generateGrid();

        // Print statistics
        grid.printStatistics();

        // Save in different formats

        std::filesystem::create_directories("out");

        grid.saveAsGrayscalePNG("out/" + baseName + "_grayscale.png");
        grid.saveAsColorPNG("out/" + baseName + "_color.png");
        grid.saveAsHeightmapPNG("out/" + baseName + "_heightmap.png");
        grid.saveAsCSV("out/" + baseName + "_elevation.csv");
    }

    return 0;
}