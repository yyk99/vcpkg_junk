//
// Requested
// example to build an octree with PCL library?
//

#define _CRT_SECURE_NO_WARNINGS 1

#ifdef _WIN32
#   pragma warning (disable:4805)
#endif

#define CHANGEDETECTIONEXAMPLE 0 /* TODO: fix the compiler */
#define VOXELCENTROIDEXAMPLE 0   /* TODO: fix it */

#include <ctime>
#include <iostream>
#include <pcl/common/time.h>
#include <pcl/io/pcd_io.h>
#include <pcl/octree/octree_pointcloud_changedetector.h>
#include <pcl/octree/octree_pointcloud_occupancy.h>
#include <pcl/octree/octree_pointcloud_voxelcentroid.h>
#include <pcl/octree/octree_search.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <random>
#include <vector>

// Type definitions for cleaner code
typedef pcl::PointXYZ PointT;
typedef pcl::PointCloud<PointT> PointCloudT;

class PCLOctreeExamples {
public:
    // Generate a random point cloud for testing
    static PointCloudT::Ptr generateRandomPointCloud(int numPoints = 1000,
                                                     float minCoord = -10.0f,
                                                     float maxCoord = 10.0f) {
        PointCloudT::Ptr cloud(new PointCloudT);
        cloud->width = numPoints;
        cloud->height = 1;
        cloud->points.resize(cloud->width * cloud->height);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<float> dis(minCoord, maxCoord);

        for (size_t i = 0; i < cloud->points.size(); ++i) {
            cloud->points[i].x = dis(gen);
            cloud->points[i].y = dis(gen);
            cloud->points[i].z = dis(gen);
        }

        return cloud;
    }

    // Example 1: Basic Octree Search
    static void basicOctreeSearch() {
        std::cout << "\n=== Basic Octree Search Example ===" << std::endl;

        // Generate point cloud
        PointCloudT::Ptr cloud = generateRandomPointCloud(1000);
        std::cout << "Generated point cloud with " << cloud->points.size()
                  << " points" << std::endl;

        // Create octree with resolution 1.0f
        float resolution = 1.0f;
        pcl::octree::OctreePointCloudSearch<PointT> octree(resolution);

        // Input point cloud to octree
        octree.setInputCloud(cloud);
        octree.addPointsFromInputCloud();

        std::cout << "Octree created with resolution: " << resolution
                  << std::endl;
        std::cout << "Number of leaf nodes: " << octree.getLeafCount()
                  << std::endl;

        // Example search point
        PointT searchPoint(0.0f, 0.0f, 0.0f);

        // 1. Voxel search (find points in same voxel)
        std::vector<int> pointIdxVec;
        if (octree.voxelSearch(searchPoint, pointIdxVec)) {
            std::cout << "\nVoxel search found " << pointIdxVec.size()
                      << " points in same voxel:" << std::endl;
            for (size_t i = 0; i < std::min(size_t(5), pointIdxVec.size());
                 ++i) {
                PointT &point = cloud->points[pointIdxVec[i]];
                std::cout << "  Point " << pointIdxVec[i] << ": (" << point.x
                          << ", " << point.y << ", " << point.z << ")"
                          << std::endl;
            }
        } else {
            std::cout << "\nNo points found in voxel containing search point"
                      << std::endl;
        }

        // 2. K-nearest neighbor search
        int K = 10;
        std::vector<int> pointIdxNKNSearch;
        std::vector<float> pointNKNSquaredDistance;

        std::cout << "\nK-nearest neighbor search (K=" << K
                  << "):" << std::endl;
        if (octree.nearestKSearch(searchPoint, K, pointIdxNKNSearch,
                                  pointNKNSquaredDistance) > 0) {
            for (size_t i = 0; i < pointIdxNKNSearch.size(); ++i) {
                PointT &point = cloud->points[pointIdxNKNSearch[i]];
                std::cout << "  " << i + 1 << ". Point " << pointIdxNKNSearch[i]
                          << ": (" << point.x << ", " << point.y << ", "
                          << point.z
                          << ") distance: " << sqrt(pointNKNSquaredDistance[i])
                          << std::endl;
            }
        }

        // 3. Radius search
        float radius = 2.0f;
        std::vector<int> pointIdxRadiusSearch;
        std::vector<float> pointRadiusSquaredDistance;

        std::cout << "\nRadius search (radius=" << radius << "):" << std::endl;
        if (octree.radiusSearch(searchPoint, radius, pointIdxRadiusSearch,
                                pointRadiusSquaredDistance) > 0) {
            std::cout << "Found " << pointIdxRadiusSearch.size()
                      << " points within radius:" << std::endl;
            for (size_t i = 0;
                 i < std::min(size_t(10), pointIdxRadiusSearch.size()); ++i) {
                PointT &point = cloud->points[pointIdxRadiusSearch[i]];
                std::cout << "  Point " << pointIdxRadiusSearch[i] << ": ("
                          << point.x << ", " << point.y << ", " << point.z
                          << ") distance: "
                          << sqrt(pointRadiusSquaredDistance[i]) << std::endl;
            }
        }
    }
#if CHANGEDETECTIONEXAMPLE
    // Example 2: Change Detection
    static void changeDetectionExample() {
        std::cout << "\n=== Change Detection Example ===" << std::endl;

        // Create two point clouds
        PointCloudT::Ptr cloudA = generateRandomPointCloud(500, -5.0f, 5.0f);
        PointCloudT::Ptr cloudB = generateRandomPointCloud(500, -5.0f, 5.0f);

        // Add some common points
        for (int i = 0; i < 100; ++i) {
            cloudB->points[i] = cloudA->points[i];
        }

        std::cout << "Cloud A: " << cloudA->points.size() << " points"
                  << std::endl;
        std::cout << "Cloud B: " << cloudB->points.size() << " points"
                  << std::endl;

        float resolution = 0.5f;
        pcl::octree::OctreePointCloudChangeDetector<PointT> octree(resolution);

        // Set input cloud and build octree structure
        octree.setInputCloud(cloudA);
        octree.addPointsFromInputCloud();

        // Switch buffers and set new input cloud
        octree.switchBuffers();
        octree.setInputCloud(cloudB);
        octree.addPointsFromInputCloud();

        // Get vector of point indices from new cloud that didn't exist in
        // previous buffer
        std::vector<int> newPointIdxVector;
        octree.getPointIndicesFromNewVoxels(newPointIdxVector);

        std::cout << "Change detection found " << newPointIdxVector.size()
                  << " new/changed points:" << std::endl;

        for (size_t i = 0; i < std::min(size_t(10), newPointIdxVector.size());
             ++i) {
            PointT &point = cloudB->points[newPointIdxVector[i]];
            std::cout << "  New point " << newPointIdxVector[i] << ": ("
                      << point.x << ", " << point.y << ", " << point.z << ")"
                      << std::endl;
        }
    }
#endif
#if VOXELCENTROIDEXAMPLE
    // Example 3: Voxel Centroid Calculation
    static void voxelCentroidExample() {
        std::cout << "\n=== Voxel Centroid Example ===" << std::endl;

        // Generate dense point cloud
        PointCloudT::Ptr cloud = generateRandomPointCloud(2000);

        float resolution = 1.0f;
        pcl::octree::OctreePointCloudVoxelCentroid<PointT> octree(resolution);

        octree.setInputCloud(cloud);
        octree.addPointsFromInputCloud();

        std::cout << "Original cloud: " << cloud->points.size() << " points"
                  << std::endl;
        std::cout << "Octree leaf nodes: " << octree.getLeafCount()
                  << std::endl;

        // Get voxel centroids
        PointCloudT::Ptr voxelCentroids(new PointCloudT);
        octree.getVoxelCentroids(*voxelCentroids);

        std::cout << "Voxel centroids: " << voxelCentroids->points.size()
                  << " points" << std::endl;
        std::cout << "Compression ratio: "
                  << float(cloud->points.size()) / voxelCentroids->points.size()
                  << ":1" << std::endl;

        // Display first few centroids
        std::cout << "\nFirst 10 voxel centroids:" << std::endl;
        for (size_t i = 0;
             i < std::min(size_t(10), voxelCentroids->points.size()); ++i) {
            PointT &point = voxelCentroids->points[i];
            std::cout << "  Centroid " << i << ": (" << point.x << ", "
                      << point.y << ", " << point.z << ")" << std::endl;
        }
    }
#endif
    // Example 4: Occupancy Grid
    static void occupancyGridExample() {
        std::cout << "\n=== Occupancy Grid Example ===" << std::endl;

        PointCloudT::Ptr cloud = generateRandomPointCloud(1000);

        float resolution = 1.0f;
        pcl::octree::OctreePointCloudOccupancy<PointT> octree(resolution);

        octree.setInputCloud(cloud);
        octree.addPointsFromInputCloud();

        std::cout << "Occupancy octree created with " << octree.getLeafCount()
                  << " occupied voxels" << std::endl;

        // Check occupancy of specific voxels
        std::vector<PointT> testPoints = {
            PointT(0.0f, 0.0f, 0.0f), PointT(5.0f, 5.0f, 5.0f),
            PointT(-3.0f, 2.0f, 1.0f),
            PointT(15.0f, 15.0f, 15.0f) // Outside range
        };

        std::cout << "\nVoxel occupancy check:" << std::endl;
        for (const PointT &point : testPoints) {
            bool occupied = octree.isVoxelOccupiedAtPoint(point);
            std::cout << "  Point (" << point.x << ", " << point.y << ", "
                      << point.z << ") -> " << (occupied ? "OCCUPIED" : "EMPTY")
                      << std::endl;
        }

        // Get all occupied voxel centers
        std::vector<PointT, Eigen::aligned_allocator<PointT>> voxelCenters;
        octree.getOccupiedVoxelCenters(voxelCenters);

        std::cout << "\nTotal occupied voxels: " << voxelCenters.size()
                  << std::endl;
        std::cout << "First 5 occupied voxel centers:" << std::endl;
        for (size_t i = 0; i < std::min(size_t(5), voxelCenters.size()); ++i) {
            std::cout << "  Voxel " << i << ": (" << voxelCenters[i].x << ", "
                      << voxelCenters[i].y << ", " << voxelCenters[i].z << ")"
                      << std::endl;
        }
    }

    // Example 5: Performance Comparison
    static void performanceComparison() {
        std::cout << "\n=== Performance Comparison ===" << std::endl;

        // Generate larger point cloud for performance testing
        PointCloudT::Ptr cloud = generateRandomPointCloud(10000);
        PointT searchPoint(0.0f, 0.0f, 0.0f);
        float radius = 2.0f;

        std::cout << "Testing with " << cloud->points.size() << " points"
                  << std::endl;
        std::cout << "Search radius: " << radius << std::endl;

        // Octree search
        pcl::octree::OctreePointCloudSearch<PointT> octree(1.0f);
        octree.setInputCloud(cloud);

        auto start = pcl::getTime();
        octree.addPointsFromInputCloud();
        auto octree_build_time = pcl::getTime() - start;

        std::vector<int> pointIdxRadiusSearch;
        std::vector<float> pointRadiusSquaredDistance;

        start = pcl::getTime();
        int found =
            octree.radiusSearch(searchPoint, radius, pointIdxRadiusSearch,
                                pointRadiusSquaredDistance);
        auto octree_search_time = pcl::getTime() - start;

        // Brute force search for comparison
        std::vector<int> bruteForceResults;
        start = pcl::getTime();

        for (size_t i = 0; i < cloud->points.size(); ++i) {
            float dx = cloud->points[i].x - searchPoint.x;
            float dy = cloud->points[i].y - searchPoint.y;
            float dz = cloud->points[i].z - searchPoint.z;
            float distance = sqrt(dx * dx + dy * dy + dz * dz);

            if (distance <= radius) {
                bruteForceResults.push_back(i);
            }
        }
        auto brute_force_time = pcl::getTime() - start;

        std::cout << "\nResults:" << std::endl;
        std::cout << "Octree build time: " << octree_build_time * 1000 << " ms"
                  << std::endl;
        std::cout << "Octree search time: " << octree_search_time * 1000
                  << " ms" << std::endl;
        std::cout << "Brute force time: " << brute_force_time * 1000 << " ms"
                  << std::endl;
        std::cout << "Octree found: " << found << " points" << std::endl;
        std::cout << "Brute force found: " << bruteForceResults.size()
                  << " points" << std::endl;
        std::cout << "Speedup factor: " << brute_force_time / octree_search_time
                  << "x" << std::endl;
    }
};

int main() {
    std::cout << "PCL Octree Examples" << std::endl;
    std::cout << "==================" << std::endl;

    try {
        // Run all examples
        PCLOctreeExamples::basicOctreeSearch();
#if CHANGEDETECTIONEXAMPLE
        PCLOctreeExamples::changeDetectionExample();
#endif
#if VOXELCENTROIDEXAMPLE
        PCLOctreeExamples::voxelCentroidExample();
#endif
        PCLOctreeExamples::occupancyGridExample();
        PCLOctreeExamples::performanceComparison();

        std::cout << "\n=== Additional Octree Operations ===" << std::endl;

        // Demonstrate additional octree capabilities
        PointCloudT::Ptr cloud =
            PCLOctreeExamples::generateRandomPointCloud(1000);
        pcl::octree::OctreePointCloudSearch<PointT> octree(1.0f);
        octree.setInputCloud(cloud);
        octree.addPointsFromInputCloud();

        std::cout << "Octree depth: " << octree.getTreeDepth() << std::endl;
        std::cout << "Octree resolution: " << octree.getResolution()
                  << std::endl;
        std::cout << "Leaf nodes: " << octree.getLeafCount() << std::endl;
        std::cout << "Branch nodes: " << octree.getBranchCount() << std::endl;

        // Get bounding box
        Eigen::Vector3d min_pt, max_pt;
        octree.getBoundingBox(min_pt[0], min_pt[1], min_pt[2], max_pt[0],
                              max_pt[1], max_pt[2]);
        std::cout << "Bounding box: min(" << min_pt.transpose() << ") max("
                  << max_pt.transpose() << ")" << std::endl;

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

// CMakeLists.txt content for building:
/*
cmake_minimum_required(VERSION 3.10)
project(pcl_octree_example)

find_package(PCL 1.8 REQUIRED)

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_executable(pcl_octree_example pcl_octree_example.cpp)
target_link_libraries(pcl_octree_example ${PCL_LIBRARIES})

set_property(TARGET pcl_octree_example PROPERTY CXX_STANDARD 14)
*/