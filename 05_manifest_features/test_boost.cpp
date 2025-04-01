//
//
//

#include <gtest/gtest.h>
#include <iostream>
#include <boost/filesystem.hpp>

TEST(t1, t1)
{
    std::cout << "Hello...\n";
}

TEST(t1, fs1)
{
    boost::filesystem::path actual(__FILE__);

    EXPECT_TRUE(boost::filesystem::exists(actual));
}

#include <CGAL/Simple_cartesian.h>
#include <iostream>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Segment_2 Segment_2;

TEST(CGAL, t0) {
    Point_2 p(1, 1), q(10, 10);
    std::cout << "p = " << p << std::endl;
    std::cout << "q = " << q.x() << " " << q.y() << std::endl;
    std::cout << "sqdist(p,q) = " << CGAL::squared_distance(p, q) << std::endl;
    Segment_2 s(p, q);
    Point_2 m(5, 9);
    std::cout << "m = " << m << std::endl;
    std::cout << "sqdist(Segment_2(p,q), m) = " << CGAL::squared_distance(s, m)
              << std::endl;
    std::cout << "p, q, and m ";
    switch (CGAL::orientation(p, q, m)) {
    case CGAL::COLLINEAR:
        std::cout << "are collinear\n";
        break;
    case CGAL::LEFT_TURN:
        std::cout << "make a left turn\n";
        break;
    case CGAL::RIGHT_TURN:
        std::cout << "make a right turn\n";
        break;
    }
    std::cout << " midpoint(p,q) = " << CGAL::midpoint(p, q) << std::endl;
}

#include <oneapi/tbb/info.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/task_arena.h>

TEST(TBB, t0) {
    // Get the default number of threads
    int num_threads = oneapi::tbb::info::default_concurrency();

    std::cout << "num_threads: " << num_threads << "\n";

    // Run the default parallelism
}