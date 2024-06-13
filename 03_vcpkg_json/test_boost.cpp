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