#include "../spin.hpp"
#include <gtest/gtest.h>



TEST(test_spin, Dotproduct) {
  
  Spin3D s1(0, 0);
    Spin3D s2(0, 1);
    Spin3D s3(std::numbers::pi, 0);

    Spin3D s4(std::numbers::pi/2, 0);
    Spin3D s5(std::numbers::pi/2, std::numbers::pi/2);

    EXPECT_EQ(s1*s2, 1);

    EXPECT_EQ(s1*s3, -1);

    EXPECT_NEAR(s1*s4, 0, 1e-6);

    EXPECT_NEAR(s1*s5, 0, 1e-6);

    EXPECT_NEAR(s4*s5, 0, 1e-6);




}

