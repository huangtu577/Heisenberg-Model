#include "../spin.hpp"
#include <gtest/gtest.h>



TEST(test_spin, Dotproduct) {
  
  Spin3D s1(0, 0);
    Spin3D s2(0, 1);
    Spin3D s3(std::numbers::pi, 0);

    Spin3D s4(std::numbers::pi/2, 0);
    Spin3D s5(std::numbers::pi/2, std::numbers::pi/2);

    EXPECT_EQ(s1*s2, 1);
    EXPECT_EQ(s1*s2, s2*s1);

    EXPECT_EQ(s1*s3, -1);
    EXPECT_EQ(s1*s3, s3*s1);


    EXPECT_NEAR(s1*s4, 0, 1e-6);
        EXPECT_EQ(s1*s4, s4*s1);


    EXPECT_NEAR(s1*s5, 0, 1e-6);
        EXPECT_EQ(s1*s5, s5*s1);



    EXPECT_EQ(s4*s5, s5*s4);

    EXPECT_EQ(s4*s3, s3*s4);

    EXPECT_EQ(s5*s3, s3*s5);

    EXPECT_EQ(s5*s2, s2*s5);

    EXPECT_EQ(s4*s2, s2*s4);

    EXPECT_EQ(s3*s2, s2*s3);



}

// TEST(test_spin, MAGNETIZATION){
//   size_t L{3};
//   double J{1.0};
//   double beta{0.5};
//   Grid3D grid(L, L, L, J, beta, 7, 1e5, "test.h5", false);

//   EXPECT_EQ(grid.magnetization(), 1);

// }
