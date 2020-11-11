#include <gtest/gtest.h>
#include <math.h>
 
double squareRoot(const double a) {
    double b = sqrt(a);
    if(b != b) { // NaN check
        return -1.0;
    }
    else{
        return sqrt(a);
    }
}

// https://www.eriksmistad.no/getting-started-with-google-test-on-ubuntu/

TEST(SquareRootTest, Positive) { 
    ASSERT_EQ(squareRoot(36.0), 6.);
    ASSERT_EQ(squareRoot(324.0), 18.);
    ASSERT_EQ(squareRoot(645.16), 25.4);
    ASSERT_EQ(squareRoot(0.0), 0.);
}
 
TEST(SquareRootTest, Negative) {
    ASSERT_EQ(squareRoot(-15.0), -1.);
    ASSERT_EQ(squareRoot(-0.2), -1.);
}
 
int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
