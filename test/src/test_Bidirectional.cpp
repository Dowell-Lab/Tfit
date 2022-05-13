/**
 * @file test_Bidirectional.cpp 
 * @author Robin Dowell
 * @brief Testing the bidirectional model class
 * @date 2022-05-13
 * 
 */
#include "gmock/gmock.h"
#include "Model.h"

TEST(Models, bidir_NormPDF)
{
    // Arrange: bring SUT to desired state
    Bidirectional sut(35., 0.5, 0.25, 0.5, 2.);

    // Act: call methods on SUT, capture output
    double result = sut.normalPDF(0.3);

    // Assert: Verify the outcome
    EXPECT_EQ(sut.write_out(), "Bidir(35.00,0.50,0.2500,0.500,2.00)");
    ASSERT_LE(abs(result - 0.381388), 0.0001);  // From R's dnorm(x=0.3)

}