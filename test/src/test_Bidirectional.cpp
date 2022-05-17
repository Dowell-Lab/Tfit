/**
 * @file test_Bidirectional.cpp 
 * @author Robin Dowell
 * @brief Testing the bidirectional model class
 * @date 2022-05-13
 * 
 */
#include "gmock/gmock.h"
#include <fstream>
#include "Model.h"

TEST(Models, bidir_standardnorms)
{
    // Arrange: bring SUT to desired state
    Bidirectional sut(35., 0.5, 0.25, 0.5, 2.);

    // Act: call methods on SUT, capture output
    double result = sut.normalPDF(0.3);
    double result2 = sut.normalCDF(0.3);

    // Assert: Verify the outcome
    EXPECT_EQ(sut.write_out(), "Bidir(35.00,0.50,0.2500,0.500,2.00)");
    ASSERT_LE(abs(result - 0.381388), 0.0001);  // From R's dnorm(x=0.3)
    ASSERT_LE(abs(result2 - 0.617911), 0.0001);  // From R's pnorm(x=0.3)
}

TEST(Models, bidir_MillsRatio)
{
    // Arrange: bring SUT to desired state
    Bidirectional sut(35., 0.5, 0.25, 0.5, 2.);

    // Act: call methods on SUT, capture output
    double result = sut.millsRatio(11);
    double result2 = sut.millsRatio(2);

    // Assert: Verify the outcome
    ASSERT_LE(abs(result - (1.0/11)), 0.0001);  // Asymptotic behavior
    // in R: exp(pnorm(x, lower.tail=FALSE, log.p=TRUE) - dnorm(x, log=TRUE))
    ASSERT_LE(abs(result2 - 0.421369), 0.0001);  // .421369 per R
}

TEST(Models, bidir_generateData)
{
    // Arrange: bring SUT to desired state
    Bidirectional sut(35., 0.5, 0.25, 0.5, 2.);

    // Act: call methods on SUT, capture output
    std::vector<double> gdata = sut.generate_data(1000);

    /*
    std::string filename = "temp.out";
    std::ofstream OUT;
    OUT.open(filename);
    std::vector<double>::iterator it;
    for (it = gdata.begin(); it != gdata.end(); it++) {
       if (*it < 0) {
         OUT << "0\t" + std::to_string(*it) << std::endl;
       } else {
         OUT << std::to_string(*it) + "\t0" << std::endl;
       }
    }
    OUT.close();
    */

    // Assert: Verify the outcome
    EXPECT_EQ(gdata.size(), 1000);
    // Should test some statistical properties of generated data.
}

/*
TEST(Models, bidir_ExpectedValues)
{
    // Arrange: bring SUT to desired state
    Bidirectional sut(35., 0.5, 0.25, 0.5, 2.);

    // Act: call methods on SUT, capture output
    double result = sut.ExpY(40, '+');
    double result2 = sut.ExpY2(40, '+');
    double result = sut.ExpX(40, '+');
    double result2 = sut.ExpX2(40, '+');

    // Assert: Verify the outcome
}
*/