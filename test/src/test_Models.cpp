/**
 * @file test_Models.cpp 
 * @author Robin Dowell
 * @brief Testing the bidirectional model class
 * @date 2022-05-13
 * 
 */
#include "gmock/gmock.h"
#include <fstream>
#include "Models.h"

/*
TEST (basicModel, updateParameters) 
{

}

TEST (basicModel, calculateRi) 
{

}

TEST (basicModel, updateExpectations) 
{

}
*/

TEST(Bidirectional, applyFootprintSubPosStrand) 
{
  // Arrange: bring SUT to desired state
  Bidirectional sut(35., 1.5, 0.25, 0.5, 30.);
  // Act
  double newZ = sut.applyFootprint(33, '+');
  // Assess
  EXPECT_EQ(newZ, 3);
}
TEST(Bidirectional, applyFootprintAddNegStrand) 
{
  // Arrange: bring SUT to desired state
  Bidirectional sut(35., 1.5, 0.25, 0.5, 30.);
  // Act
  double newZ = sut.applyFootprint(33, '-');
  // Assess
  EXPECT_EQ(newZ, 63);
}


TEST(Bidirectional, MillsRatioCorrect)
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

TEST(Bidirectional, generateDataHasCorrectMean)
{
    // Arrange: bring SUT to desired state
    Bidirectional sut(35., 0.5, 0.25, 1.0, 2.);

    // Act: call methods on SUT, capture output
    std::vector<double> gdata = sut.generate_data(100000);

    double total = 0;
    std::vector<double>::iterator it;
    for (it = gdata.begin(); it != gdata.end(); it++) {
      total += (*it);
    }
    double mean = total / gdata.size();

    // Assert: Verify the outcome
    ASSERT_LE(mean - (sut.loading.mu + (1/sut.initiation.lambda)), 0.01);  
}

TEST(Bidirectional, GrushkaEQKalambetPositive)
{
    // Arrange: bring SUT to desired state
    Bidirectional sut(35., 1.5, 0.25, 0.5, 0.);

    // Act: call methods on SUT, capture output
    double grushka = sut.pdf(37.,'+');
    double kalambet = sut.pdf_alt(37.,'+');

    // Assert: Verify the outcome
    ASSERT_LE(abs(grushka - kalambet), 0.0001);  
}

TEST(Bidirectional, GrushkaEQKalambetNegative)
{
    // Arrange: bring SUT to desired state
    Bidirectional sut(35., 1.5, 0.25, 0.5, 0.);

    // Act: call methods on SUT, capture output
    double grushka = sut.pdf(37.,'-');
    double kalambet = sut.pdf_alt(37.,'-');

    // Assert: Verify the outcome
    ASSERT_LE(abs(grushka - kalambet), 0.0001);  
}

// double ExpX(double z, char strand);
TEST(Bidirectional, ExpX)
{
    // Arrange: bring SUT to desired state
    Bidirectional sut(35., 1.5, 0.25, 0.5, 0.);

    // Act: call methods on SUT, capture output
    double expectedValX = sut.ExpX(33, '+');

    // Assert: Verify the outcome
    ASSERT_LE(abs(expectedValX - 32.3861), 0.0001);  
}

// double ExpY(double z, char s);
TEST(Bidirectional, ExpY) 
{
  // Arrange: bring SUT to desired state
  Bidirectional sut(35., 1.5, 0.25, 0.5, 0.);

  // Act: call methods on SUT, capture output
  double expectedValY = sut.ExpY(33, '+');

  // Assert: Verify the outcome
  ASSERT_LE(abs(expectedValY - 0.613866), 0.0001);  
}

// double ExpX2(double z, char strand);
TEST(Bidirectional, ExpX2) 
{
  // Arrange: bring SUT to desired state
  Bidirectional sut(35., 1.5, 0.25, 0.5, 0.);

  // Act: call methods on SUT, capture output
  double X2expected = sut.ExpX2(33, '+');

  // Assert: Verify the outcome
  ASSERT_LE(abs(X2expected - 1048.9248), 0.0001);  
}

// double ExpY2(double z, char s);
TEST(Bidirectional, ExpY2) 
{
  // Arrange: bring SUT to desired state
  Bidirectional sut(35., 1.5, 0.25, 0.5, 0.);

  // Act: call methods on SUT, capture output
  double Y2expected = sut.ExpY2(33, '+');

  // Assert: Verify the outcome
  ASSERT_LE(abs(Y2expected - 0.676969), 0.0001);  
}

// Models, bidir_setPriors
// Models, bidir_calculate EXP
// Models, bidir_updateExpectations
// Models, bidir_calcExpectedVals
// Models, bidir_updateParameters


// Models, uniform_pdf
// Models, uniform_updateParameters
// Models, uniform_setPi
// Models, uniform_calculateLiklihood
// Models, uniform_setBounds


// Models, fullModel_pdf
// Models, fullModel_resetSufficiency
// Models, fullModel_getResponsibility
// Models, fullModel_updateParameters
// Models, fullModel_calculateRi
// Models, fullModel_updateExpectations
