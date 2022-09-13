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

TEST (BasicModel, updateParameters) 
{
  //Arrange
  BasicModel sut;
  //Act
  sut.updateParameters(10.,3.);
  //Assert
  EXPECT_EQ(sut.getPi(), 0.5);
  double expect = 1./28.;
  EXPECT_EQ(sut.getWeight(), expect);
}

/**
 * @brief Need to double check -- this is a dumb test.
 */
TEST (BasicModel, calculateRi) 
{
  //Arrange
  perStrandInfo point(3.,5.);
  perStrandInfo newRi;
  BasicModel sut;

  //Act
  newRi = sut.calculateRi(10, point);
  
  // Assert
  EXPECT_EQ(newRi.forward, 1.);
  EXPECT_EQ(newRi.reverse, 1.);
}

/**
 * @brief Need to double check -- this is a dumb test.
 */
TEST (BasicModel, updateExpectations) 
{
  //Arrange
  perStrandInfo point(3.,5.);
  perStrandInfo normalizedRi(15.,20.);
  BasicModel sut;

  //Act
  sut.updateExpectations(point, normalizedRi);
  
  // Assert
  EXPECT_EQ(sut.sufficiencyStats.Rk.forward, 0.);
  EXPECT_EQ(sut.sufficiencyStats.Rk.reverse, 0.);
}

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

/**
 * @brief Need to double check -- this is a dumb test.
 */
TEST(Bidirectional, calcExpectedVals)
{
  //Arrange
  Bidirectional sut;
  //Act
  sut.calcExpectedVals(6., '+', 36.);
  //Assert
  EXPECT_EQ(sut.sumOverN.sumRExpX, 216.);
  EXPECT_EQ(sut.sumOverN.sumRExpY, 0.);
  EXPECT_EQ(sut.sumOverN.sumRExpX2, 1332.);
  EXPECT_EQ(sut.sumOverN.sumRXExpY, 9288.);
}

/**
 * @brief Need to double check -- this is a dumb test.
 */
TEST(Bidirectional, updateExpectations) 
{
  //Arrange
  perStrandInfo point(3.,5.);
  perStrandInfo normalizedRi(15.,20.);
  Bidirectional sut;

  //Act
  sut.updateExpectations(6., point, normalizedRi);
  
  // Assert
  EXPECT_EQ(sut.sufficiencyStats.Rk.forward, 0.);
  EXPECT_EQ(sut.sufficiencyStats.Rk.reverse, 0.);
}

TEST(Bidirectional, updateParameters) 
{
  //Arrange
  Bidirectional sut;
  //Act
  sut.updateParameters(10.,3.);
  //Assert
  EXPECT_EQ(sut.getPi(), 0.5);
  double expect = 1./28.;
  EXPECT_EQ(sut.getWeight(), expect);
  EXPECT_EQ(sut.getMu(), 0.);
  EXPECT_LE(abs(sut.getSigma() - 0.707107), 0.01);
  EXPECT_EQ(sut.getLambda(), 1.);

  // Should footprint be tested?
}

TEST(UniformModel, pdfUnweighted)
{
  // Arrange
  UniformModel sut(5,10);
  sut.setWeight(1.0);
  //Act
  double ans = sut.pdf(6,'+');
  //Assert
  EXPECT_LE(abs(ans - 0.2), 0.0001);  // R: dunif(6,min=5,max=10);
}

TEST(UniformModel, pdfWeighted)
{
  // Arrange
  UniformModel sut(5,10);
  sut.setWeight(3.0);
  //Act
  double ans = sut.pdf(6,'+');
  //Assert
  EXPECT_LE(abs(ans - (3*0.2)), 0.0001);  // R: 3*dunif(6,min=5,max=10);
}

/**
 * @brief Need to double check -- this is a dumb test.
 */
TEST(UniformModel, calculateLikelihood)
{
  // Arrange: bring SUT to desired state
  RawData data;
  data.addDataPoints(10, 14, 4);   
  data.addDataPoints(14, 15, 1);
  data.addDataPoints(8, 10, -2);

  dInterval roi(&data,2,1);
  UniformModel sut;
  sut.initalizeBounds(5,20,1.0,1.0);

  // Act
  double ans = sut.calculateLikelihood(&roi);

  //Assert
  EXPECT_DOUBLE_EQ(ans, -33.08047253394033);
}

/**
 * @brief Need to double check -- this is a dumb test.
 */
TEST(FullModel, pdf)
{
  //Arrange
  FullModel sut;
  sut.initBounds(7., 2., 1., 0.33, 5., 30.);
  sut.bidir.setFootprint(2.);
  // std::cout << sut.write_out() << std::endl;

  //Act
  double ans = sut.pdf(3.,'+');
  //Assert
  EXPECT_DOUBLE_EQ(ans,0.00014099188829407296); 
}

/**
 * @brief Need to double check -- this is a dumb test.
 * Need to setup Responsibilities for this to be a real test.
 */
TEST(FullModel, getResponsibility)
{  //Arrange
  FullModel sut;
  sut.initBounds(7., 2., 1., 0.33, 5., 30.);
  sut.bidir.setFootprint(2.);
  // std::cout << sut.write_out() << std::endl;

  //Act
  double ans = sut.getResponsibility();

  //Assert
  EXPECT_EQ(ans, 0.);
}

/**
 * @brief Need to double check -- this is a dumb test.
 * Need to setup Responsibilities for this to be a real test.
 */
TEST(FullModel, updateParameters)
{  
  //Arrange
  FullModel sut;
  sut.initBounds(7., 2., 1., 0.33, 5., 30.);
  sut.bidir.setFootprint(2.);

  //Act
  sut.updateParameters(2.,2.);
  // std::cout << sut.write_out() << std::endl;

  //Assert
  EXPECT_EQ(sut.bidir.getMu(), 0.);
}


/**
 * @brief Need to double check -- this is a dumb test.
 */
TEST(FullModel, calculateRi)
{ 
  //Arrange
  FullModel sut;
  sut.initBounds(7., 2., 1., 0.33, 5., 30.);
  sut.bidir.setFootprint(2.);

  perStrandInfo coverage(15,25);
  perStrandInfo output;
  
  // Act
  // std::cout << output.write_out() << std::endl;
  output = sut.calculateRi(10., coverage);
  // std::cout << output.write_out() << std::endl;

  // Assert
  EXPECT_EQ(output.forward, 3.0);
  EXPECT_EQ(output.reverse, 3.0);
}