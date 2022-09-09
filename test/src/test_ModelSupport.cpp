/**
 * @file test_ModelSupport.cpp 
 * @author Robin Dowell
 * @brief Testing the distribution classes
 * @date 2022-06-16
 * 
 */
#include "gmock/gmock.h"
#include "ModelSupport.h"

/* Priors class is dumb container, no logic 
TEST (Priors,) { }
*/

// perStrandInfo_sumBothStrands
TEST (PerStrandInfo, sumBothStrands)
{
    // Arrange: bring SUT to desired state
    perStrandInfo sut = perStrandInfo(5,9);
    // Act: call methods on SUT, capture output
    double answer = sut.sumBothStrands();

    // Assert: Verify the outcome
    EXPECT_EQ(answer, 14);
}

// Responsibilities_resetRi
TEST(Responsibilities, resetRiCorrect) 
{
    // Arrange: bring SUT to desired state
    Responsibilities sut;
    sut.Ri.forward = 3.;
    sut.Ri.reverse = 6.;
    sut.Rk.forward = 3.;
    sut.Rk.reverse = 6.;

    // Act: call methods on SUT, capture output
    sut.resetRi();

    // Assert: Verify the outcome
    EXPECT_EQ(sut.Ri.forward, 0);
    EXPECT_EQ(sut.Ri.reverse, 0);
    EXPECT_EQ(sut.Rk.forward, 3);
    EXPECT_EQ(sut.Rk.reverse, 6);
}

// Responsibilities_getResponsibilities
TEST(Responsibilities, getResponsibilities) 
{
    // Arrange: bring SUT to desired state
    Responsibilities sut;
    sut.Ri.forward = 1.;
    sut.Ri.reverse = 2.;
    sut.Rk.forward = 3.;
    sut.Rk.reverse = 5.;

    // Act: call methods on SUT, capture output
    double sum = sut.getResponsibility();

    // Assert: Verify the outcome
    EXPECT_EQ(sum, 8);      // Should be sum of Rk
}

class sumExpectedTest: public :: testing:: Test {
  protected:
  void SetUp() override {
     sut.sumRExpY = 3;
     sut.sumRExpX = 5;
     sut.sumRExpX2 = 4;
     sut.sumRXExpY = 8;
  }
  sumExpected sut;
};

// sumEXpected 
TEST_F(sumExpectedTest, resetCorrect) 
{
    // Act: call methods on SUT, capture output
    sut.reset();

    // Assert: Verify the outcome
    EXPECT_EQ(sut.sumRExpY, 0);
    EXPECT_EQ(sut.sumRExpX, 0);
    EXPECT_EQ(sut.sumRExpX2, 0);
    EXPECT_EQ(sut.sumRXExpY, 0);
}

TEST_F(sumExpectedTest, SumRExpYModify) 
{
    // Act: call methods on SUT, capture output
    sut.addToSumRExpY(5.);

    // Assert: Verify the outcome
    EXPECT_EQ(sut.getSumRExpY(), 8);
}

TEST_F(sumExpectedTest, SumRExpXModify) 
{
    // Act: call methods on SUT, capture output
    sut.addToSumRExpX(5.);

    // Assert: Verify the outcome
    EXPECT_EQ(sut.getSumRExpX(), 10);
}

TEST_F(sumExpectedTest, SumRExpX2Modify) 
{
    // Act: call methods on SUT, capture output
    sut.addToSumRExpX2(5.);

    // Assert: Verify the outcome
    EXPECT_EQ(sut.getSumRExpX2(), 9);
}

TEST_F(sumExpectedTest, SumRXExpYModify) 
{
    // Act: call methods on SUT, capture output
    sut.addToSumRXExpY(5.);

    // Assert: Verify the outcome
    EXPECT_EQ(sut.getSumRXExpY(), 13);
}

/* bidirConstraints is dumb container with no logic
TEST(bidirConstraints, ) 
{
    // Arrange: bring SUT to desired state
    // Act: call methods on SUT, capture output
    // Assert: Verify the outcome
}
*/

