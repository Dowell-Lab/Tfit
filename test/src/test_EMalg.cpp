/**
 * @file test_EMalg.cpp 
 * @author Robin Dowell
 * @brief Testing the Bedfile and Bedgraph classes: \ref Bedfile Bedgraph
 * @date 2022-07-26
 * 
 */
#include "gmock/gmock.h"
#include "EMalg.h"

TEST(EMalg, fitFailswithNoData) 
{
    // Arrange: bring SUT to desired state
    EMalg sut;

    // Act: call methods on SUT, capture output
    bool success = sut.fit();

    // Assert: Verify the outcome
    EXPECT_FALSE(success);
}

/*
TEST(EMalg, computeBackground_NoiseOnly) 
{
  // Arrange: bring SUT to desired state
  RawData data;
  data.addDataPoints(10, 14, 4);   
  data.addDataPoints(14, 15, 1);
  data.addDataPoints(8, 10, -2);
  dInterval interval(&data,2,1);
 
  EMalg sut;
  sut.data = &interval;

  // Act: call methods on SUT, capture output
  bool success = sut.fit();

  // Assert: Verify the outcome
  EXPECT_TRUE(success);
  EXPECT_DOUBLE_EQ(sut.models.ll, -51.089280028918225);
}
*/

TEST(EMalg, Itest_fit_oneModel) 
{
  // Arrange: bring SUT to desired state
  //Act
  //Assert
}