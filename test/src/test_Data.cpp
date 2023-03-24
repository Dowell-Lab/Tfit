/**
 * @file test_Data.cpp 
 * @author Robin Dowell
 * @brief Testing the data interval classes: \ref RawData and \ref dInterval
 * @date 2022-04-15
 * 
 */
#include "gmock/gmock.h"
#include "gInterval.h"
#include "Data.h"

/* Currently PointCov is a dumb container, no logic
TEST(PointCov, ) { }
*/

TEST(RawData, addDataPoints)
{
    // Arrange: bring SUT to desired state
    RawData sut;

    // Act: call methods on SUT, capture output
    sut.addDataPoints(10, 14, 8);
    sut.addDataPoints(8, 10, -2);

    // std::cout << sut.write_out() << std::endl;
    // std::cout << sut.data_dump() << std::endl;

    // Assert: Verify the outcome
    // Number of data points = number of positions
    EXPECT_EQ(sut.forward.size(), 4);
    EXPECT_EQ(sut.reverse.size(), 2);
}

TEST(RawData, removeDup)
{
    // Arrange: bring SUT to desired state
    RawData sut;

    // Act: call methods on SUT, capture output
    sut.addDataPoints(10, 14, 4);   // Position 13 measured as both 4 and 1!
    sut.addDataPoints(13, 15, 1);
    sut.addDataPoints(12, 14, -2); // This is on reverse so should be kept.

    // std::cout << sut.data_dump() << std::endl;
    sut.RemoveDuplicates();
    //std::cout << sut.data_dump() << std::endl;

    // Assert: Verify the outcome
    EXPECT_EQ(sut.forward.size(), 5);
    EXPECT_EQ(sut.reverse.size(), 2);
}

TEST(RawData, LengthCorrect_NoRepeats)
{
    // Arrange: bring SUT to desired state
    RawData sut;
    sut.addDataPoints(10, 14, 4);   // Position 13 measured as both 4 and 1!
    sut.addDataPoints(8, 10, -2);

    // Assert: Verify the outcome
    EXPECT_EQ(sut.Length(), 6);
}

//  void Sort();    //!< Sorts both forward and reverse by coordinate
TEST(RawData, CorrectlySorted)
{
    // Arrange: bring SUT to desired state
    RawData sut;
    sut.addDataPoints(14, 16, 1);
    sut.addDataPoints(10, 14, -2);  
    sut.addDataPoints(2, 10, 3);
    sut.addDataPoints(3, 6, -1);

    // Act: function we are testing 
    sut.Sort();

    // Assert: Verify the outcome
    bool outoforder = false;
    for (unsigned int i = 0; i < sut.forward.size(); i++) {
      if (i > 0) { // Ingores the first one
         if (sut.forward[i-1].coordinate > sut.forward[i].coordinate) {
             outoforder = true;
         }
      }
    }
    EXPECT_FALSE(outoforder);

    outoforder = false;     // Now for reverse strand
    for (unsigned int i = 0; i < sut.reverse.size(); i++) {
      if (i > 0) { // Ingores the first one
         if (sut.reverse[i-1].coordinate > sut.reverse[i].coordinate) {
             outoforder = true;
         }
      }
    }
    EXPECT_FALSE(outoforder);
}

/********* Test dInterval: ********/

TEST (dInterval, LengthCorrectNoData)
{
    // Arrange: bring SUT to desired state
    dInterval sut = dInterval();;

    // Assert: Verify the outcome
    EXPECT_EQ(sut.getLength(), 0);
}

/*
TEST (dInterval, LengthCorrectWithRaw)
{
    // Arrange: bring SUT to desired state
    RawData data;
    data.addDataPoints(10, 14, 4);   
    data.addDataPoints(14, 15, 1);
    data.addDataPoints(8, 10, -2);

    // Act: call methods on SUT, capture output
    dInterval sut(&data,2,1);
    // std::cout << sut.data_dump() << std::endl;

    // Assert: Verify the outcome
    EXPECT_EQ(sut.getLength(), 7);
}

TEST (dInterval, LengthCorrectNoRaw)
{
    // Arrange: bring SUT to desired state
    RawData data;
    data.addDataPoints(10, 14, 4);   
    data.addDataPoints(14, 15, 1);
    data.addDataPoints(8, 10, -2);
    dInterval sut(&data,2,1);

    // Assert: Verify the outcome
    EXPECT_EQ(sut.getLength(), 7);
}
*/


class dIntervalSumTests: public :: testing::Test {
  protected:
  void SetUp() override {
    data.addDataPoints(4, 14, 4);   
    data.addDataPoints(14, 23, 3);
    data.addDataPoints(4, 14, -2);   
    data.addDataPoints(14, 23, -6);
  }
  RawData data;
};

// double sumForward();
TEST_F(dIntervalSumTests, sumForward)
{  
  // Arrange: bring SUT to desired state
  dInterval sut(2,1);
  sut.createFromRaw(&data);

  // Assert: Verify the outcome
  EXPECT_EQ(sut.sumForward(), 67);  // 10*4 + 9*3
}

// double sumReverse();
TEST_F(dIntervalSumTests, sumReverse)
{
  // Arrange: bring SUT to desired state
  dInterval sut(2,1);
  sut.createFromRaw(&data);

    // Assert: Verify the outcome
    EXPECT_EQ(sut.sumReverse(), 74); 
}

// double sumInterval(int, int, char);  // on single strand
TEST_F(dIntervalSumTests, sumIntervalPosStrand)
{
  // Arrange: bring SUT to desired state
  dInterval sut(2,1);
  sut.createFromRaw(&data);

  // Note these are indexes into dinterval, not genomic positions
  double sum = sut.sumInterval(2, 4, '+');

  // Assert: Verify the outcome
  EXPECT_EQ(sum, 16);  // 
}

TEST_F (dIntervalSumTests, sumIntervalNegStrand)
{
  // Arrange: bring SUT to desired state
  dInterval sut(2,1);
  sut.createFromRaw(&data);

  // Note these are indexes into dinterval, not genomic positions
  double sum = sut.sumInterval(2, 4, '-');

    // Assert: Verify the outcome
    EXPECT_EQ(sum, 8);  // 
}

//  double sumAlldata();  // both strands
TEST_F (dIntervalSumTests, sumAlldata)
{
  // Arrange: bring SUT to desired state
  dInterval sut(2,1);
  sut.createFromRaw(&data);

  // Note these are indexes into dinterval, not genomic positions
  double sum = sut.sumAlldata();

  // Assert: Verify the outcome
  EXPECT_EQ(sum, 141);  // 
}

