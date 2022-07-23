/**
 * @file test_Data.cpp 
 * @author Robin Dowell
 * @brief Testing the data interval classes: \ref RawData and \ref dInterval
 * @date 2022-04-15
 * 
 */
#include "gmock/gmock.h"
#include "Data.h"

/* Currently PointCov is a dumb container, no logic
TEST(PointCov, ) 
{
}
*/

  // void addDataPoints(double st, double sp, double cov);
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

//   void RemoveDuplicates();  //!< Removes any coordinate with zero reads on both strands
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

  // double Length();    
TEST(RawData, LengthCorrect_NoRepeats)
{
    // Arrange: bring SUT to desired state
    RawData sut;
    sut.addDataPoints(10, 14, 4);   // Position 13 measured as both 4 and 1!
    sut.addDataPoints(8, 10, -2);

    // Assert: Verify the outcome
    EXPECT_EQ(sut.Length(), 6);
}

/*
//  void Sort();    //!< Sorts both forward and reverse by coordinate
TEST(RawData, CorrectlySorted)
{
    // Arrange: bring SUT to desired state
    RawData sut;
    sut.addDataPoints(14, 16, 1);
    sut.addDataPoints(10, 14, 4);  
    sut.addDataPoints(2, 10, 3);

    // std::cout << sut.data_dump() << std::endl;
    sut.Sort();
    //std::cout << sut.data_dump() << std::endl;

    // Assert: Verify the outcome
    // Coordinates on forward should be in numerical order.
}
*/


/********* Test dInterval: ********/
/*
// double getLength();
TEST (dInterval, LengthCorrect)
{

}

// double sumForward();
TEST (dInterval, sumForward)
{

}

// double sumReverse();
TEST (dInterval, sumReverse)
{

}

// double sumInterval(int, int, char);  // on single strand
TEST (dInterval, sumInterval)
{

}

//  double sumAlldata();  // both strands
TEST (dInterval, sumAlldata)
{

}

*/

TEST(dInterval, CreateBinScaleWithDuplicates)
{
    // Arrange: bring SUT to desired state
    RawData data;
    data.addDataPoints(10, 14, 4);   
    data.addDataPoints(13, 15, 1);
    data.addDataPoints(8, 10, -2);
    data.addDataPoints(11,17, 0);

    // Act: call methods on SUT, capture output
    dInterval sut(&data,2,1);
    // std::cout << sut.data_dump() << std::endl;

    // Assert: Verify the outcome
    EXPECT_EQ(sut.bins, 4);
    EXPECT_EQ(sut.sumAlldata(), 21);
}

class CoordinateTransformTest: public :: testing::Test {
    protected:
    void SetUp() override {
      data.addDataPoints(10, 14, 4);   
      data.addDataPoints(14, 15, 1);
      data.addDataPoints(8, 10, -2);
      data.addDataPoints(1,8, 0);
    }
    // void TearDown() override { }
    RawData data;
};

TEST_F(CoordinateTransformTest, Genome2Index2Genome)
{
    // Arrange: bring SUT to desired state
    dInterval sut(&data,2,1);

    // std::cout << sut.data_dump() << std::endl;

    // Act: call methods on SUT, capture output
    double result1 = sut.getGenomeCoordfromIndex(2);
    double result2 = sut.getIndexfromGenomic(result1);

    // Assert: Verify the outcome
    EXPECT_EQ(result2, 2);      // index -> genomic -> index
}

TEST_F(CoordinateTransformTest, Index2Genome2Index)
{
    // Arrange: bring SUT to desired state
    dInterval sut(&data,2,1);

    // Act: call methods on SUT, capture output
    double result4 = sut.getIndexfromGenomic(13);
    double result3 = sut.getGenomeCoordfromIndex(result4);

    // Assert: Verify the outcome
    EXPECT_EQ(result3, 13);     // genomic -> index -> genomic
}

//  int getIndexfromData(double);  // given data coordinate, get closest index
//  double getGenomefromData(double);
// double getDataCoordfromIndex(int);
// double getDataCoordfromGenomeCoord(double);
TEST_F(CoordinateTransformTest, coordData2Index)
{
    // Arrange: bring SUT to desired state
    dInterval sut(&data,2,1);

    // Act: call methods on SUT, capture output
    int result1 = sut.getIndexfromData(4);

    // Assert: Verify the outcome
    EXPECT_EQ(result1, 2);      // index -> genomic -> index
}

/*

// int getWithinRangeofPosition(double position, double dist);
TEST (dInterval, getWithinRangeofPosition)
{

}

// void DeallocateX();   // Deallocates X, leaves other variables intact.

*/