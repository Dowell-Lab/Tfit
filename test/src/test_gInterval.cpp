/**
 * @file test_gInterval.cpp 
 * @author Robin Dowell
 * @brief Testing for gInterval container
 * @date 2023-01-18
 * 
 */
#include "gmock/gmock.h"
#include "Bed.h"

TEST(gInterval, basic) {
    // Arrange: bring SUT to desired state
    gInterval container;
    // std::cout << name << std::endl;

    // Act: call methods on SUT, capture output

    // Assert: Verify the outcome
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
*/

/*
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
*/

/*  Recover when this function comes back.  
class getWithinRangePosition: public :: testing::Test {
   protected:
   void SetUp() override {
        // Pos strand
    data.addDataPoints(1, 10, 4);   
    data.addDataPoints(15, 18, 1);
    data.addDataPoints(22, 24, 2);
        // Neg strand
    data.addDataPoints(1, 7, -1);   
    data.addDataPoints(11, 13, -2);
    data.addDataPoints(26, 28, -3);
   }   
   RawData data;
};

// int getWithinRangeofPosition(double position, double dist);
TEST_F (getWithinRangePosition, PosStrandInRange)
{
    // Arrange: bring SUT to desired state
   dInterval sut(&data,2,1);

    // Act: call methods on SUT, capture output
    double position = sut.getWithinRangeofPosition(11, 2);

    // Assert: Verify the outcome
    EXPECT_EQ(sut.reverse.size(), 2);
}
*/

