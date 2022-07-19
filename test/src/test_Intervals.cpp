/**
 * @file test_Intervals.cpp 
 * @author Robin Dowell
 * @brief Testing the data interval classes: \ref gInterval \ref bed6
 * @date 2022-01-27
 * 
 */
#include "gmock/gmock.h"
#include "Intervals.h"

/* Need to test I/O for bed3, bed4, bed6 and bed12.  */

TEST(gInterval, writeBedEQreadBed)
{
    // Arrange: bring SUT to desired state
    gInterval temp = gInterval("TestName", 100, 1000, "chrTest"); 
    std::string name = temp.write_asBED();  // Writes as bed4

    gInterval sut = gInterval();

    // Act: call methods on SUT, capture output
    sut.setfromBedLine(name);       // Reads as bed4

    // Assert: Verify the outcome
    EXPECT_THAT(temp.write_out(), sut.write_out());
}

class Bed6ContainsTest: public :: testing::Test {
  protected:
  void SetUp() override { }
  // void TearDown() override 
  bed6 sut = bed6("TestName", 100, 1000, "chrTest", 30, "."); 
};

TEST_F(Bed6ContainsTest, containsEdgeCase) {
  EXPECT_TRUE(sut.Contains((double)100)); // Edge case
}

TEST_F(Bed6ContainsTest, containsValue) {
  EXPECT_TRUE(sut.Contains((double)300)); // Edge case
}

TEST_F(Bed6ContainsTest, containsEdgeCaseHalfOpen) {
  EXPECT_FALSE(sut.Contains((double)1000)); // Edge case
}

TEST_F(Bed6ContainsTest, containsBelowRange) {
  EXPECT_FALSE(sut.Contains((double)50));   // Below range 
}

TEST_F(Bed6ContainsTest, containsAboveRange) {
  EXPECT_FALSE(sut.Contains((double)2000));   // Above range 
}

TEST(bed6, overlapTestingEdge)
{
    // Arrange: bring SUT to desired state
    bed6 OR1 = bed6("TestChr", 100, 1000, "OR1", 30, "."); 
    bed6 NR1 = bed6("TestChr", 300, 1000, "NR1", 30, ".");  // half open, so no overlap!
    bed6 ER1 = bed6("TestChr", 299, 500, "ER1", 30, "."); 

    bed6 sut = bed6("TestChr", 50, 300, "R2", 30, "."); 

    // Assert: Verify the outcome
    EXPECT_TRUE(sut.Overlap((gInterval *)&OR1));
    EXPECT_TRUE(sut.Overlap((gInterval *)&ER1));
    EXPECT_FALSE(sut.Overlap((gInterval *)&NR1));
}

TEST(bed6, writeBedEQreadBed)
{
    // Arrange: bring SUT to desired state
    bed6 temp = bed6("TestName", 100, 1000, "chrTest", 30, "."); 
    std::string name = temp.write_asBED();      // Write as bed6
    // std::cout << name << std::endl;

    bed6 sut = bed6();

    // Act: call methods on SUT, capture output
    sut.setfromBedLine(name);   // Read as bed6

    // Assert: Verify the outcome
    EXPECT_THAT(temp.write_out(), sut.write_out());
}
