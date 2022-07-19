/**
 * @file test_Intervals.cpp 
 * @author Robin Dowell
 * @brief Testing the data interval classes: \ref gInterval \ref bed6
 * @date 2022-01-27
 * 
 */
#include "gmock/gmock.h"
#include "Intervals.h"

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

TEST(gInterval, addDataPoint_contained)
{
    // Arrange: bring SUT to desired state
    gInterval sut = gInterval("TestName", 100, 1000, "chrTest"); 

    // Act: call methods on SUT, capture output
    sut.addDataPoint(120,150,4,false);      // Note calls rawData::addDataPoints, test elsewhere

    // Assert: Verify the outcome
    // Confirming that the coordinates were properly adjusted (when allowed)
    // Check that the RawData contains what its supposed to should be tested by RawData!
    EXPECT_TRUE(sut.data != NULL);  // addDataPoint created a link
    EXPECT_EQ(sut.start, (double)100);
    EXPECT_EQ(sut.stop, (double)1000);
}

TEST(gInterval, addDataPoint_edgeCaseNoExpand)
{
    // Arrange: bring SUT to desired state
    gInterval sut = gInterval("TestName", 100, 1000, "chrTest"); 

    // Act: call methods on SUT, capture output
    sut.addDataPoint(80,150,4,false);      // Note calls rawData::addDataPoints, test elsewhere

    // Assert: Verify the outcome
    // Confirming that the coordinates were properly adjusted (when allowed)
    EXPECT_EQ(sut.start, (double)100);
}

TEST(gInterval, addDataPoint_edgeCaseStartExpand)
{
    // Arrange: bring SUT to desired state
    gInterval sut = gInterval("TestName", 100, 1000, "chrTest"); 

    // Act: call methods on SUT, capture output
    sut.addDataPoint(80,150,4,true);      // Note calls rawData::addDataPoints, test elsewhere

    // Assert: Verify the outcome
    // Confirming that the coordinates were properly adjusted (when allowed)
    EXPECT_EQ(sut.start, (double)80);
}

TEST(gInterval, addDataPoint_edgeCaseStopExpand)
{
    // Arrange: bring SUT to desired state
    gInterval sut = gInterval("TestName", 100, 1000, "chrTest"); 

    // Act: call methods on SUT, capture output
    sut.addDataPoint(900,1500,4,true);      // Note calls rawData::addDataPoints, test elsewhere

    // Assert: Verify the outcome
    // Confirming that the coordinates were properly adjusted (when allowed)
    EXPECT_EQ(sut.stop, (double)1500);
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

/*
TEST(bed12, writeBedEQreadBed) {
    // Arrange: bring SUT to desired state
    bed12 temp = bed12("TestName", 100, 1000, "chrTest", 30, "."); 
    std::string name = temp.write_asBED();      // Write as bed6
    // std::cout << name << std::endl;

    bed12 sut = bed12();

    // Act: call methods on SUT, capture output
    sut.setfromBedLine(name);   // Read as bed6

    // Assert: Verify the outcome
    EXPECT_THAT(temp.write_out(), sut.write_out());
}
*/