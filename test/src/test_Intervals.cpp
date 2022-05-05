/**
 * @file test_Intervals.cpp 
 * @author Robin Dowell
 * @brief Testing the data interval classes: \ref gInterval \ref bed6
 * @date 2022-01-27
 * 
 */
#include "gmock/gmock.h"
#include "Intervals.h"

TEST(Intervals, basicWrite)
{
    // Arrange: bring SUT to desired state
    gInterval sut = gInterval("chrTest", 100, 1000, "TestName"); 

    // Act: call methods on SUT, capture output
    std::string name = sut.write_out(); 

    // Assert: Verify the outcome
    EXPECT_THAT(name, "#chrTest:100-1000,TestName");
}

TEST(Intervals, basicIO)
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

TEST(Intervals, bed6IO)
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

TEST(Intervals, contains)
{
    // Arrange: bring SUT to desired state
    bed6 sut = bed6("TestName", 100, 1000, "chrTest", 30, "."); 
    // std::cout << name << std::endl;

    // Assert: Verify the outcome
    EXPECT_TRUE(sut.Contains((double)100)); // Edge case
    EXPECT_TRUE(sut.Contains((double)300));   // Middle value
    EXPECT_FALSE(sut.Contains((double)1000));   // Edge value, half open system!
    EXPECT_FALSE(sut.Contains((double)50));   // Below range 
    EXPECT_FALSE(sut.Contains((double)2000));   // Above range 
}

TEST(Intervals, overlaps)
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

/*
TEST(Intervals, addData)
{
    // Arrange: bring SUT to desired state
    bed6 temp = bed6("TestChr", 100, 1000, "ID", 30, "+"); 

    bed6 sut = bed6();

    // Act: call methods on SUT, capture output
    sut.setfromBedLine(name);

    // Assert: Verify the outcome
    EXPECT_THAT(temp.write_out(), sut.write_out());
}
*/