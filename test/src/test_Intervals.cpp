/**
 * @file test_Intervals.cpp 
 * @author Robin Dowell
 * @brief Testing the data interval class/objects:  dInterval, gInterval
 * @date 2022-01-27
 * 
 */
#include "gmock/gmock.h"
#include "Intervals.h"

TEST(Intervals, gInt_basicBED4)
{
    // Arrange: bring SUT to desired state
    std::string chromosome, name;
    name = "TestName";
    chromosome = "chrTest";
    gInterval sut = gInterval(chromosome, 100, 1000, name); 

    // Act: call methods on SUT, capture output
    name = sut.write_out(); 

    // Assert: Verify the outcome
    EXPECT_THAT(name, "#chrTest:100-1000,TestName");
}

TEST(Intervals, gInt_basicBED6)
{
    // Arrange: bring SUT to desired state
    std::string chromosome, name, strand;
    name = "TestName";
    chromosome = "chrTest";
    strand = ".";
    bed6 sut = bed6(chromosome, 100, 1000, name, 30, strand); 

    // Act: call methods on SUT, capture output
    name = sut.write_out(); 

    // Assert: Verify the outcome
    EXPECT_THAT(name, "#chrTest:100-1000,TestName,30,+");
}

// Stupid simplest first case:  Just to get skeleton class setup
TEST(Intervals, dInt_RecoverIdentity)
{
    // Arrange: bring SUT to desired state
    std::string name;
    name = "TestName";
    dInterval sut = dInterval(name, 100, 1000); 

    // Act: call methods on SUT, capture output
    name = sut.write_out(); 

    // Assert: Verify the outcome
    EXPECT_THAT(name, "#TestName:d:100-1000,1,0");
}
