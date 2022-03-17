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

TEST(Intervals, gInt_BED4io)
{
    // Arrange: bring SUT to desired state
    std::string chromosome, name;
    name = "TestName";
    chromosome = "chrTest";
    gInterval temp = gInterval(chromosome, 100, 1000, name); 
    gInterval sut = gInterval();

    // Act: call methods on SUT, capture output
    name = temp.write_asBED();
    sut.setfromBedLine(name);

    // Assert: Verify the outcome
    EXPECT_THAT(temp.write_out(), sut.write_out());
}

TEST(Intervals, gInt_basicBED6)
{
    // Arrange: bring SUT to desired state
    std::string chromosome, name, strand;
    name = "TestName";
    chromosome = "chrTest";
    strand = "+";
    bed6 sut = bed6(chromosome, 100, 1000, name, 30, strand); 

    // Act: call methods on SUT, capture output
    name = sut.write_out(); 

    // Assert: Verify the outcome
    EXPECT_THAT(name, "#chrTest:100-1000,TestName,30,+");
}

TEST(Intervals, gInt_BED6io)
{
    // Arrange: bring SUT to desired state
    std::string chromosome, name, strand;
    name = "TestName";
    chromosome = "chrTest";
    strand = ".";
    bed6 temp = bed6(chromosome, 100, 1000, name, 30, strand); 
    bed6 sut = bed6();

    // Act: call methods on SUT, capture output
    name = temp.write_asBED();
    // std::cout << name << std::endl;
    sut.setfromBedLine(name);

    // Assert: Verify the outcome
    EXPECT_THAT(temp.write_out(), sut.write_out());
}

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
