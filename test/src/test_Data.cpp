/**
 * @file test_dInterval.cpp 
 * @author Robin Dowell
 * @brief Testing the data interval class/objects:  dInterval, gInterval
 * @date 2022-04-15
 * 
 */
#include "gmock/gmock.h"
#include "Data.h"

TEST(Data, rawData_addDataPoints)
{
    // Arrange: bring SUT to desired state
    RawData sut;

    // Act: call methods on SUT, capture output
    sut.addDataPoints(10, 14, 4);
    sut.addDataPoints(13, 15, 1);
    sut.addDataPoints(8, 10, -2);

    // std::cout << sut.write_out() << std::endl;
    // std::cout << sut.data_dump() << std::endl;

    // Assert: Verify the outcome
    EXPECT_EQ(sut.forward.size(), 6);
    EXPECT_EQ(sut.Length(), 7);
}

TEST(Data, rawData_removeDup)
{
    // Arrange: bring SUT to desired state
    RawData sut;

    // Act: call methods on SUT, capture output
    sut.addDataPoints(10, 14, 4);   // Position 13 measured as both 4 and 1!
    sut.addDataPoints(13, 15, 1);
    sut.addDataPoints(8, 10, -2);

    // std::cout << sut.data_dump() << std::endl;
    sut.RemoveDuplicates();
    //std::cout << sut.data_dump() << std::endl;

    // Assert: Verify the outcome
    EXPECT_EQ(sut.minX, 8);
    EXPECT_EQ(sut.maxX, 15);
    EXPECT_EQ(sut.forward.size(), 5);
}

TEST(Data, dInt_bin)
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
}
