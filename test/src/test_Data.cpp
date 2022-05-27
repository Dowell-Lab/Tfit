/**
 * @file test_Data.cpp 
 * @author Robin Dowell
 * @brief Testing the data interval classes: \ref RawData and \ref dInterval
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

TEST(Data, dInt_bin_scale)
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
    EXPECT_EQ(sut.N, 21);
}

TEST(Data, coordinateTranslation)
{
    // Arrange: bring SUT to desired state
    RawData data;
    data.addDataPoints(10, 14, 4);   
    data.addDataPoints(13, 15, 1);
    data.addDataPoints(8, 10, -2);
    data.addDataPoints(11,17, 0);

    dInterval sut(&data,2,1);
    std::cout << sut.data_dump() << std::endl;

    // Act: call methods on SUT, capture output
    double result1 = sut.getGenomeCoordfromIndex(2);
    double result2 = sut.getIndexfromGenomic(result1);

    double result4 = sut.getIndexfromGenomic(13);
    double result3 = sut.getGenomeCoordfromIndex(result4);

    // Assert: Verify the outcome
    EXPECT_EQ(result2, 2);      // index -> genomic -> index
    EXPECT_EQ(result3, 13);     // genomic -> index -> genomic
}
