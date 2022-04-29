/**
 * @file test_dInterval.cpp 
 * @author Robin Dowell
 * @brief Testing the data interval class/objects:  dInterval, gInterval
 * @date 2022-04-15
 * 
 */
#include "gmock/gmock.h"
#include "Data.h"

TEST(Intervals, dInt_RecoverIdentity)
{
    // Arrange: bring SUT to desired state
    dInterval sut = dInterval();

    // Act: call methods on SUT, capture output
    std::string name = sut.write_out(); 

    // Assert: Verify the outcome
    EXPECT_THAT(name, "\t-1,1,1");
}
