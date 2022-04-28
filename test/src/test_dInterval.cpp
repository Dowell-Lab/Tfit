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
    std::string name;
    name = "TestName";
    dInterval sut = dInterval(name);

    // Act: call methods on SUT, capture output
    name = sut.write_out(); 

    // Assert: Verify the outcome
    EXPECT_THAT(name, "#TestName:min:0:max:0:num_bins:0:total:0");
}
