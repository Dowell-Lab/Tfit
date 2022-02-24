/**
 * @file test_dIntervals.cpp 
 * @author Robin Dowell
 * @brief Testing the data interval class/object 
 * @date 2022-01-27
 * 
 */
#include "gmock/gmock.h"
#include "dInterval.h"

// Stupid simplest first case:  Just to get skeleton class setup
TEST(dInterval, RecoverIdentity)
{
    // Arrange: bring SUT to desired state
    std::string name;
    name = "TestName";
    dInterval sut = dInterval(name, 100, 1000); 

    // Act: call methods on SUT, capture output
    name = sut.write_out(); 

    // Assert: Verify the outcome
    EXPECT_THAT(name, "#TestName:100-1000,1,0");
}

TEST(dInterval, IterateValues)
{
  // Arrange
  // Act
  // Assert
}
