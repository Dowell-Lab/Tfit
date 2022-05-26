/**
 * @file reproduce_segment.cpp 
 * @author Robin Dowell
 * @brief Testing the old (Joey's) \ref segment class (to be refactored) 
 * @date 2022-02-01
 * 
 */
#include "gmock/gmock.h"
#include "load.h"  // contains segment class

TEST(segment, LengthEquivalent)
{
    // Arrange: bring SUT to desired state
    segment sut = segment("chrTest", 100, 1000, 3, "+"); 

    // Act: call methods on SUT, capture output
    double length = sut.getXLength(); 
    double oldLength = sut.maxX - sut.minX;

    // Assert: Verify the outcome
    EXPECT_EQ(length, oldLength);
}
