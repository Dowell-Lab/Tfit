/**
 * @file test_SetROI.cpp 
 * @author Robin Dowell
 * @brief Testing the set of ROI data class
 * @date 2022-04-15
 * 
 */
#include "gmock/gmock.h"
#include "Intervals.h"
#include "Data.h"
#include "Regions.h"

TEST(SetROI, constructAdd)
{
    // Arrange: bring SUT to desired state
    SetROI sut;
    bed6 testregion("Example1", 2000, 4000, "ID1", 300, "+");

    // Act: call methods on SUT, capture output
    sut.addRegionToSet(&testregion);

    // Assert: Verify the outcome
    EXPECT_EQ(sut.chr_names.num_elements, 1);
}
