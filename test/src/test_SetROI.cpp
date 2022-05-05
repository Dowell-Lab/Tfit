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

TEST(SetROI, search)
{
    // Arrange: bring SUT to desired state
    SetROI sut;
    bed6 testregion("Example1", 2000, 4000, "ID1", 300, "+");
    sut.addRegionToSet(&testregion);
    bed6 query("Example1", 3000, 3500, "test", 100, "+");

    std::vector<gInterval *> results;

    // Act: call methods on SUT, capture output
    results = sut.findOverlapIntervals(&query);

    // std::cout << results.size() << std::endl;
    // std::vector<gInterval *>::iterator it;
    // for (it = results.begin(); it != results.end(); it++) {
     //  std::cout << (*it)->write_out() << std::endl;
    // }

    // Assert: Verify the outcome
    EXPECT_EQ(results.size(), 1);
}

TEST(SetROI, addDataToExistingROI)
{
    // Arrange: bring SUT to desired state
    SetROI sut;
    bed6 testregion("Chr1", 1000, 4000, "ID1", 300, "+");
    sut.addRegionToSet(&testregion);

    sut.addDataToExistingROI("Chr1", 1000, 2000, 3);
    sut.addDataToExistingROI("Chr1", 1500, 4000, -3);
    sut.addDataToExistingROI("Chr1", 2000, 3000, 1);

    // Act: call methods on SUT, capture output

    // Assert: Verify the outcome
    EXPECT_EQ(sut.regions.size(), 0);
}

