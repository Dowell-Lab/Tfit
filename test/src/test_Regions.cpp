/**
 * @file test_Regions.cpp 
 * @author Robin Dowell
 * @brief Testing the set of ROI data class: \ref SetROI
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
    bed6 *testregion;
    testregion = new bed6("Example1", 2000, 4000, "ID1", 300, "+");

    // Act: call methods on SUT, capture output
    sut.addRegionToSet(testregion);

    // Assert: Verify the outcome
    EXPECT_EQ(sut.chr_names.num_elements, 1);

}

TEST(SetROI, search)
{
    // Arrange: bring SUT to desired state
    SetROI sut;
    bed6 *testregion;
    testregion = new bed6("Example1", 2000, 4000, "ID1", 300, "+");
    sut.addRegionToSet(testregion);

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

    sut.clearTrees();
}

TEST(SetROI, addDataToExistingROI)
{
    // Arrange: bring SUT to desired state
    SetROI sut;
    bed6 *testregion;
    testregion = new bed6("Chr1", 10, 40, "ID1", 300, "+");
    sut.addRegionToSet(testregion);        // The interval (bed)

    sut.addDataToExistingROI("Chr1", 10, 20, 3);        // The data (bedGraph)

    int idx = sut.chr_names.lookupIndex("Chr1");    // Index for Chr1 
    std::vector<gInterval *> outRegions;
    outRegions = sut.regions[idx];

    //std::cout << outRegions[0]->data->data_dump() << std::endl;

    // Act: call methods on SUT, capture output
    sut.ConditionDataSet(2,1);      // Creates dIntervals for all regions.
    //std::cout << outRegions[0]->data->cdata->data_dump() << std::endl;

    // Assert: Verify the outcome
    EXPECT_EQ(outRegions[0]->chromosome, "Chr1");
    if (outRegions[0]->data != NULL) {
      EXPECT_EQ(outRegions[0]->data->minX, 10);
      if (outRegions[0]->data->cdata != NULL) {
          EXPECT_EQ(outRegions[0]->data->cdata->sumAlldata(), 30);
      } else {
        FAIL() << "We should have conditioned data!";
      }
    } else {
        FAIL() << "We should have data!";
    }
    sut.clearTrees();
}

TEST(SetROI, addDataCreateROI)
{
    // Arrange: bring SUT to desired state
    SetROI sut;

    sut.addDataCreateROI("Chr1", 10, 20, -3);        // The data (bedGraph)
            // Note now using negative strand!

    int idx = sut.chr_names.lookupIndex("Chr1");    // Index for Chr1 
    std::vector<gInterval *> outRegions;
    outRegions = sut.regions[idx];

    //std::cout << outRegions[0]->data->data_dump() << std::endl;

    // Act: call methods on SUT, capture output
    sut.ConditionDataSet(2,1);      // Creates dIntervals for all regions.
    //std::cout << outRegions[0]->data->cdata->data_dump() << std::endl;

    // Assert: Verify the outcome
    EXPECT_EQ(outRegions[0]->chromosome, "Chr1");
    if (outRegions[0]->data != NULL) {
      EXPECT_EQ(outRegions[0]->data->maxX, 20);     // This time check max!
      if (outRegions[0]->data->cdata != NULL) {
          EXPECT_EQ(outRegions[0]->data->cdata->sumAlldata(), 30);
      } else {
        FAIL() << "We should have conditioned data!";
      }
    } else {
        FAIL() << "We should have data!";
    }
}

