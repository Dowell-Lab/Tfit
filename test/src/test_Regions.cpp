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

/* 
  std::string print_tree_at_chromosome(std::string chromo);
  std::string write_out();
   
  void createSearchIndex();   // Builds the interval trees given the regions 
  std::vector<gInterval *>findOverlapIntervals(gInterval *);
  void clearTrees();

  void addDataCreateROI(std::string chr, double start, double stop, double coverage);
  bool addDataToExistingROI(std::string chr, double start, double stop, double coverage);
  void ConditionDataSet(int v_delta, int v_scale);
  */

TEST(SetROI, addRegionToSet_IncrementNumElements)
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


TEST(SetROI, findOverlapInterval_NoContents)
{
    // Arrange: bring SUT to desired state
    SetROI sut;
    std::vector<gInterval *> results;
    bed6 query = bed6("NewChr", 3200, 3500, "test", 100, "+");

    // Act: call methods on SUT, capture output
    results = sut.findOverlapIntervals(&query);

    // Assert: Verify the outcome
    EXPECT_EQ(results.size(), 0);
}

class findOverlapIntervalTests: public ::testing::Test {
  protected:
  void SetUp() override {
    region1 = new bed6("Example1", 2000, 4000, "ID1", 300, "+");
    region2 = new bed6("Example2", 1000, 5000, "ID2", 300, "-");
    region3 = new bed6("Example1", 3000, 5000, "ID3", 300, "+");
    sut.addRegionToSet(region1);
    sut.addRegionToSet(region2);
    sut.addRegionToSet(region3);
  }
  void TearDown() override {
    sut.clearTrees();
  }
  SetROI sut;
  bed6 *region1, *region2, *region3;
  std::vector<gInterval *> results;
};

TEST_F(findOverlapIntervalTests, SingleOverlap)
{
    // Arrange: bring SUT to desired state
    bed6 query = bed6("Example1", 1000, 2500, "test", 100, "+");

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

TEST_F(findOverlapIntervalTests, MultipleOverlap)
{
    // Arrange: bring SUT to desired state
    bed6 query = bed6("Example1", 3200, 3500, "test", 100, "+");

    // Act: call methods on SUT, capture output
    results = sut.findOverlapIntervals(&query);

    // Assert: Verify the outcome
    EXPECT_EQ(results.size(), 2);
}

TEST_F(findOverlapIntervalTests, ChromDoesNotExist)
{
    // Arrange: bring SUT to desired state
    bed6 query = bed6("NewChr", 3200, 3500, "test", 100, "+");

    // Act: call methods on SUT, capture output
    results = sut.findOverlapIntervals(&query);

    // Assert: Verify the outcome
    EXPECT_EQ(results.size(), 0);
}


/*
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
*/

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

