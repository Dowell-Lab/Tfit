/**
 * @file test_Regions.cpp 
 * @author Robin Dowell
 * @brief Testing the set of ROI data class: \ref SetROI
 * @date 2022-04-15
 * 
 */
#include "gmock/gmock.h"
#include "Data.h"
#include "Regions.h"

TEST(SetROI, addRegionToSet_IncrementNumElements)
{
    // Arrange: bring SUT to desired state
    SetROI sut;
    bed12 *testregion;
    testregion = new bed12("Example1", 2000, 4000, "ID1");

    // Act: call methods on SUT, capture output
    sut.addRegionToSet(testregion);

    // Assert: Verify the outcome
    EXPECT_EQ(sut.chr_names.num_elements, 1);
}

TEST(SetROI, findOverlapInterval_NoContents)
{
    // Arrange: bring SUT to desired state
    SetROI sut;
    std::vector<bed12 *> results;
    bed12 query("NewChr", 3200, 3500, "test");

    // Act: call methods on SUT, capture output
    results = sut.findOverlapIntervals(&query);

    // Assert: Verify the outcome
    EXPECT_EQ(results.size(), 0);
}

class SetROITest: public ::testing::Test {
  protected:
  void SetUp() override {
    region1 = new bed12("Example1", 2000, 4000, "ID1");
    region2 = new bed12("Example2", 1000, 5000, "ID2");
    region3 = new bed12("Example1", 3000, 5000, "ID3");
    sut.addRegionToSet(region1);
    sut.addRegionToSet(region2);
    sut.addRegionToSet(region3);
  }
  void TearDown() override {
    sut.clearTrees();
  }
  SetROI sut;
  bed12 *region1, *region2, *region3;
  std::vector<bed12 *> results;
};

TEST_F(SetROITest, findOverlapIntervalTests_SingleOverlap)
{
    // Arrange: bring SUT to desired state
    bed12 query("Example1", 1000, 2500, "test");

    // Act: call methods on SUT, capture output
    results = sut.findOverlapIntervals(&query);

    // Assert: Verify the outcome
    EXPECT_EQ(results.size(), 1);
}

TEST_F(SetROITest, findOverlapIntervalTests_MultipleOverlap)
{
    // Arrange: bring SUT to desired state
    bed12 query("Example1", 3200, 3500, "test");

    // Act: call methods on SUT, capture output
    results = sut.findOverlapIntervals(&query);

    // Assert: Verify the outcome
    EXPECT_EQ(results.size(), 2);
}

TEST_F(SetROITest, findOverlapIntervalTests_ChromDoesNotExist)
{
    // Arrange: bring SUT to desired state
    bed12 query("NewChr", 3200, 3500, "test");

    // Act: call methods on SUT, capture output
    results = sut.findOverlapIntervals(&query);

    // Assert: Verify the outcome
    EXPECT_EQ(results.size(), 0);
}

TEST_F(SetROITest, addDataToExistingROI_contained)
{
  // Act: call methods on SUT, capture output
  bool foundexisting = sut.addDataToExistingROI("Example2", 1010, 1020, 3);  // The data (bedGraph)

  // Assert: Verify the outcome
  EXPECT_TRUE(foundexisting);
}

TEST_F(SetROITest, addDataToExistingROI_edgeCase)
{
  // Act: call methods on SUT, capture output
  bool foundexisting = sut.addDataToExistingROI("Example2", 900, 1020, 3);  // The data (bedGraph)

  // Assert: Verify the outcome
  EXPECT_TRUE(foundexisting);
}

TEST_F(SetROITest, addDataToExistingROI_notFound)
{
  // Act: call methods on SUT, capture output
  bool foundexisting = sut.addDataToExistingROI("Example2", 500, 550, 3);  // The data (bedGraph)

  // Assert: Verify the outcome
  EXPECT_FALSE(foundexisting);
}

TEST_F(SetROITest, addDataCreateROI_contained)
{
  // Arrange: bring SUT to desired state
  bed12 query("Example2", 1000, 2500, "test");

  // Act: call methods on SUT, capture output
  // Example2 is 1000 to 5000
  sut.addDataCreateROI("Example2", 1010, 1020, -3);        // The data (bedGraph)

  results = sut.findOverlapIntervals(&query);

  // Assert: Verify the outcome, assumes only one result!
  // Contained won't change boundaries
  EXPECT_EQ((double)1000, results[0]->start);
}

TEST_F(SetROITest, addDataCreateROI_newIdentifier)
{
  // Act: call methods on SUT, capture output
  sut.addDataCreateROI("NewChr", 10, 20, -3);        // The data (bedGraph)

  // Assert: Verify the outcome
  // Fixture has 2 (Example1, Example2) this should add a third (Chr1)
  EXPECT_EQ(sut.chr_names.num_elements,3);
  // And identifier should exist in bimap
  EXPECT_TRUE(sut.chr_names.lookupIndex("NewChr") >= 0);
}

TEST_F(SetROITest, addDataCreateROI_edgecase)
{
  // Arrange: bring SUT to desired state
  bed12 query("Example2", 1000, 2500, "test");

  // Act: call methods on SUT, capture output
  // Example2 is 1000 to 5000
  sut.addDataCreateROI("Example2", 950, 1020, -3);        // The data (bedGraph)

  results = sut.findOverlapIntervals(&query);

  // Assert: Verify the outcome, assumes only one result!
  // Should expand boundary to include data added. 
  EXPECT_EQ((double)950, results[0]->start);
}
