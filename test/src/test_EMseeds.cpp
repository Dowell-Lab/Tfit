/**
 * @file test_EMseeds.cpp 
 * @author Robin Dowell
 * @brief Testing the distribution classes
 * @date 2022-06-16
 * 
 */
#include "gmock/gmock.h"
#include "split.h"
#include "EMseeds.h"

TEST(Seeds, writeEmpty) {
  // Arrange
  Seeds sut;
  // Act
  std::string output = sut.writeSeedsAsBedFields();
  //Assert
  EXPECT_EQ(output, "0.\t\t");

}

TEST(Seeds, EmptySize) {
  // Arrange
  Seeds sut;
  // Act
  //Assert
  EXPECT_EQ(sut.getNumSeeds(), 0);
}

class SeedsTest: public :: testing::Test {
    protected:
        void SetUp() override {
            PointCov seed1((double)3., (double).13);
            PointCov seed2((double)18., (double).23);
            PointCov seed3((double)25., (double).03);
            sut.mu_seeds.push_back(seed1);
            sut.mu_seeds.push_back(seed2);
            sut.mu_seeds.push_back(seed3);
        }
    // void TearDown() override { }
    Seeds sut;
};

TEST_F(SeedsTest, getNumSeeds)
{
    // Act: call methods on SUT, capture output
    int numS = sut.getNumSeeds();

    // Assert: Verify the outcome
    EXPECT_EQ(numS, 3);
}

TEST_F(SeedsTest, getMaxSeed)
{
    // Act: call methods on SUT, capture output
    double max = sut.getMaxWeight();

    // Assert: Verify the outcome
    EXPECT_EQ(max, 0.23);
}

TEST_F(SeedsTest, getMinSeed)
{
    // Act: call methods on SUT, capture output
    double min = sut.getMinWeight();

    // Assert: Verify the outcome
    EXPECT_EQ(min, .03);
}

TEST_F(SeedsTest, sortByWeight)
{
    // Act: call methods on SUT, capture output
    sut.SortByWeights();

    // Assert: Verify the outcome
    EXPECT_EQ(sut.mu_seeds[0].coverage, .03);;
}

TEST_F(SeedsTest, sortByCoordinate)
{
    // Act: call methods on SUT, capture output
    sut.SortByPositions();

    // Assert: Verify the outcome
    EXPECT_EQ(sut.mu_seeds[0].coordinate, 3.);;
}

TEST_F(SeedsTest, writeAsBedFields)
{
    // Act: call methods on SUT, capture output
    std::string output = sut.writeSeedsAsBedFields();

    // Assert: Verify the outcome
    EXPECT_EQ(output, "3.\t13.,23.,3.\t3.,18.,25.");

}

TEST_F(SeedsTest, IOonSeeds)
{
    std::string output = sut.writeSeedsAsBedFields();

    std::vector<std::string> lineArray; // Contents of line, split on tab (\t)
    lineArray = string_split(output, '\t');
    std::string seedNum = lineArray[0];
    std::string seedWeights = lineArray[1];
    std::string seedPos = lineArray[2];

    Seeds newSUT;

    // Act: call methods on SUT, capture output
    newSUT.getSeedsfromBedFields(seedNum, seedWeights, seedPos);

    // Assert: Verify the outcome
    EXPECT_EQ(newSUT.getNumSeeds(), 3);
    EXPECT_EQ(newSUT.getMaxWeight(), 0.23);
}







