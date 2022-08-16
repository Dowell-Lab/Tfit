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

TEST_F(SeedsTest, writeAsHalfBed12)
{
    // Act: call methods on SUT, capture output
    std::string output = sut.writeHalfBed12(1, 100);

    // Assert: Verify the outcome
    EXPECT_EQ(output, "\t1.\t100.\t00,00,00\t3.\t13.,23.,3.\t3.,18.,25.");

}

TEST_F(SeedsTest, IOonHalfBed12)
{
    std::string output = sut.writeHalfBed12(1, 100);
    std::string firstSix = "chr1\t1\t100\ttest\t+\t300";
    std::string line = firstSix + output;
    std::vector<std::string> lineArray; // Contents of line, split on tab (\t)
    lineArray = string_split(line, '\t');

    Seeds newSUT;
    // Act: call methods on SUT, capture output
    newSUT.grabSeedsfromBed12(lineArray);

    // Assert: Verify the outcome
    EXPECT_EQ(newSUT.getNumSeeds(), 3);
    EXPECT_EQ(newSUT.getMaxWeight(), 0.23);

}







