/**
 * @file test_EMseeds.cpp 
 * @author Robin Dowell
 * @brief Testing the distribution classes
 * @date 2022-06-16
 * 
 */
#include "gmock/gmock.h"
#include "EMseeds.h"


class SeedsTest: public :: testing::Test {
    protected:
        void SetUp() override {
            PointCov seed1((double)3., (double)13.);
            PointCov seed2((double)18., (double)23.);
            PointCov seed3((double)25., (double)3.);
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
    EXPECT_EQ(max, 23.);
}

TEST_F(SeedsTest, getMinSeed)
{
    // Act: call methods on SUT, capture output
    double min = sut.getMinWeight();

    // Assert: Verify the outcome
    EXPECT_EQ(min, 3.);
}

TEST_F(SeedsTest, sortByWeight)
{
    // Act: call methods on SUT, capture output
    sut.SortByWeights();

    // Assert: Verify the outcome
    EXPECT_EQ(sut.mu_seeds[0].coverage, 3.);;
}

TEST_F(SeedsTest, sortByCoordinate)
{
    // Act: call methods on SUT, capture output
    sut.SortByPositions();

    // Assert: Verify the outcome
    EXPECT_EQ(sut.mu_seeds[0].coordinate, 3.);;
}


