/**
 * @file test_Random.cpp
 * @author Robin Dowell
 * @brief Unit Testing: testing Tfit/src/helper.cpp Random object
 * @version 0.1
 * @date 2022-02-28
 * 
 */
#include "gmock/gmock.h"
#include "helper.h"

TEST(Random, fetchProbabilityRandom)
{
    // Arrange
    Random sut;

    // Act
    double ran_num = sut.fetchProbability();

    // Assert
    EXPECT_EQ(sut.testing, 0);
    ASSERT_GE(ran_num, 0);
    ASSERT_LE(ran_num, 1);
}

TEST(Random, fetchProbabilitySeed)
{
    // Arrange
    Random sut(1974);

    // Act
    double ran_num = sut.fetchProbability();

    // Assert
    EXPECT_EQ(sut.testing, 1);
    ASSERT_GE(ran_num, 0);
    ASSERT_LE(ran_num, 1);
}

TEST(Random, fetchNormalRandom)
{
    // Arrange
    Random sut;
	std::normal_distribution<double> ndist(8, 2);

    // Act
    double ran_num = sut.fetchNormal(8,2);

    // Assert
    EXPECT_EQ(sut.testing, 0);
    ASSERT_GE(ran_num, ndist.min());
    ASSERT_LE(ran_num, ndist.max());
}

TEST(Random, fetchNormalSeed)
{
    // Arrange
    Random sut(1974);

    // Act
    double ran_num = sut.fetchNormal(8,2);

    // Assert
    EXPECT_EQ(sut.testing, 1);
    ASSERT_DOUBLE_EQ(ran_num, 7.3837135785183463);
}