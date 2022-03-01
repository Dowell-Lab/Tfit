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
    // std::cout << ran_num << std::endl;
    
    // Assert
    ASSERT_LE(abs(ran_num - 0.4610), 0.0001);  // Using a "tolerance" of 0.0001
}

TEST(Random, fetchNormalRandom)
{
    // Arrange
    Random sut;
	std::normal_distribution<double> ndist(8, 2);

    // Act
    double ran_num = sut.fetchNormal(8,2);
    // std::cout << ran_num << std::endl;

    // Assert
    ASSERT_LE(abs(ran_num - 11.1246), 0.0001);  // Using a "tolerance" of 0.0001
}