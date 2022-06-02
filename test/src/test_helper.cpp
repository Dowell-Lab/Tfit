/**
 * @file test_helper.cpp
 * @author Robin Dowell
 * @brief Unit Testing: Helper functions
 * @version 0.1
 * @date 2022-02-28
 * 
 */
#include "gmock/gmock.h"
#include "helper.h"

TEST(helper, Random_fetchProbabilityRandom)
{
    // Arrange
    Random sut;

    // Act
    double ran_num = sut.fetchProbability();
    // std::cout << ran_num << std::endl;
    
    // Assert
    ASSERT_LE(abs(ran_num - 0.4610), 0.0001);  // Using a "tolerance" of 0.0001
}

TEST(helper, Random_fetchNormalRandom)
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

TEST(helper, Bimap_createAndSearch)
{
    // Arrange
    Bimap sut;

    // Act
    
    // Assert
   EXPECT_EQ(sut.num_elements, 0);
}