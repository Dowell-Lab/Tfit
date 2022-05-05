/**
 * @file test_Bimap.cpp
 * @author Robin Dowell
 * @brief Unit Testing: testing the \ref Bimap object (in helper.h)
 * @version 0.1
 * @date 2022-04-07
 * 
 */
#include "gmock/gmock.h"
#include "helper.h"

TEST(BiMap, createAndSearch)
{
    // Arrange
    Bimap sut;

    // Act
    
    // Assert
   EXPECT_EQ(sut.num_elements, 0);
}

