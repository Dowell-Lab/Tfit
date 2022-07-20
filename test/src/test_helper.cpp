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

// Bimap: addIdentifier, lookupIndex, lookupName
TEST(Bimap, addIdentifierNew)
{
    // Arrange
    Bimap sut;

    // Act
    sut.addIdentifier("chromosome1");
    
    // Assert
   EXPECT_EQ(sut.num_elements, 1);
}

TEST(Bimap, addIdentifierExisting)
{
    // Arrange
    Bimap sut;
    sut.addIdentifier("chromosome1");

    // Act
    sut.addIdentifier("chromosome1");
    
    // Assert
   EXPECT_EQ(sut.num_elements, 1);
}

TEST(Bimap, lookupIndex)
{
    // Arrange
    Bimap sut;
    sut.addIdentifier("chromosome1");
    sut.addIdentifier("chr2");

    // Act
    int index = sut.lookupIndex("chr2");
    
    // Assert
   EXPECT_EQ(index, 1);
}

TEST(Bimap, lookupIndexNoExist)
{
    // Arrange
    Bimap sut;
    sut.addIdentifier("chromosome1");
    sut.addIdentifier("chr2");

    // Act
    int index = sut.lookupIndex("chrX");
    
    // Assert
   EXPECT_EQ(index, -1);
}

TEST(Bimap, lookupName)
{
    // Arrange
    Bimap sut;
    sut.addIdentifier("chromosome1");
    sut.addIdentifier("chr2");

    // Act
    std::string name = sut.lookupName(0);
    
    // Assert
   EXPECT_EQ(name, "chromosome1");
}

TEST(Bimap, lookupNameNoExist)
{
    // Arrange
    Bimap sut;
    sut.addIdentifier("chromosome1");
    sut.addIdentifier("chr2");

    // Act
    std::string name = sut.lookupName(3);
    
    // Assert
   EXPECT_EQ(name, "");
}