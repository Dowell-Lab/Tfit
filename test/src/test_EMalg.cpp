/**
 * @file test_EMalg.cpp 
 * @author Robin Dowell
 * @brief Testing the Bedfile and Bedgraph classes: \ref Bedfile Bedgraph
 * @date 2022-07-26
 * 
 */
#include "gmock/gmock.h"
#include "EMalg.h"

TEST(EMalg, fitFailswithNoData) 
{
    // Arrange: bring SUT to desired state
    EMalg sut;

    // Act: call methods on SUT, capture output
    bool success = sut.fit();

    // Assert: Verify the outcome
    EXPECT_FALSE(success);
}

// TEST(EMalg, NoiseOnlyTest)

// TEST(EMalg, fit)

// Test(EMalg, computeBackgroundModel)

// Test(EMalg, Initialize)

// TEST(EMalg, Estep)

// TEST(EMalg, Mstep)

// TEST(EMalg, adjustBounds)
