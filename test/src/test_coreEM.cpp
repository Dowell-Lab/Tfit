/**
 * @file test_coreEM.cpp 
 * @author Robin Dowell
 * @brief Exploring the core EM algorithm
 * @date 2022-02-24
 * 
 */
#include "gmock/gmock.h"
#include "load.h"  // contains segment_fits class
#include "model.h" // contains fit2 implementation (classifier object)

TEST(classifier, exerciseFit2)
{
  // Arrange: bring SUT to desired state
  segment *data;
  vector<double> seeds;  
  classifier sut; // fit2 is a subroutine of a classifier

  // Act: call methods on SUT, capture output
  sut.fit2(data,seeds,0,0);

  // Assert: Verify the outcome
  //something about the fits.
}

