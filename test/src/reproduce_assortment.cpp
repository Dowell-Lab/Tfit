/**
 * @file reproduce_assortment.cpp 
 * @author Robin Dowell
 * @brief Assorted reproduce tests for Joey's code.
 * @date 2022-02-01
 * 
 */
#include "gmock/gmock.h"
#include "load.h"  // contains segment class

TEST(assorted, SegmentLengthEquivalent)
{
    // Arrange: bring SUT to desired state
    segment sut = segment("chrTest", 100, 1000, 3, "+"); 

    // Act: call methods on SUT, capture output
    double length = sut.getXLength(); 
    double oldLength = sut.maxX - sut.minX;

    // Assert: Verify the outcome
    EXPECT_EQ(length, oldLength);
}

TEST(assorted, reproduceNestedBug)
{
  // Arrange: bring SUT to desired state
  string file_name = "../examples/testNested.K-models.tsv"; 
  params * P  = new params();
  vector<segment_fits *> fits 		= load::load_K_models_out(file_name);

  // Act: call methods on SUT, capture output
  load::write_out_bidirectionals_with_penalty(fits, P, 2, 1);

  // Assert: Verify the outcome
  EXPECT_EQ(fits[0]->model, 5);
  delete(P);
}
