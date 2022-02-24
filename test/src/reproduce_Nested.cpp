/**
 * @file reproduce_Nested.cpp 
 * @author Robin Dowell
 * @brief Recreate the nested regions bug 
 * @date 2022-02-01
 * 
 */
#include "gmock/gmock.h"
#include "load.h"  // contains segment_fits class

TEST(nested, reproduceNested)
{
  // Arrange: bring SUT to desired state
  string file_name = "../testNested.K-models.tsv"; 
  params * P  = new params();
  vector<segment_fits *> fits 		= load::load_K_models_out(file_name);

  // Act: call methods on SUT, capture output
  load::write_out_bidirectionals_with_penalty(fits, P, 2, 1);

  // Assert: Verify the outcome
  EXPECT_EQ(fits[0]->model, 5);
}
