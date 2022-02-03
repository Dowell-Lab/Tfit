/**
 * @file test_segment.cpp 
 * @author Robin Dowell
 * @brief Testing the old segment class (to be refactored) 
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
  load::write_out_bidirectionals_ms_pen(fits, P, 2, 1);

  // Act: call methods on SUT, capture output

  // Assert: Verify the outcome
  EXPECT_EQ(fits[0]->model, 4);
}
