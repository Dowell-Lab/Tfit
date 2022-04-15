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

TEST(Load, loadBedGraph)
{
  // Lets load a simple relatively small bedgraph file.
  segment *data;    // What needs to be set in here?
  string emptyfilename;
  string joint_bedgraph = "../examples/reduced.bg";
	map<string, int> chrom_to_ID;
	map<int, string> ID_to_chrom;

  // Act: call methods on SUT, capture output
	vector<segment *> 	segments 	= load::load_bedgraphs_total(emptyfilename, 
			emptyfilename, joint_bedgraph, 25, 100, "chr21", chrom_to_ID, ID_to_chrom );

  // Assert: Verify the outcome
  // cout << segments[0]->write_withData() << std::endl;
  EXPECT_EQ(segments[0]->write_interval(),"#chr21:33400500-33409597,56761");

}

