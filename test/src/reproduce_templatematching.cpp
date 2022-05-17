/**
 * @file reproduce_templatematching.cpp 
 * @author Robin Dowell
 * @brief recreate the template matching algorithm
 * @date 2022-05-17
 * 
 */
#include "gmock/gmock.h"

#include "load.h"  // contains segment_fits class
#include "model.h" // contains fit2 implementation (classifier object)
#include "FDR.h" // slice_ratio
#include "template_matching.h" // run_global_template_matching

TEST(Bidir, templatematch)
{
  // Lets load a simple relatively small bedgraph file.
  segment *data;    // What needs to be set in here?
  string emptyfilename;
  string joint_bedgraph = "../examples/reduced.bg";
	map<string, int> chrom_to_ID;
	map<int, string> ID_to_chrom;
  params * P  = new params();

	vector<segment *> 	segments 	= load::load_bedgraphs_total(emptyfilename, 
			emptyfilename, joint_bedgraph, 25, 100, "chr21", chrom_to_ID, ID_to_chrom );

	slice_ratio SC;
	SC.mean = 0.78, SC.std = 0.08; //this dependent on -w 0.9 !!!
	SC.set_2(0.95);

  // Act: call methods on SUT, capture output
	double threshold 	= run_global_template_matching(segments, P, SC);	

  // Assert: Verify the outcome
  // cout << segments[0]->write_withData() << std::endl;
  EXPECT_EQ(segments[0]->write_interval(),"#chr21:33400500-33409597,56761");
  EXPECT_EQ(threshold, 1);

}

