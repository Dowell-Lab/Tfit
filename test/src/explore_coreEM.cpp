/**
 * @file explore_coreEM.cpp 
 * @author Robin Dowell
 * @brief Exploring the core EM algorithm, not a test so much as a framework 
 * for debugging.
 * @date 2022-02-24
 * 
 */
#include "gmock/gmock.h"
#include "load.h"  // contains segment_fits class
#include "model.h" // contains fit2 implementation (classifier object)

/* This is a very heavy "test" -- more like trying to peel away most of Tfit's
cruft to get to the core Em alg2rithm (fit2).  */
TEST(classifier, exerciseFit2)
{
  string emptyfilename;

  // Lets load a typical small bedgraph file.
  // string joint_bedgraph = "../examples/reduced.bg";    //chr21 	
  string joint_bedgraph = "../examples/typical_region.bg";    //chr21 33401693 33407411	

  // loading sets up these mappings, we are going to ignore them.
	map<string, int> chrom_to_ID;
	map<int, string> ID_to_chrom;

  // The assorted hard coded values are set by defaults in read_in_parameters

	vector<segment *> 	segments 	= load::load_bedgraphs_total(emptyfilename, 
			emptyfilename, joint_bedgraph, 25, 100, "chr21", chrom_to_ID, ID_to_chrom );  
  segment *data = segments[0];  // We just want the first one

  // std::cout << data->write_Data_Columns() << std::endl;

  classifier sut = classifier(2, 0.0001, 2000, 0.05, 0, 1, 1, 1, 1, 1, 1, 0);

  double center = 20;
  data->centers.push_back(center);
  center = 45;
  data->centers.push_back(center);

  // Act on sut (run EM!)
  cout << "Before: " +  sut.write_classifier_status() << std::endl;
  sut.fit2(data, data->centers, 0, 0);    
  cout << "components: " + sut.write_components() << std::endl;
  cout << "After: " +  sut.write_classifier_status() << std::endl;

  // Assert: Verify the outcome
  EXPECT_EQ(sut.K, 2);
  // EXPECT_EQ(sut.ll, -10);  Why is this -inf here?
}

