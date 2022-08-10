/**
 * @file reproduce_EMcomponents.cpp 
 * @author Robin Dowell
 * @brief Checking that revised subroutines reproduce the portions of 
 * the original EM code that are designed to replace.
 * @date 2022-07-28
 * 
 */
#include "gmock/gmock.h"
#include "load.h"  // contains segment_fits class
#include "model.h" // contains fit3 implementation (modifiable version of fit2)

#include "Bedfile.h"
#include "EMalg.h"

/*
TEST(reproduceEM, dataloading)
{
  // Lets load a typical small bedgraph file.
  string joint_bedgraph = "../examples/typical_region.bg";    //chr21 33401693 33407411	

  // Joey's code
	map<string, int> chrom_to_ID;
	map<int, string> ID_to_chrom;
  string emptyfilename;
  cout << "Setup Variables" << std::endl;
	vector<segment *> 	segments 	= load::load_bedgraphs_total(emptyfilename, 
			emptyfilename, joint_bedgraph, 25, 100, "chr21", chrom_to_ID, ID_to_chrom );  
  cout << "Load Joey" << std::endl;
  segment *data = segments[0];  // We just want the first one
  classifier sut = classifier(0, 0.0001, 2000, 0.05, 0, 1, 1, 1, 1, 1, 1, 0);
  double center = 21;
  data->centers.push_back(center);
  cout << "Pick Segment and Centers" << std::endl;

  Bedgraph genData;
  genData.load_file(joint_bedgraph,false);
  cout << "New Method Load Data" << std::endl;
  EMalg v2model;

  // Act on sut (run EM!)
  // cout << "Before: " +  sut.write_classifier_status() << std::endl;
  // int returned = sut.fit3(data, data->centers, 0, 0);    
  // cout << "components: " + sut.write_components() << std::endl;
  // cout << "After: " +  sut.write_classifier_status() << std::endl;

  // Assert: Verify the outcome
  EXPECT_LE(abs(sut.ll + 218769.31), 0.01);  
}
*/
