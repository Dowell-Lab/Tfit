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
  // Lets load a simple relatively small bedgraph file.
  string emptyfilename;
  string joint_bedgraph = "../examples/reduced.bg";
	map<string, int> chrom_to_ID;
	map<int, string> ID_to_chrom;
	vector<segment *> 	segments 	= load::load_bedgraphs_total(emptyfilename, 
			emptyfilename, joint_bedgraph, 25, 100, "chr21", chrom_to_ID, ID_to_chrom );
  segment *data = segments[0];

  // Populate data centers?
  for (int b = 0; b < data->bidirectional_bounds.size(); b++) {
    double center = data->bidirectional_bounds[b][0] + data->bidirectional_bounds[b][1];
    center /= 2.;
    center -= data->start;
    center /= 100;    // scale
    data->centers.push_back(center);
  }

  // How is the classifier setup before fit2 is run?
	map<int, vector<classifier> > A;
  //expection of binomal(0.05, n) where 0.05 probability mapping noise
	double noise   	= abs(data->stop - data->start)*0.05/(data->fN + data->rN); 
  // Using defaults from parameters
	A[0].push_back( classifier(0, 0.0001, 2000, noise, 0, 1, 1, 1, 1, 1, 1, 0));
  for (int r = 0; r < 10; r++) {
    A[1].push_back(classifier(1, 0.0001, 2000, 0.05, 0, 1, 1, 1, 1, 1, 1, 0));
  }

  typedef map<int, vector<classifier> > ::iterator it_type;
  for (it_type k = A.begin(); k!= A.end(); k++){
			int N 	=  k->second.size();    // I think this is K
			for (int r = 0; r < N; r++ ){   // This is rounds
				A[k->first][r].fit2(data, data->centers,0,0);
		}
	}


  // Act: call methods on SUT, capture output
  // sut.fit2(data,seeds,0,0);  // I believe seeds and two alg modifiers can all be empty/zero

  // Assert: Verify the outcome
  //something about the fits.
}

