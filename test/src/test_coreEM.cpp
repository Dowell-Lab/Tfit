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


/* This is a very heavy "test" -- more like trying to peel away most of Tfit's
cruft to get to the core Em algorithm (fit2).  */
TEST(classifier, exerciseFit2)
{
  // Lets load a simple relatively small bedgraph file.
  string emptyfilename;
  string joint_bedgraph = "../examples/reduced.bg";
  // loading sets up these mappings
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

  // Now we have to setup the classifier.  Unclear we need all 
  // this cruft, but doing it anyway to start. 
	map<int, vector<classifier> > A;
  //expection of binomal(0.05, n) where 0.05 probability mapping noise
	double noise   	= abs(data->stop - data->start)*0.05/(data->fN + data->rN); 
  // Using defaults from parameters
  // A[0] is background noise/uniform?
	A[0].push_back( classifier(0, 0.0001, 2000, noise, 0, 1, 1, 1, 1, 1, 1, 0));
  for (int r = 0; r < 10; r++) {
    A[1].push_back(classifier(1, 0.0001, 2000, 0.05, 0, 1, 1, 1, 1, 1, 1, 0));
  }

  /* Unclear that all this iteration is needed.*/
  typedef map<int, vector<classifier> > ::iterator it_type;
  for (it_type k = A.begin(); k!= A.end(); k++){
			int N 	=  k->second.size();    // I think this is K
			for (int r = 0; r < N; r++ ){   // This is rounds
        cout << "Before: " +  A[k->first][r].write_classifier_status() << std::endl;
        // Act: call methods on SUT, capture output
        // So we run this and r rounds per K
				A[k->first][r].fit2(data, data->centers,0,0);
        cout << "After: " +  A[k->first][r].write_classifier_status() << std::endl;
		}
	}

  // Assert: Verify the outcome
  EXPECT_EQ(A[0][0].K, 1);
}

