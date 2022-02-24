/**
 * @file test_main.cpp
 * @author Robin Dowell 
 * @copyright Copyright 2022 Dowell Lab 
 * @brief This is for testing and debugging only.
 * @version 0.1
 * @date 2022-02-21
 * 
 */

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

#include <chrono>  // This is an unapproved C++11 header
#include <iostream>
#include <limits>
#include <map>
#include <thread>  // This is an unapproved C++11 header

#include "bidir_main.h"
#include "bootstrap.h"
#include "density_profiler.h"
#include "error_stdo_logging.h"
#include "load.h"
#include "MPI_comm.h"
#include "model.h"
#include "model_main.h"
#include "model_selection.h"
#include "read_in_parameters.h"
#include "select_main.h"
#include "template_matching.h"
#include "across_segments.h"
#include "model_single.h"

void setupParams (params *P);
vector<map<int, vector<simple_c_free_mode> >> run_model(vector<segment *> FSI, params * P);

/**
 * @brief Program main.  
 * @param argc
 * @param argv
 * @return 
 */
int main(int argc, char* argv[]) {
  params *P = new params();
  setupParams(P);

  std::string data_file = "../end2end/reduced.bg";
  std::string interval_file = "../end2end/singleregion.bed";

  /* Read in a segment of data */
  map<int, string> IDS;
  vector<segment *> FSI;
  map<string, vector<segment *>> GG;

  FSI = load::load_intervals_of_interest(interval_file, IDS, P, 0);
  for (auto &element : FSI) {
	  cout << element->write_interval() << std::endl;
	  GG[element->chrom].push_back(element);
  }
  vector<segment *> integrated_segments = load::insert_bedgraph_to_segment_joint(GG,
				"", "", data_file, 0);
  /* Centering and scaling */
  load::BIN(integrated_segments, 25, 100, true); // Note 25 and 100 are defaults from params
  for (auto &element : integrated_segments) {
	  cout << element->write_withData() << std::endl;
  }

  /* Attempt to fit a single model (K=1) */
  vector<map<int, vector<simple_c_free_mode>>> FITS = run_model(integrated_segments, P);
}

/* Setup as needed for this run */
void setupParams (params *P) {
  P->threads = 1;
  P->model = 1;
  P->p["-chr"] = "chr21";
  P->p["-pad"] = "2000";
  P->p["-minK"] = "1";
  P->p["-maxK"] = "1";
  P->p["-rounds"] = "5";

}

vector<map<int, vector<simple_c_free_mode> >> run_model (vector<segment *> FSI, params * P) {
	vector<map<int, vector<simple_c_free_mode> >> D;
	typedef map<int, vector<classifier> > ::iterator it_type;

  /* defaults */
	double scale 	= 100; int num_proc 				= 1;
	int verbose 	= 1; double N 		= FSI.size();
	double percent 	= 0; 

	//printf("FSI.size: %d\n", FSI.size());
	for (int i = 0 ; i < FSI.size(); i++){
		if ((i / N) > (percent+0.05)){
			percent 	= (i / N);
		}

		//first need to populate data->centers
		for (int b = 0 ; b < FSI[i]->bidirectional_bounds.size(); b++){
			double center = FSI[i]->bidirectional_bounds[b][0] +  FSI[i]->bidirectional_bounds[b][1] ;
			center/=2.;
			center-=FSI[i]->start;
			center/=scale;
			cout << center << std::endl;
			FSI[i]->centers.push_back(center);
		}

		segment * data 	= FSI[i];
		map<int, vector<classifier> > A 	= make_classifier_struct_free_model(P, FSI[i]);
		for (it_type k = A.begin(); k!= A.end(); k++){
			int N 	=  k->second.size(); // second ... num classifiers
			for (int r = 0; r < N; r++ ){	 // This is "-rounds"
				cout << "RUN EM r:" << r << std::endl;
	  			cout << data->write_centers() << std::endl;
				A[k->first][r].fit2(data, data->centers,0,0);
			}
		}
		D.push_back(get_max_from_free_mode(A, FSI[i], i));
	}
	return D;
}
