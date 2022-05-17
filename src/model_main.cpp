/**
 * @file model_main.cpp
 * @author Joey Azofeifa
 * @brief This is the primary executable for running the model module. 
 * @version 0.1
 * @date 2016-05-20
 * 
 */
#include "model_main.h"

#include <omp.h>

#include "density_profiler.h"
#include "MPI_comm.h"
#include "template_matching.h"

using namespace std;

/**
 * @brief The "model" module code.
 * Runs the EM model on a collection of regions.  
 * @param P		The parameters with which the program was run.
 * @param rank	MPI processor number
 * @param nprocs	MPI total processes available
 * @param density	(when called by main, this is zero)
 * @param job_ID	MPI job identifier
 * @param LG		a log file associated with this process
 * @return 	int  Is this a status code?
 */
int model_run(params * P, int rank, int nprocs, double density, int job_ID, Log_File * LG){
	int verbose 	= stoi(P->p["-v"]);
	LG->write("\ninitializing model module...............................done\n\n",verbose);
	int threads 	= P->threads; // omp_get_max_threads();//number of OpenMP threads that are available for use	
	string job_name = P->p["-N"];

	//=======================================================================================
	//input file paths
	string forward_bed_graph_file 	= P->p["-i"];
	string reverse_bed_graph_file 	= P->p["-j"];
	string joint_bed_graph_file 	= P->p["-ij"];
	string interval_file 			= P->p["-k"];
	string out_file_dir 			= P->p["-o"];
	string spec_chrom 				= P->p["-chr"];

	//=======================================================================================
	//(1a) load intervals and keep track of their associated IDS
	map<int, string> IDS;
	vector<segment *> FSI;
	LG->write("loading intervals of interest...........................",verbose);
	FSI 	= load::load_intervals_of_interest(interval_file, IDS, P,0 );
	if (FSI.empty()){
		if (rank==0){
			printf("exiting...\n");
		}
		return 1;
	}
	LG->write("done\n",verbose);
	//(1b) now broadcast the intervals of interest to individual MPI processes
	LG->write("sending interval assignments............................",verbose);
	map<string, vector<segment *> > GG 	= MPI_comm::send_out_single_fit_assignments(FSI, rank, nprocs);
	LG->write("done\n",verbose);

	//=======================================================================================
	//(2a) load bedgraph files and insert them into intervals of interest (interval tree...)
	LG->write("inserting bedgraph data.................................",verbose);
	vector<segment*> integrated_segments= load::insert_bedgraph_to_segment_joint(GG, 
		forward_bed_graph_file, reverse_bed_graph_file, joint_bed_graph_file, rank);

    std::cout << "Loaded: " << std::endl;
    vector<segment*>::iterator it;
	for(it=integrated_segments.begin(); it!=integrated_segments.end(); it++) {
       std::cout << (*it)->write_withData() << std::endl;
	}

	//(2b) for each segment we are going to bin and scale and center, numerical stability
	LG->write("done\n",verbose);
	LG->write("binning, centering, scaling.............................",verbose);
	load::BIN(integrated_segments, stod(P->p["-br"]), stod(P->p["-ns"]),true);	
	LG->write("done\n",verbose);

	for(it=integrated_segments.begin(); it!=integrated_segments.end(); it++) {
       std::cout << (*it)->write_withData() << std::endl;
	}

	//=======================================================================================
	//(3a) now run template matching for seeding the EM  
	LG->write("running template matching...............................",verbose);
	slice_ratio SC;
  // WHY are these hard coded here?
	SC.mean = 0.78, SC.std = 0.08; //this dependent on -w 0.9 !!!
	SC.set_2(stod(P->p["-bct"]));
	
	run_global_template_matching(integrated_segments, P, SC);	
    std::cout << SC.write_contents() << std::endl;
	for(it=integrated_segments.begin(); it!=integrated_segments.end(); it++) {
       std::cout << (*it)->write_bidirectional_bounds() << std::endl;
	}
	LG->write("done\n",verbose);

	//=======================================================================================
	//(4a) now going to run the model across all segments
	vector<map<int, vector<simple_c_free_mode> >> FITS 		= run_model_across_free_mode(integrated_segments,
		 P,LG);
	//(4b) gather all the model fits
	LG->write("gathering all model fits................................",verbose);
	map<int, map<int, vector<simple_c_free_mode>  > > GGG 	= MPI_comm::gather_all_simple_c_free_mode(FITS, rank, 
		nprocs);
	LG->write("done\n",verbose);
	//(4c) now write out model fits
	if (rank==0){//write_out_to_MLE, //out_dir+  P->p["-N"] + "-" + to_string(job_ID)+  "_K_models_MLE.tsv"
		LG->write("writing out results (MLE)...............................",verbose);
		string file_name = "";
		load::write_out_models_from_free_mode(GGG, P, job_ID, IDS, density, file_name);
		LG->write("done\n",verbose);
		LG->write("loading results (MLE)...................................",verbose);
		vector<segment_fits *> fits 		= load::load_K_models_out(file_name);
		LG->write("done\n",verbose);		
		LG->write("writing out results (model selection)...................",verbose);
		load::write_out_bidirectionals_with_penalty(fits, P, job_ID, density );
		LG->write("done\n",verbose);
	}
	LG->write("\nexiting model module....................................done\n\n",verbose);
	//wait for everybody to catch up
	MPI_comm::wait_on_root(rank, nprocs);

	return 1;	

}
