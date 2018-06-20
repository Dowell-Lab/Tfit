#include "bidir_main.h"
#include "template_matching.h"
#include "density_profiler.h"
#include "MPI_comm.h"
#include <mpi.h>
#include <omp.h>
#include "model_main.h"
#include "select_main.h"
#include "error_stdo_logging.h"
#include "FDR.h"
#include "BIC.h"
#include "old_template_matching.hpp"
using namespace std;

/** Performs bidirectional predictions given an old style parameter structure; use bidir_run_pwrapper instead.
 * @depreciated
 * @param P params struct containing parsed command line arguments.
 * @param rank MPI rank parameter. This is used to determine the level of output a given instance of this function will generate.
 * @param nprocs Number of processors available. This should, in practice, be an MPI runtime parameter.
 * @param job_ID Current job identifier. This should, in practice, be an MPI runtime parameter.
 * @param LG Log_File object pointing to a valid output log file.
 * @return An output status code. 0 indicates no errors during operation. Any other value indicates errors.
 */
int bidir_run(params * P, int rank, int nprocs, int job_ID, Log_File * LG){

	int verbose 	= stoi(P->p["-v"]);
	P->p["-merge"] 	= "1";
	
	LG->write("\ninitializing bidir module...............................done\n", verbose);
	int threads 	= omp_get_max_threads();//number of OpenMP threads that are available for use
	
	//===========================================================================
	//get job_ID and open file handle for log files
	string job_name = P->p["-N"];
	
	//===========================================================================
	//input files and output directories
	string forward_bedgraph 	= P->p["-i"]; //forward strand bedgraph file
	string reverse_bedgraph 	= P->p["-j"]; //reverse strand bedgraph file
	string joint_bedgraph 		= P->p["-ij"]; //joint forward and reverse strand bedgraph file
	string tss_file 		= P->p["-tss"];
	string out_file_dir 		= P->p["-o"] ;//out file directory
	//===========================================================================
	//template searching parameters
	double sigma, lambda, foot_print, pi, w;
	double ns 	= stod(P->p["-ns"]);
	sigma 		= stod(P->p["-sigma"]), lambda= stod(P->p["-lambda"]);
	foot_print 	= stod(P->p["-foot_print"]), pi= stod(P->p["-pi"]), w= stod(P->p["-w"]);

	//(2a) read in bedgraph files 
	map<string, int> chrom_to_ID;
	map<int, string> ID_to_chrom;
       
	vector<double> parameters 	= {sigma, lambda, foot_print,pi, w};
	if (not tss_file.empty() and rank == 0){
		vector<segment *> FSI;
		LG->write("loading TSS intervals...................................",verbose);
		map<int, string> IDS;
		vector<segment *> tss_intervals 	= load::load_intervals_of_interest(tss_file, 
			IDS, P, 1 );
		map<string, vector<segment *>> GG 	= MPI_comm::convert_segment_vector(tss_intervals);
		LG->write("done\n", verbose);
		
		LG->write("inserting coverage data.................................",verbose);
		vector<segment*> integrated_segments= load::insert_bedgraph_to_segment_joint(GG, 
			forward_bedgraph, reverse_bedgraph, joint_bedgraph, rank);
		LG->write("done\n", verbose);

		LG->write("Binning/Normalizing TSS intervals.......................",verbose);
		load::BIN(integrated_segments, stod(P->p["-br"]), stod(P->p["-ns"]),true);	
		LG->write("done\n", verbose);
		
		LG->write("Computing Average Model.................................",verbose);
		parameters  	= compute_average_model(integrated_segments, P);
		//these values are on the genome scale (non normalize to ns)
		LG->write("done\n", verbose);
		LG->write("\nAverage Model Parameters\n", verbose);
		LG->write("-sigma      : " + to_string(parameters[0])+ "\n", verbose);
		LG->write("-lambda     : " + to_string(parameters[1])+ "\n", verbose);
		LG->write("-foot_print : " + to_string(parameters[2])+ "\n", verbose);
		LG->write("-pi         : " + to_string(parameters[3])+ "\n", verbose);
		LG->write("-w          : " + to_string(parameters[4])+ "\n\n", verbose);
	}
	parameters 				= MPI_comm::send_out_parameters( parameters, rank, nprocs);
	P->p["-sigma"] 	       = to_string(parameters[0]);
	P->p["-lambda"]        = to_string(parameters[1]);
	P->p["-foot_print"]    = to_string(parameters[2]);
	P->p["-pi"] 	       = to_string(parameters[3]);
	P->p["-w"] 	       = to_string(parameters[4]);

	LG->write("loading bedgraph files..................................", verbose);
	vector<segment *> 	segments 	= load::load_bedgraphs_total(forward_bedgraph, reverse_bedgraph, joint_bedgraph,
		stoi(P->p["-br"]), stof(P->p["-ns"]), P->p["-chr"], chrom_to_ID, ID_to_chrom );

	if (segments.empty()){
		printf("exiting...\n");
		return 1;
	}
	LG->write("done\n", verbose);




	slice_ratio SC;
	if (stoi(P->p["-FDR"] ) ){
	  LG->write("getting likelihood score distribution...................", verbose);
	  SC                      = get_slice(segments, pow(10,6) , pow(10,4) ,P  );
	  LG->write("done\n\n", verbose);
	  if (not SC.converged){
	    LG->write("converged            : False (restoring default values)\n"  ,verbose );
	  }else{
	    LG->write("converged            : True\n" ,verbose );
	  }
	  LG->write("score mean           : "+to_string(SC.mean) + "\n" ,verbose );
	  LG->write("standard Deviation   : "+to_string(SC.std ) + "\n" ,verbose );
	  LG->write("h                    : "+to_string(SC.w ) + "\n" ,verbose );
	  LG->write("threshold            : "+to_string(SC.threshold) + "\n\n" ,verbose );
	}
	else{
	  SC.mean = 0.6, SC.std = 0.001; //this dependent on -w 0.9 !!!
	  SC.set_2(stod(P->p["-bct"]));
	}

	//================================================
	//computing this SC via multiple times of the EM
	//run, should probably write something to combind
	//them all; i.e. MPI

	//=================================================
	//(2b) so segments is indexed by inidividual chromosomes, want to broadcast 
	//to sub-processes and have each MPI call run on a subset of segments
	vector<segment*> all_segments  	= segments;
	LG->write("slicing segments........................................", verbose);
	segments 						= MPI_comm::slice_segments(segments, rank, nprocs);	
	LG->write("done\n", verbose);

	//===========================================================================
	//(3a) now going to run the template matching algorithm based on pseudo-
	//moment estimator and compute BIC ratio (basically penalized LLR)
	LG->write("running template matching algorithm.....................", verbose);
	double threshold 	= run_global_template_matching(segments, out_file_dir, P, SC);	
	//(3b) now need to send out, gather and write bidirectional intervals 
	LG->write("done\n", verbose);
	


	LG->write("scattering predictions to other MPI processes...........", verbose);
	int total =  MPI_comm::gather_all_bidir_predicitions(all_segments, 
							     segments , 1, nprocs, out_file_dir, job_name, job_ID,P,0);
	MPI_Barrier(MPI_COMM_WORLD); //make sure everybody is caught up!

	LG->write("done\n", verbose);
	if (rank==0){
	  LG->write("\nThere were " +to_string(total) + " prelimary bidirectional predictions\n\n", verbose);
	}
	
	//===========================================================================
	//this should conclude it all
	LG->write("clearing allocated segment memory.......................", verbose);	
	load::clear_segments(all_segments);
	LG->write("done\n", verbose);
	//===========================================================================
	//(4) if MLE option was provided than need to run the model_main::run()
	//
	if (stoi(P->p["-MLE"])){
		P->p["-k"] 	= P->p["-o"]+ job_name+ "-" + to_string(job_ID)+ "_prelim_bidir_hits.bed";
		model_run(P, rank, nprocs,0, job_ID, LG);
		
	}
	LG->write("exiting bidir module....................................done\n\n", verbose);
	return 1;
}

/** Performs bidirectional predictions given newer style ParamWrapper inputs.
 * @param pw ParamWrapper containing parsed command line arguments.
 * @param rank MPI rank parameter. This is used to determine the level of output a given instance of this function will generate.
 * @param nprocs Number of processors available. This should, in practice, be an MPI runtime parameter.
 * @param job_ID Current job identifier. This should, in practice, be an MPI runtime parameter.
 * @param LG Log_File object pointing to a valid output log file.
 * @return An output status code. 0 indicates no errors during operation. Any other value indicates errors.
 */
int bidir_run_pwrapper(ParamWrapper *pw, int rank, int nprocs, int job_ID, Log_File * LG){


	int verbose 	= pw->verbose;
	//pw->merge 	= "1";
	
	LG->write("\ninitializing bidir module...............................done\n", verbose);
	int threads 	= omp_get_max_threads();//number of OpenMP threads that are available for use
	
	//===========================================================================
	//get job_ID and open file handle for log files
	string job_name = pw->jobName;
	
	//===========================================================================
	//input files and output directories
	string forward_bedgraph 	= pw->forwardStrand; //forward strand bedgraph file
	string reverse_bedgraph 	= pw->reverseStrand; //reverse strand bedgraph file
	string joint_bedgraph 		= pw->mergedStrand; //joint forward and reverse strand bedgraph file
	string tss_file 		= pw->promoterTSS;
	string out_file_dir 		= pw->outputDir;//out file directory
	//===========================================================================
	//template searching parameters
	double sigma, lambda, foot_print, pi, w;
	double ns 	= pw->ns;
	sigma 		= pw->sigma, lambda= pw->lambda;
	foot_print 	= pw->footPrint, pi=pw->pi, w=pw->w;

	//(2a) read in bedgraph files 
	map<string, int> chrom_to_ID;
	map<int, string> ID_to_chrom;
    printf("Start footprint: %f\n", pw->footPrint);
       
	vector<double> parameters 	= {sigma, lambda, foot_print,pi, w};
    printf("Value of not tss_file.empty: %d\n", !tss_file.empty());
    printf("Tss filename: %s\n", tss_file.c_str());
	if (tss_file!="" and rank == 0){
		vector<segment *> FSI;
		LG->write("loading TSS intervals...................................",verbose);
		map<int, string> IDS;
		vector<segment *> tss_intervals 	= load::load_intervals_of_interest_pwrapper(tss_file, 
			IDS, pw, 1 );
		map<string, vector<segment *>> GG 	= MPI_comm::convert_segment_vector(tss_intervals);
		LG->write("done\n", verbose);
		
		LG->write("inserting coverage data.................................",verbose);
		vector<segment*> integrated_segments= load::insert_bedgraph_to_segment_joint(GG, 
			forward_bedgraph, reverse_bedgraph, joint_bedgraph, rank);
		LG->write("done\n", verbose);

		LG->write("Binning/Normalizing TSS intervals.......................",verbose);
		load::BIN(integrated_segments, pw->br, pw->ns,true);	
		LG->write("done\n", verbose);
		
		LG->write("Computing Average Model.................................",verbose);
		parameters  	= compute_average_model_pwrapper(integrated_segments, pw);
		//these values are on the genome scale (non normalize to ns)
		LG->write("done\n", verbose);
		LG->write("\nAverage Model Parameters\n", verbose);
		LG->write("-sigma      : " + to_string(parameters[0])+ "\n", verbose);
		LG->write("-lambda     : " + to_string(parameters[1])+ "\n", verbose);
		LG->write("-foot_print : " + to_string(parameters[2])+ "\n", verbose);
		LG->write("-pi         : " + to_string(parameters[3])+ "\n", verbose);
		LG->write("-w          : " + to_string(parameters[4])+ "\n\n", verbose);
	}
	parameters 				= MPI_comm::send_out_parameters( parameters, rank, nprocs);
    //Todo: consider the implications of this:
    pw->sigma=parameters[0];
    pw->lambda=parameters[1];
    pw->footPrint=parameters[2];
    pw->pi=parameters[3];
    pw->w=parameters[4];
    printf("Stop footprint: %f\n", pw->footPrint);
	//P->p["-sigma"] 	       = to_string(parameters[0]);
	//P->p["-lambda"]        = to_string(parameters[1]);
	//P->p["-foot_print"]    = to_string(parameters[2]);
	//P->p["-pi"] 	       = to_string(parameters[3]);
	//P->p["-w"] 	       = to_string(parameters[4]);

	LG->write("loading bedgraph files..................................", verbose);
	vector<segment *> 	segments 	= load::load_bedgraphs_total(forward_bedgraph, reverse_bedgraph, joint_bedgraph,
		pw->br, pw->ns, pw->chromosome, chrom_to_ID, ID_to_chrom );

	if (segments.empty()){
		printf("exiting...\n");
		return 1;
	}
	LG->write("done\n", verbose);

	slice_ratio SC;
	if (pw->fdr){
	  LG->write("getting likelihood score distribution...................", verbose);
	  SC                      = get_slice_pwrapper(segments, pow(10,6) , pow(10,4) , pw);
	  LG->write("done\n\n", verbose);
	  if (not SC.converged){
	    LG->write("converged            : False (restoring default values)\n"  ,verbose );
	  }else{
	    LG->write("converged            : True\n" ,verbose );
	  }
	  LG->write("score mean           : "+to_string(SC.mean) + "\n" ,verbose );
	  LG->write("standard Deviation   : "+to_string(SC.std ) + "\n" ,verbose );
	  LG->write("h                    : "+to_string(SC.w ) + "\n" ,verbose );
	  LG->write("threshold            : "+to_string(SC.threshold) + "\n\n" ,verbose );
      
      SC.dump();
	}
	else{
	  SC.mean = 0.6, SC.std = 0.001; //this dependent on -w 0.9 !!!
	  SC.set_2(pw->llrthresh);
	}

	//================================================
	//computing this SC via multiple times of the EM
	//run, should probably write something to combind
	//them all; i.e. MPI

	//=================================================
	//(2b) so segments is indexed by inidividual chromosomes, want to broadcast 
	//to sub-processes and have each MPI call run on a subset of segments
	vector<segment*> all_segments  	= segments;
	LG->write("slicing segments........................................", verbose);
	segments 						= MPI_comm::slice_segments(segments, rank, nprocs);	
	LG->write("done\n", verbose);

	//===========================================================================
	//(3a) now going to run the template matching algorithm based on pseudo-
	//moment estimator and compute BIC ratio (basically penalized LLR)
	LG->write("running template matching algorithm.....................", verbose);
	double threshold 	= run_global_template_matching_pwrapper(segments, out_file_dir, pw, SC);	
	//(3b) now need to send out, gather and write bidirectional intervals 
	LG->write("done\n", verbose);
	


	LG->write("scattering predictions to other MPI processes...........", verbose);
	int total =  MPI_comm::gather_all_bidir_predicitions_pwrapper(all_segments, 
							     segments , rank, nprocs, out_file_dir, job_name, job_ID, pw, 0);
	MPI_Barrier(MPI_COMM_WORLD); //make sure everybody is caught up!

	LG->write("done\n", verbose);
	if (rank==0){
	  LG->write("\nThere were " +to_string(total) + " prelimary bidirectional predictions\n\n", verbose);
	}
	
	//===========================================================================
	//this should conclude it all
	LG->write("clearing allocated segment memory.......................", verbose);	
	load::clear_segments(all_segments);
	LG->write("done\n", verbose);
	//===========================================================================
	//(4) if MLE option was provided than need to run the model_main::run()
	//
	if (pw->mle){
		pw->regionsOfInterest 	= pw->outputDir+job_name+ "-" + to_string(job_ID)+ "_prelim_bidir_hits.bed";
        printf("Regions of interest file..............................%s\n", pw->regionsOfInterest.c_str());
        
		model_run_pwrapper(pw, rank, nprocs,0, job_ID, LG);
		
	}
	LG->write("exiting bidir module....................................done\n\n", verbose);
	return 1;
}

/** Performs bidirectional predictions given a ParamWrapper utilizing an older algorithm.
 * During its development history, tfit's bidir module and backend were changed multiple times. This
 * function exists to emulate the functionality present in an older instance of the module, since some of the
 * changes in behavior since have been undesirable for various tasks.
 * @param pw ParamWrapper containing parsed command line arguments.
 * @param rank MPI rank parameter. This is used to determine the level of output a given instance of this function will generate.
 * @param nprocs Number of processors available. This should, in practice, be an MPI runtime parameter.
 * @param job_ID Current job identifier. This should, in practice, be an MPI runtime parameter.
 * @param LG Log_File object pointing to a valid output log file.
 * @return An output status code. 0 indicates no errors during operation. Any other value indicates errors.
 */
int bidir_run_old_long_pwrapper(ParamWrapper *pw, int rank, int nprocs, int job_ID, Log_File * LG){


	int verbose 	= pw->verbose;
	//P->p["-merge"] 	= "1";
	LG->write("\ninitializing bidir module...............................done\n", verbose);
	int threads 	= omp_get_max_threads();//number of OpenMP threads that are available for use
	
	//===========================================================================
	//get job_ID and open file handle for log files
	string job_name = pw->jobName;
	
	//===========================================================================
	//input files and output directories
	string forward_bedgraph 	= pw->forwardStrand; //forward strand bedgraph file
	string reverse_bedgraph 	= pw->reverseStrand; //reverse strand bedgraph file
	string joint_bedgraph 		= pw->mergedStrand; //joint forward and reverse strand bedgraph file
	string tss_file 		= pw->promoterTSS;
	string out_file_dir 		= pw->outputDir;//out file directory
	//===========================================================================
	//template searching parameters
	double sigma, lambda, foot_print, pi, w;
	double ns 	= pw->ns;
	sigma 		= pw->sigma, lambda=(double)pw->lambda;
    printf("%lf, %lf", lambda, pw->lambda);
	foot_print 	= pw->footPrint, pi=pw->pi, w=pw->w;



	//(2a) read in bedgraph files 
	map<string, int> chrom_to_ID;
	map<int, string> ID_to_chrom;
	vector<double> parameters;
	if (tss_file!="" and rank == 0){
        printf("tss_file: %s\n", tss_file.c_str());
		vector<segment *> FSI;
		LG->write("loading TSS intervals...................................",verbose);
		map<int, string> IDS;
		vector<segment *> tss_intervals 	= load::load_intervals_of_interest_pwrapper(tss_file, 
			IDS, pw, 1 );
		map<string, vector<segment *>> GG 	= MPI_comm::convert_segment_vector(tss_intervals);
		LG->write("done\n", verbose);
	
		vector<segment*> integrated_segments= load::insert_bedgraph_to_segment_joint(GG, 
			forward_bedgraph, reverse_bedgraph, joint_bedgraph, rank);
		LG->write("Binning/Normalizing TSS intervals.......................",verbose);
		load::BIN(integrated_segments, pw->br, pw->ns,true);	
		LG->write("done\n", verbose);
		
		LG->write("Computing Average Model.................................",verbose);
		parameters  	= compute_average_model_pwrapper(integrated_segments, pw);
		LG->write("done\n", verbose);
		sigma 	= parameters[1], lambda= parameters[2];
		foot_print= parameters[3], pi= parameters[4], w= parameters[5];
		LG->write("\nAverage Model Parameters\n", verbose);
		LG->write("-sigma      : " + to_string(sigma*ns)+ "\n", verbose);
        LG->write("-ns         : " + to_string(ns)+"\n", verbose);
        LG->write("-lambda (pw): " + to_string(lambda)+"\n", verbose);
		LG->write("-lambda     : " + to_string(ns/lambda)+ "\n", verbose);
		LG->write("-foot_print : " + to_string(foot_print*ns)+ "\n", verbose);
		LG->write("-pi         : " + to_string(pi)+ "\n", verbose);
		LG->write("-w          : " + to_string(w)+ "\n\n", verbose);
	}
	
	else
    {
        parameters.push_back(sigma);
        parameters.push_back(lambda);
        parameters.push_back(foot_print);
        parameters.push_back(pi);
        parameters.push_back(w);LG->write("-sigma      : " + to_string(sigma*ns)+ "\n", verbose);
        LG->write("-ns         : " + to_string(ns)+"\n", verbose);
        LG->write("-lambda (pw): " + to_string(lambda)+"\n", verbose);
		LG->write("-lambda     : " + to_string(ns/lambda)+ "\n", verbose);
		LG->write("-foot_print : " + to_string(foot_print*ns)+ "\n", verbose);
		LG->write("-pi         : " + to_string(pi)+ "\n", verbose);
		LG->write("-w          : " + to_string(w)+ "\n\n", verbose);
    }
    
	pw->sigma 			= sigma;
	pw->lambda=lambda;
	pw->footPrint=foot_print;
	pw->pi=pi;
	pw->w=w;
	MPI_comm::send_out_parameters( parameters, rank, nprocs);

	LG->write("loading bedgraph files..................................", verbose);
	vector<segment *> 	segments 	= load::load_bedgraphs_total(forward_bedgraph, reverse_bedgraph, joint_bedgraph,
		pw->br, pw->ns, pw->chromosome, chrom_to_ID, ID_to_chrom );

	if (segments.empty()){
		printf("exiting...\n");
		return 1;
	}
	LG->write("done\n", verbose);
	//(2b) so segments is indexed by inidividual chromosomes, want to broadcast 
	//to sub-processes and have each MPI call run on a subset of segments
	vector<segment*> all_segments  	= segments;
	LG->write("slicing segments........................................", verbose);
	segments 						= MPI_comm::slice_segments(segments, rank, nprocs);	
	LG->write("done\n", verbose);
	//===========================================================================
	//(3a) now going to run the template matching algorithm based on pseudo-
	//moment estimator and compute BIC ratio (basically penalized LLR)

	LG->write("running moment estimator algorithm......................", verbose);
	run_global_template_matching_old_long(segments, out_file_dir, pw);	
	//(3b) now need to send out, gather and write bidirectional intervals 
	LG->write("done\n", verbose);
	LG->write("scattering predictions to other MPI processes...........", verbose);
	int total =  MPI_comm::gather_all_bidir_predicitions_pwrapper(all_segments, 
			segments , rank, nprocs, out_file_dir, job_name, job_ID,pw,0);
	MPI_Barrier(MPI_COMM_WORLD); //make sure everybody is caught up!

	LG->write("done\n", verbose);
	if (rank==0){
		LG->write("\nThere were " +to_string(total) + " prelimary bidirectional predictions\n\n", verbose);
	}
	
	//===========================================================================
	//this should conclude it all
	LG->write("clearing allocated segment memory.......................", verbose);	
	load::clear_segments(all_segments);
	LG->write("done\n", verbose);
	//===========================================================================
	//(4) if MLE option was provided than need to run the model_main::run()
	//
	if (pw->mle){
		pw->regionsOfInterest 	= pw->outputDir+job_name+ "-" + to_string(job_ID)+ "_prelim_bidir_hits.bed";
        printf("Regions of interest file..............................%s\n", pw->regionsOfInterest.c_str());
        
		model_run_pwrapper(pw, rank, nprocs,0, job_ID, LG);
	}
	LG->write("exiting bidir module....................................done\n\n", verbose);
	return 1;
}

/** Performs bidirectional predictions given a ParamWrapper utilizing an older algorithm.
 * During its development history, tfit's bidir module and backend were changed multiple times. This
 * function exists to emulate the functionality present in an older instance of the module, since some of the
 * changes in behavior since have been undesirable for various tasks.
 * 
 * Please note that this is an even older and potentially incomplete representation of the bidir algorithm.
 * As such, it should not be used until extensive testing is conducted.
 * @param pw ParamWrapper containing parsed command line arguments.
 * @param rank MPI rank parameter. This is used to determine the level of output a given instance of this function will generate.
 * @param nprocs Number of processors available. This should, in practice, be an MPI runtime parameter.
 * @param job_ID Current job identifier. This should, in practice, be an MPI runtime parameter.
 * @param LG Log_File object pointing to a valid output log file.
 * @return An output status code. 0 indicates no errors during operation. Any other value indicates errors.
 */
int bidir_old_run_pwrapper(ParamWrapper * pw, int rank, int nprocs, int job_ID, Log_File * LG){
	int verbose 	= pw->verbose;
	//P->p["-merge"] 	= "1";
	LG->write("\ninitializing bidir module...............................done\n", verbose);
	int threads 	= omp_get_max_threads();//number of openMP threads that are available for use
	
	//===========================================================================
	//get job_ID and open file handle for log files
	string job_name = pw->jobName; //P->p["-N"];
	
	//===========================================================================
	//input files and output directories
	string forward_bedgraph 	= pw->forwardStrand;//P->p["-i"]; //forward strand bedgraph file
	string reverse_bedgraph 	= pw->reverseStrand; //reverse strand beddgraph file
	string combined_bedgraph    = pw->mergedStrand;
	string out_file_dir 		= pw->outputDir;//out file directory
	//TODO: copy bedgraph splitting code from Fstitch to accomodate this function.
	
	//(2a) read in bedgraph files 
	map<string, int> chrom_to_ID;
	map<int, string> ID_to_chrom;
	LG->write("loading bedgraph files..................................", verbose);
      
	vector<segment *> 	segments 	= load::load_bedgraphs_total(forward_bedgraph, reverse_bedgraph, combined_bedgraph,
		pw->br, pw->ns, pw->chromosome, chrom_to_ID, ID_to_chrom );
	LG->write("done\n", verbose);
	//(2b) so segments is indexed by inidividual chromosomes, want to broadcast 
	//to sub-processes and have each MPI call run on a subset of segments
	vector<segment*> all_segments  	= segments;
	LG->write("slicing segments........................................", verbose);
	segments 						= MPI_comm::slice_segments(segments, rank, nprocs);	
	LG->write("done\n", verbose);
	//===========================================================================
	//(3a) now going to run the template matching algorithm based on pseudo-
	//moment estimator and compute BIC ratio (basically penalized LLR)

	LG->write("running moment estimator algorithm......................", verbose);
    //NOTE as of 2018-05-02: This has been replaced with newer code in the run_global_template_matching_old_long
    //function. This should implement the functionality seen in time for the 2017-2018 paper based on 
    //commit history. 
	//run_global_template_matching_old(segments, out_file_dir, 4, 
	//		0.,pw->ns,pw->llrthresh, threads,0. ,0 );	
    run_global_template_matching_old_long(segments, out_file_dir, pw);
	//(3b) now need to send out, gather and write bidirectional intervals 
	LG->write("done\n", verbose);
	LG->write("scattering predictions to other MPI processes...........", verbose);
	int total =  MPI_comm::gather_all_bidir_predicitions_pwrapper(all_segments, 
			segments , rank, nprocs, out_file_dir, job_name, job_ID,pw,0);
	LG->write("done\n", verbose);
	if (rank==0){
		LG->write("\nThere were " +to_string(total) + " prelimary bidirectional predictions\n\n", verbose);
	}
	
	//===========================================================================
	//this should conclude it all
	LG->write("clearing allocated segment memory.......................", verbose);	
	load::clear_segments(all_segments);
	LG->write("done\n", verbose);
	//===========================================================================
	//(4) if MLE option was provided than need to run the model_main::run()
	//
    /*
	if (pw->mle){
		string k = pw->outputDir + job_name+ "-" + to_string(job_ID)+ "_prelim_bidir_hits.bed";
		model_run(P, rank, nprocs,0, job_ID, LG);
		if (stoi(P->p["-select"])){
			string bidir_ms_file 	= P->p["-o"]+  P->p["-N"] + "-" + to_string(job_ID)+  "_divergent_classifications.bed";
			string file_name 		= P->p["-o"]+  P->p["-N"] + "-" + to_string(job_ID)+  "_K_models_MLE.tsv";
			remove( bidir_ms_file.c_str() );
			//query out_dir +  P->p["-N"] + "-" + to_string(job_ID)+  "_K_models_MLE.tsv"
			P->p["-q"] 	= P->p["-o"] +  P->p["-N"] + "-" + to_string(job_ID)+  "_K_models_MLE.tsv";
			select_run(P, rank, nprocs, job_ID, LG);
			LG->write("loading results (MLE)...................................",verbose);
			vector<segment_fits *> fits 		= load::load_K_models_out(file_name);
			LG->write("done\n",verbose);		
			
			LG->write("writing out results (model selection)...................",verbose);
			load::write_out_bidirectionals_ms_pen(fits, P,job_ID, 0 );
			LG->write("done\n",verbose);
	

		}
	}*/
	LG->write("exiting bidir module....................................done\n\n", verbose);
	return 1;
}
