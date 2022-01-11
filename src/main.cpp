/**
 * @file main.cpp
 * @author Joey Azofeifa
 * @brief This is the primary executable file.  
 * It contains the mpi code, reading parameters, then forking to one of three functions:
 * \ref bidir_run, \ref model_run or \ref select_run -- which are all separate files. 
 * @version 0.1
 * @date 2016-05-20
 * 
 */
#include <mpi.h>
#include "load.h"
#include "model.h"
#include <iostream>
#include "across_segments.h"
#include <limits>
#include <math.h>   
#include <errno.h>
#include "error_stdo_logging.h"
#include <time.h>  
#include <stdio.h>   
#include <chrono>
#include <map>
#include "read_in_parameters.h"
#include "model_selection.h"
#include <thread>
#include "template_matching.h"
#include "MPI_comm.h"
#include "bootstrap.h"
#include "density_profiler.h"
#include <omp.h>
#include "bidir_main.h"
#include "model_main.h"
#include "select_main.h"
using namespace std;

/**
 * @brief Program main.  
 * It contains the mpi code, reading parameters, then forking to one of three modules. 
 * @param argc
 * @param argv
 * @return 
 */
int main(int argc, char* argv[]){
  MPI::Init(argc, argv);
  // the nprocs is the total number of processors available
  int nprocs		= MPI::COMM_WORLD.Get_size();
  // The rank ranges from 0 to (nprocs-1) and is current process number
  int rank 		= MPI::COMM_WORLD.Get_rank();
  int threads  	= omp_get_max_threads();

  params * P 	= new params();
  read_in_parameters(argv, P, rank);
  if (P->EXIT){
    if (rank == 0){
      printf("exiting...\n");
    }
    delete P;
    MPI::Finalize();
    return 0;
  }
  int job_ID 		=  MPI_comm::get_job_ID(P->p["-log_out"], P->p["-N"], rank, nprocs);
  
  int verbose 	= stoi(P->p["-v"]);
  Log_File * LG 	= new  Log_File(rank, job_ID, P->p["-N"], P->p["-log_out"]);
  if (verbose > 0) { 
     // printf("This should be output when verbose is set!\n"); 
   } 
 if (verbose and rank==0){ 
    P->display(nprocs,threads);
    LG->write(P->get_header(1),0);
  }
  if (P->bidir){
    bidir_run(P, rank, nprocs, job_ID,LG);
  } else if (P->model){
    model_run(P, rank, nprocs, 0, job_ID,LG);
  }else if (P->select){
    // This one is commented out -- ie. this function does nothing.
    select_run(P, rank, nprocs, job_ID,LG);	
  }
  if (rank == 0){
    load::collect_all_tmp_files(P->p["-log_out"], P->p["-N"], nprocs, job_ID);
  }
  
  MPI::Finalize();
  
  return 0;
}
