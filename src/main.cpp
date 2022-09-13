/**
 * @file main.cpp
 * @author Joey Azofeifa, Robin Dowell
 * @copyright Copyright 2016 Dowell Lab 
 * @brief This is the primary executable file.  
 * It contains the mpi code, reading parameters, then forking to one of three functions:
 * \ref bidir_run, \ref model_run or \ref select_run -- which are all separate files. 
 * @version 0.1
 * @date 2016-05-20
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

#include "across_segments.h"
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

int g_testing {};  // Ugh. Global variable to say we're NOT testing code here.

/**
 * @brief Program main.  
 * It contains the mpi code, reading parameters, then forking to one of three modules. 
 * @param argc
 * @param argv
 * @return 
 */
int main(int argc, char* argv[]) {
  MPI::Init(argc, argv);
  // the nprocs is the total number of processors available
  int nprocs  = MPI::COMM_WORLD.Get_size();
  // The rank ranges from 0 to (nprocs-1) and is current process number
  int rank    = MPI::COMM_WORLD.Get_rank();
  int threads = omp_get_max_threads();

  params * P  = new params();
  read_in_parameters(argv, P, rank);
  if (P->EXIT) {
    if (rank == 0) {
      printf("exiting...\n");
    }
    delete P;
    MPI::Finalize();
    return 0;
  }
  int job_ID = MPI_comm::get_job_ID(P->p["-log_out"], P->p["-N"], rank, nprocs);

  int verbose   = stoi(P->p["-v"]);

  int reqthreads = stoi(P->p["-threads"]);
  printf("THREADS: %d requested; %d available\n", reqthreads, threads);
  if ((reqthreads != 0) && (reqthreads <= threads)) {
    threads = reqthreads;
    P->threads = reqthreads;
  } else {
    P->threads = threads;
  }

  Log_File * LG = new  Log_File(rank, job_ID, P->p["-N"], P->p["-log_out"]);
  if (verbose > 0) {
  // printf("This should be output when verbose is set!\n");
  }
  if (verbose && rank == 0) {
    P->display(nprocs, threads);
    LG->write(P->get_header(1), 0);
  }
  if (P->testing) {
    // For testing behavior of subparts of Tfit.
    // std::cout << "Hijacked!" << std::endl;
    // bidir_rdd(P, rank, nprocs, job_ID, LG);
    model_rdd(P, rank, nprocs, 0, job_ID, LG);
  } else if (P->bidir) {
    bidir_run(P, rank, nprocs, job_ID, LG);
  } else if (P->model) {
    model_run(P, rank, nprocs, 0, job_ID, LG);
  } else if (P->select) {
    // Contents of select_run are commented out:
    // i.e. this function does nothing.
    select_run(P, rank, nprocs, job_ID, LG);
  }
  if (rank == 0) {
    load::collect_all_tmp_files(P->p["-log_out"], P->p["-N"], nprocs, job_ID);
  }

  MPI::Finalize();

  return 0;
}
