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

/**
 * @brief Program main.  
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

  nprocs = 1;  rank = 0; threads = 1;  // Force MPI to minimum

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

  Log_File * LG = new  Log_File(rank, job_ID, P->p["-N"], P->p["-log_out"]);
  if (verbose > 0) {
  // printf("This should be output when verbose is set!\n");
  }
  if (verbose && rank == 0) {
    P->display(nprocs, threads);
    LG->write(P->get_header(1), 0);
  }
  if (P->bidir) {
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
