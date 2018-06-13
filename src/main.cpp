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

/** Main function for tfit-revisions.
 */
int main(int argc, char* argv[]){
  int rank, nprocs;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  
  
  int threads  	= omp_get_max_threads();

  //Todo: consider adding a rank parameter, since it appears to change how error messages are generated.
  ParamWrapper pw(argc, argv);
  //rank=0;
  
  if(pw.exit)
  {
      if(rank==0)
      {
          printf("Exiting...\n");
      }
      MPI_Finalize();
      return 0;
  }
  
  /*
  params * P 	= new params();
  read_in_parameters(argv, P, rank );
  if (P->EXIT){
    if (rank == 0){
      printf("exiting...\n");
    }
    delete P;
    MPI_Finalize();
    return 0;
  }*/
  int job_ID 		=  MPI_comm::get_job_ID(pw.logDir, pw.jobName, rank, nprocs); //P->p["-log_out"], P->p["-N"], rank, nprocs);
  
  int verbose 	= pw.verbose; //stoi(//P->p["-v"]);
  Log_File * LG 	= new Log_File(rank, job_ID, pw.jobName, pw.logDir); //new  Log_File(rank, job_ID, P->p["-N"], P->p["-log_out"]);
  if (verbose and rank==0){
      pw.display(nprocs, threads);
    //P->display(nprocs,threads);
    LG->write(pw.getHeader(1), 0);
  }
  //The code seems to be breaking riiiight here.
  if (pw.bidir){
      printf("About to launch bidir module.\n");
    //bidir_run(P, rank, nprocs, job_ID,LG);
    bidir_run_pwrapper(&pw, 0, nprocs, job_ID, LG);
  }
  else if(pw.bidirOld)
  {
      printf("About to launch bidir module with older behavior.\n");
      //bidir_old_run_pwrapper(&pw, rank, nprocs, job_ID, LG);
      bidir_run_old_long_pwrapper(&pw, rank, nprocs, job_ID, LG);
  }
  else if (pw.model){
      printf("About to launch model module.\n");
    //model_run(P, rank, nprocs,0,job_ID,LG);
    model_run_pwrapper(&pw, rank, nprocs, 0, job_ID, LG);
  }else if (pw.select){
    //select_run(P, rank, nprocs, job_ID,LG);	
      printf("About to launch select module.\n");
    select_run_pwrapper(&pw, rank, nprocs, job_ID, LG);
  }
  
  else
  {
      printf("Something exceptionally bad just happened!\n");
  }
  if (rank == 0){
    load::collect_all_tmp_files(pw.logDir, pw.jobName, nprocs, job_ID);
  }
  
  MPI_Finalize();
  
  return 0;
}
