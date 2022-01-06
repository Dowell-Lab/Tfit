/**
 * @file error_stdo_logging.cpp
 * @author Joey Azofeifa
 * @brief 
 * @version 0.1
 * @date 2016-05-20
 * 
 */

/* Constructors: Log_File 
 * Author: Joey Azofeifa 
 *
 * Purpose: Helper functions for logging output to a designate file.
 *
 */
// Empty constructor
Log_File::Log_File(){}

// rank, job_ID, job_name, output_directory 
Log_File::Log_File(int R, int JI, string JN, string OUT){
	rank 	= R, job_ID = JI, job_name 	= JN, log_out_dir = OUT;
   	string log_out 	= OUT + "tmp_" + JN+ "-"+ to_string(job_ID)+ "_" + to_string(rank) + ".log"  ;
	FHW.open(log_out);
	FHW<<"Application out and error file for processes: " + to_string(rank)+ "\n";
}

void Log_File::write(string LINE, int verbose){
	if (verbose and rank==0){
		printf("%s", LINE.c_str() );
		cout.flush();
	}
	FHW<<LINE;
	FHW.flush();
}
