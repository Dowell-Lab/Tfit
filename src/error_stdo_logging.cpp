#include "error_stdo_logging.h"

/** Default constructor for the Log_File class. Does nothing.
 */
Log_File::Log_File(){}
/** Full constructor for the Log_File class.
 * @param R Rank parameter obtained through the MPI runtime.
 * @param JI Job ID obtained through the MPI runtime.
 * @param JN Job name passed to tfit on the command line.
 * @param OUT Log directory passed to tfit on the command line.
 */
Log_File::Log_File(int R, int JI, string JN, string OUT){
	rank 	= R, job_ID = JI, job_name 	= JN, log_out_dir = OUT;
   	string log_out 	= OUT + "tmp_" + JN+ "-"+ to_string(job_ID)+ "_" + to_string(rank) + ".log"  ;
	FHW.open(log_out);
	FHW<<"Application out and error file for processes: " + to_string(rank)+ "\n";
}
/** Writes the specified message to the log file.
 * If verbose is set, then the output will also be written to stdout.
 * This check prevents unnecessary messages from being generated while Tfit runs.
 * 
 * @param LINE String to write to the log
 * @param verbose Whether or not to write the line to standard output in addition to the log file.
 */
void Log_File::write(string LINE, int verbose){
	if (verbose and rank==0){
		printf("%s", LINE.c_str() );
		cout.flush();
	}
	FHW<<LINE;
	FHW.flush();
}
