#ifndef model_main_H
#define model_main_H
#include "read_in_parameters.h"
#include "error_stdo_logging.h"
#include "ParamWrapper.hpp"

/** Utilizes an old style params struct and data computed from the bidir module to make predictions with respect to regions of transcription in a given input bedgraph file. This function is depreciated in favor of model_run_pwrapper, as the old style params struct interface is depreciated.
 * @depreciated
 * @param P params struct containing parsed command line arguments.
 * @param rank MPI rank parameter. This is used to determine the level of output a given instance of this function will generate.
 * @param nprocs Number of processors available. This should, in practice, be an MPI runtime parameter.
 * @param density Obstensibly used to set various output modes; in practice, this parameter is unused and does nothing.
 * @param job_ID Current job identifier. This should, in practice, be an MPI runtime parameter.
 * @param LG Log_File object pointing to a valid output log file.
 * @return An output status code. 0 indicates no errors during operation. Any other value indicates errors.
 */
int model_run(params *, int, int, double, int, Log_File * );

/** Utilizes a ParamWrapper and data computed from the bidir module to make predictions with respect to regions of transcription in a given input bedgraph file. This function is depreciated in favor of model_run_pwrapper, as the old style params struct interface is depreciated.
 * @depreciated
 * @param P params struct containing parsed command line arguments.
 * @param rank MPI rank parameter. This is used to determine the level of output a given instance of this function will generate.
 * @param nprocs Number of processors available. This should, in practice, be an MPI runtime parameter.
 * @param density Obstensibly used to set various output modes; in practice, this parameter is unused and does nothing.
 * @param job_ID Current job identifier. This should, in practice, be an MPI runtime parameter.
 * @param LG Log_File object pointing to a valid output log file.
 * @return An output status code. 0 indicates no errors during operation. Any other value indicates errors.
 */
int model_run_pwrapper(ParamWrapper *, int, int, double, int, Log_File *);

#endif 
