#ifndef model_main_H
#define model_main_H
#include "read_in_parameters.h"
#include "error_stdo_logging.h"
#include "ParamWrapper.hpp"


int model_run(params *, int, int, double, int, Log_File * );
int model_run_pwrapper(ParamWrapper *, int, int, double, int, Log_File *);

#endif 
