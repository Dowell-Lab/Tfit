#ifndef bidir_main_H
#define bidir_main_H
#include "read_in_parameters.h"
#include "error_stdo_logging.h"
#include "ParamWrapper.hpp"

int bidir_run(params *, int, int,int, Log_File *);
int bidir_run_pwrapper(ParamWrapper *, int, int, int, Log_File *);


#endif
