#ifndef select_main_H
#define select_main_H
#include "read_in_parameters.h"
#include "error_stdo_logging.h"
#include "ParamWrapper.hpp"

int select_run(params *, int, int,int, Log_File *);
int select_run_pwrapper(ParamWrapper *, int, int, int, Log_File *);

#endif
