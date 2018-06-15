#ifndef bootstrap_H
#define bootstrap_H
#include "load.h"
#include <vector>
#include "read_in_parameters.h"
#include <iostream>
#include <fstream>

using namespace std;

/** Performs initial data acquisition and computations given a set of segments. This includes computing internal variances as well as setting all relevant segment parameters to match those passed in the (depreciated) params object passed. This function is unused in the present Tfit codebase.
 * @param segments Set of segments to initialize.
 * @param P Command line parameters. 
 * @param log_file Output log file represented as an ofstream.
 */
void run_bootstrap_across(vector<segment *>, params *, ofstream&  );

#endif
