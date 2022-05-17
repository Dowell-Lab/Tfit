/**
 * @file template_matching.h
 * @author Joey Azofeifa
 * @brief 
 * @version 0.1
 * @date 2016-05-20
 * 
 */
#ifndef template_H
#define template_H

#include <iostream>
#include <math.h>
#include <vector>

#include "FDR.h"
#include "load.h"
#include "read_in_parameters.h"

using namespace std;

int sample_centers(vector<double>, double);

double run_global_template_matching(vector<segment*> , string,  params * ,slice_ratio );
double RF_run_global_template_matching(vector<segment*> segments, string out_dir,  
                                    params * P, slice_ratio SC);


// Unused cruft?: 
// void EX(vector<segment*> , double, double , double & , double &);
// void noise_global_template_matching(vector<segment*>, double);
// vector<double> peak_bidirs(segment * );

extern double INF;
extern double  nINF;

#endif
