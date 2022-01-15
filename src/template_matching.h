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

vector<double> peak_bidirs(segment * );
int sample_centers(vector<double>, double);
void noise_global_template_matching(vector<segment*>, double);

double run_global_template_matching(vector<segment*> , string,  params * ,slice_ratio );
void EX(vector<segment*> , double, double , double & , double &);

extern double INF;
extern double  nINF;

#endif
