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

double run_global_template_matching(vector<segment*> segments, 
            params * P, slice_ratio SC);
double RF_run_global_template_matching(vector<segment*> segments, 
                                    params * P, slice_ratio SC);

void RF_BIC_template(segment * data,  double * BIC_values, double * densities, double * densities_r, double window, 
		  double sigma, double lambda, double foot_print, double pi, double w, int thr);

double RF_BIC3(double ** X, int j, int k, int i,
	    double N_pos, double N_neg,  double sigma, double lambda, double fp, double pi, double w);


// Unused cruft?: 
// void EX(vector<segment*> , double, double , double & , double &);
// void noise_global_template_matching(vector<segment*>, double);
// vector<double> peak_bidirs(segment * );

extern double INF;
extern double  nINF;

#endif
