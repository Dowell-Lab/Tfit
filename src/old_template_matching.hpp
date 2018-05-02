#ifndef OLD_TEMPLATE_MATCHING_HPP
#define OLD_TEMPLATE_MATCHING_HPP

#include "load.h"
#include "model.h"
#ifdef USING_ICC
#include <mathimf.h>
#else
#include <math.h>   
#endif
#include <limits>
#include <iostream>
#include <algorithm>
#include "template_matching.h"
#include <fstream>
#include <random>
#include "omp.h"
#ifdef USING_ICC
#include <aligned_new>
#endif

#include "ParamWrapper.hpp"

//vector<vector<double>> bubble_sort3(vector<vector<double>> X);
double BIC2(double ** X,  double * avgLL, double * variances,double * lambdas, 
	double ** skews, double N_pos, double N_neg, double S_pos, 
		double S_neg, double S2_pos, double S2_neg,double mu, int j,int k,int i, double scale );
double get_ll(double ** X, double mu, double w, double pi, double l, int j, int k);
double BIC3(double ** X, int j, int k, int i,
	double N_pos, double N_neg, double * avgLL, double * variances,double * lambdas, 
	double ** skews);
void BIC_template(segment * data, double * avgLL, double * BIC_values, double * densities, double * densities_r,
	double * variances,double * lambdas, double ** skews ,double window, int np, int single,double foot_res,double scale);
void run_global_template_matching_old(vector<segment*> segments, 
	string out_dir,  double res, double density,
	double scale, double ct, int np, double skew, int single);
void noise_global_template_matching_old(vector<segment*> segments, double scale);
void run_global_template_matching_old_long(vector<segment*> segments, string out_dir, ParamWrapper *pw);

#endif
