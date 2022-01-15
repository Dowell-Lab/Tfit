/**
 * @file FDR.h
 * @author Joey Azofeifa
 * @brief 
 * @version 0.1
 * @date 2016-11-22
 * 
 */
#ifndef FDR_H
#define FDR_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include "omp.h"

#include "load.h"
#include "model.h"
#include "read_in_parameters.h"

using namespace std;

class normal{
 public:
  double mean, std, x,threshold, c2;
  normal();
  normal(double, double, double);
  double cdf(double);
  double pdf(double); 
};

class exponential{
 public:
  double lambda,mu;
  exponential();
  exponential(double, double);
  double pdf(double );
};

class pareto{
 public:
  double alpha,mu;
  pareto();
  pareto(double, double);
  double pdf(double );
};




class slice_ratio{
 public:
  double start, stop ; //these should be base ten
  int bins ; //the number of segments
  double mean , std , w,c,threshold ;
  bool converged;
  vector<vector<double> > XY ; //bins X 3
  vector<normal> NORMS;
  normal norm_all;
  slice_ratio(double, double, int);
  slice_ratio();
  void set(double); //this computes the mean/std of each and makes the normal class
  void set_2(double);
  void insert(double);
  double pvalue(double);
  int get_closest(double);
};
slice_ratio get_slice(vector<segment *> , int,double,params * P );

#endif
