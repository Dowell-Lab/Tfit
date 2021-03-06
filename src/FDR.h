#ifndef FDR_H
#define FDR_H
#include <limits>
#include "load.h"
#include "model.h"
#include <iostream>
#include <algorithm>
#include <random>
#include "omp.h"
#include <cmath>
#include <vector>
#include <string>
#include "read_in_parameters.h"
#include "ParamWrapper.hpp"
#include <cstring>
using namespace std;

/** Represents a normal distribution.
 * This class provides methods that enable a user to sample from the distribution represented.
 */
class normal{
 public:
  double mean, std, x,threshold, c2;
  normal();
  normal(double, double, double);
  double cdf(double);
  double pdf(double); 
  void dump(char *header);
};

/** Represents an exponential distribution.
 * This class provides methods that enable a user to sample from an exponential distribution.
 */
class exponential{
 public:
  double lambda,mu;
  exponential();
  exponential(double, double);
  double pdf(double );
};

/** Represents a pareto distribution.
 * This class provides methods that enable a user to sample from a pareto distribution.
 */
class pareto{
 public:
  double alpha,mu;
  pareto();
  pareto(double, double);
  double pdf(double );
};

/** Represents an individual model component over a set of values.
 */
class slice_ratio{
 public:
  double start, stop ; //these should be base ten
  int bins ; //the number of segments
  double mean , std , w,c,threshold ;
  bool converged;
  double score;
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
  //TODO: Add log file support for dump():
  void dump();
};
slice_ratio get_slice(vector<segment *> , int,double,params * P );
slice_ratio get_slice_pwrapper(vector<segment *>, int, double, ParamWrapper *);

#endif
