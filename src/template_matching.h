#ifndef template_H
#define template_H
#include "load.h"
#include "read_in_parameters.h"
#include <math.h>
#include <vector>
#include <iostream>
#include "FDR.h"
#include "ParamWrapper.hpp"
using namespace std;

vector<double> peak_bidirs(segment * );

/** Samples a random integer from a unifrom distribution.
 * @param centers Set of distribution centers used to calculate the range of the returned value.
 * @param p Unused.
 * @return Random integer drawn from a uniform distribution.
 */
int sample_centers(vector<double> centers, double p);
void noise_global_template_matching(vector<segment*>, double);

/** Runs template matching given a depreciated params object.
 * This should effectively implement the functionality seen within the bidir module.
 * @param segments Set of input segments read from an input bedgraph.
 * @param out_dir Output directory in which to save the trained model.
 * @param P Depreciated params object.
 * @param SC slice_ratio used to keep track of various thresholds.
 */
double run_global_template_matching(vector<segment*> segments, string out_dir,  params *P, slice_ratio SC);

/** Runs template matching given a ParamWrapper object.
 * This should effectively implement the functionality seen within the bidir module.
 * @param segments Set of input segments read from an input bedgraph.
 * @param out_dir Output directory in which to save the trained model.
 * @param pw ParamWrapper from which the function reads a number of parameters.
 * @param SC slice_ratio used to keep track of various thresholds.
 */
double run_global_template_matching_pwrapper(vector<segment*> segments, string out_dir, ParamWrapper *pw, slice_ratio SC, int);
void EX(vector<segment*> , double, double , double & , double &);

/** Sorts all values within the first vector of vectors by the values in the second vector of vectors.
 * @param X Vector of vectors to sort.
 * @return Vector of sorted vectors.
 */
vector<vector<double>> bubble_sort3(vector<vector<double>> X);

extern double INF;
extern double  nINF;

#endif
