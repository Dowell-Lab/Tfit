#ifndef template_H
#define template_H
#include "load.h"
#include <math.h>
#include <vector>
#include <iostream>
using namespace std;
vector<double> peak_bidirs(segment * );
int sample_centers(vector<double>, double);
void run_global_template_matching(vector<segment*>, 
	string, double, double, double, double, int, double, int,ofstream& );
void optimize(map<string, interval_tree *>, 
	vector<segment*>,double, int, 
	string,string, int );
extern double INF;
extern double  nINF;

#endif
