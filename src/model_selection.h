/**
 * @file model_selection.h
 * @author Joey Azofeifa
 * @brief 
 * @version 0.1
 * @date 2016-05-20
 * 
 */
#ifndef model_selection_H
#define model_selection_H
#include <vector>
#include "load.h"
using namespace std;

double ROC(vector<segment_fits *>, vector<segment_fits *>, double&, double&, double&, double&, double&, double& );

#endif
