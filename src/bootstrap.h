/**
 * @file bootstrap.h
 * @author Joey Azofeifa
 * @brief 
 * @version 0.1
 * @date 2016-05-20
 * 
 */
#ifndef bootstrap_H
#define bootstrap_H

#include <fstream>
#include <iostream>
#include <vector>

#include "load.h"
#include "read_in_parameters.h"

using namespace std;


void run_bootstrap_across(vector<segment *>, params *, ofstream&  );

#endif
