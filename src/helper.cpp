/**
 * @file helper.cpp
 * @author Robin Dowell 
 * @brief 
 * @version 0.1
 * @date 2022-02-28
 * 
 */
#include "helper.h"

#include <iostream>
#include <sstream>
#include <random>	// random numbers within distributions

Random::Random() { 
   testing = 0;
   std::random_device rd;

   mt.seed(rd());
   srand(100);
}

Random::Random(int t_seed) {
   testing = 1;
   std::random_device rd;

   mt.seed(rd());
   srand(t_seed);
}

double Random::fetchUniform(double lower, double upper) {
	if (!testing) {
		std::uniform_real_distribution<double> udist(lower, upper);
		return udist(mt);
	} else {
      double rn =  lower + static_cast <double> (rand()) / 
	  	(static_cast <double> (RAND_MAX/(upper-lower)));
      return rn;
	}
}

double Random::fetchNormal(double mean, double std) {
	if (!testing) {
		std::normal_distribution<double> ndist(mean, std);
		return ndist(mt);
	} else {
	  // restrict to pseudo-uniform within 2 sd of mean
      double lower = mean - 2*std;
	  double upper = mean + 2*std;
	  return Random::fetchUniform(lower,upper);
	}
}

double Random::fetchProbability() {
	return Random::fetchUniform(0,1);
}

