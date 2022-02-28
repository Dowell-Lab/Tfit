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

double Random::FetchUniform(double lower, double upper) {
	if (!testing) {
		std::uniform_real_distribution<double> udist(lower, upper);
		return udist(mt);
	} else {
      return rand(); 
	}
}

double Random::FetchNormal(double lower, double upper) {
	if (!testing) {
		std::normal_distribution<double> ndist(lower, upper);
		return ndist(mt);
	} else {
      return rand(); 
	}
}

double Random::FetchProbability() {
	return Random::FetchUniform(0,1);
}

