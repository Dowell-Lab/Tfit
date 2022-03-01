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

extern int g_testing;

Random::Random() { 
   int truerandom = 1;
   if (g_testing) { truerandom = 0; }		// ugh a global!

   std::random_device rd;

   if (truerandom) {
    mt.seed(rd());
   } else {
    mt.seed(1974);
   }
}

double Random::fetchUniform(double lower, double upper) {
	std::uniform_real_distribution<double> udist(lower, upper);
	return udist(mt);
}

double Random::fetchNormal(double mean, double std) {
	std::normal_distribution<double> ndist(mean, std);
	return ndist(mt);
}

double Random::fetchProbability() {
	return Random::fetchUniform(0,1);
}

