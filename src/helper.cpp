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

double Random::fetchExponential(double lambda) {
	std::exponential_distribution<double> edist(lambda);
	return edist(mt);
}

double Random::fetchProbability() {
	return Random::fetchUniform(0,1);
}

Bimap::Bimap() {
   num_elements = 0;
}

void Bimap::addIdentifier(std::string name) {
   if (!(str2index.count(name))) {  // a new chromosome
    // Lets add a new chromosome to the index list
    str2index[name] = num_elements; 
    // Reverse index
    index2str[num_elements] = name;
    num_elements++;
   }
}

int Bimap::lookupIndex(std::string name) {
   if (str2index.count(name)) {  // exists check
     return str2index[name];
   } else {
      return -1;
   }
}

std::string Bimap::lookupName(int index) {
   return index2str[index];
}

std::string Bimap::print_index_names() {
  std::string output;
  for (int i=0; i< num_elements; i++) {
    if (i > 0) { output += " ";}
    output += index2str[i];
  }
  return output;
}

std::string tfit::prettyDecimal (double x, double sigfig) {
   int sfig = sigfig + 1;      // Necessary to get decimal + signfig!
  return std::to_string(x).substr(0,std::to_string(x).find(".") + sfig);
}
