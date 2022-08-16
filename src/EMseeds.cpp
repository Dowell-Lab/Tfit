/**
 * @file EMseeds.cpp
 * @author Robin Dowell 
 * @brief Contains the EM algorithm for model inference
 * @version 0.1
 * @date 2022-06-02
 * 
 */
#include "EMseeds.h"

#include <algorithm>
#include <string>
#include <vector>
#include <iostream>

#include "split.h"
#include "helper.h" // Random
#include "Intervals.h"   // gInterval


Seeds::Seeds() {
  // initialize seeds 
}

std::string Seeds::write_out() {
  return "No seeding yet!";
}

double Seeds::getMaxWeight() {
  double max = -1;
  for (int i = 0; i < getNumSeeds(); i++) {
    if (mu_seeds[i].coverage > max) {
       max = mu_seeds[i].coverage;
    }
  }
  return max; 
}

double Seeds::getMinWeight() {
  double min = INT32_MAX;
  for (int i = 0; i < getNumSeeds(); i++) {
    if (mu_seeds[i].coverage < min) {
       min = mu_seeds[i].coverage;
    }
  }
  return min; 
}

int Seeds::getNumSeeds() {
  return (int)mu_seeds.size();
}

void Seeds::grabSeedsfromBed12 (std::vector<std::string> lineArray) {
  if (lineArray.size() < 12) {  // accept BED3 or BED4
    return;   // This is a Bed3, bed4 or bed6 file!
  }
  if (lineArray.size() >= 12) { // Or should this be 11?
    if (mu_seeds.size() > 0) mu_seeds.clear();  // Removes old seeds if exist.
/*
      // Check these are correct for this interval?
      double start = stod(lineArray[6]);
      double stop = stod(lineArray[7]);
    */

    int numseeds = stod(lineArray[9]);
    std::vector<std::string> weightsArray; // Contents of field[11], split on comma (,)
    std::vector<std::string> seedsArray; // Contents of field[11], split on comma (,)

    weightsArray = string_split(lineArray[10], ',');
    seedsArray = string_split(lineArray[11], ',');
    if (weightsArray.size() != numseeds) {
      std::cout << "ERROR: Invalid BED12, field 11 should contain # elements as specified in field 10!\n" << std::endl;
    }
    if (seedsArray.size() != numseeds) {
      std::cout << "ERROR: Invalid BED12, field 12 should contain # elements as specified in field 10!\n" << std::endl;
    }

    mu_seeds.reserve(numseeds);
    for (int i = 0; i < numseeds; i++) {
      PointCov singleSeed(stod(seedsArray[i]), stod(weightsArray[i]));
      mu_seeds.push_back(singleSeed);
    }
  }
}

std::string Seeds::writeHalfBed12(double start, double stop) {
  std::string output = "";
  // No seeds, no second half, i.e. stick with Bed6!
  if (mu_seeds.size() == 0) {   return output; }
  // Otherwise we're going to output this as 2nd half of a bed12:
  output += "\t" + tfit::prettyDecimal(start, 0);              // field 7
  output += "\t" + tfit::prettyDecimal(stop, 0);            // field 8
  output += "\t00,00,00";                                   // field 9, color = black #000000
  output += "\t" + tfit::prettyDecimal(mu_seeds.size(),0);  // field 10: # seeds
  output += "\t";
  if (getNumSeeds() > 0) {  // field 11: weights
    output += tfit::prettyDecimal(mu_seeds[0].coverage*100, 0);
    for (int i = 1; i < mu_seeds.size(); i++) { 
      output += "," + tfit::prettyDecimal(mu_seeds[i].coverage*100, 0); 
    }
  }
  output += "\t";
  if (getNumSeeds() > 0) {  // field 12: all positions relative to start
    output += tfit::prettyDecimal(mu_seeds[0].coordinate, 0);
    for (int i = 1; i < mu_seeds.size(); i++) { 
      output += "," + tfit::prettyDecimal(mu_seeds[i].coordinate,0); 
    }
  }
  return output;
}

void Seeds::SortByWeights() {
  std::sort(mu_seeds.begin(),mu_seeds.end(), PointCov::sortOnCovComp);
}

void Seeds::SortByPositions() {
  std::sort(mu_seeds.begin(),mu_seeds.end(), PointCov::sortOnCoordComp);
}

/***** SeedManager *****/

SeedManager::SeedManager(): numgen() {
  // initialize SeedManager
}

std::string SeedManager::write_out() {
  return "nothing yet";
}

double SeedManager::grabSeed() {
  std::discrete_distribution<int> ddistro{setSeeds->getMinWeight(), setSeeds->getMaxWeight()};
  return setSeeds->mu_seeds[ddistro(numgen.mt)].coordinate;
}

void SeedManager::setupRandomSeeds(int numseeds, gInterval *region) {
  double relativeStop = (region->stop - region->start)-1; // Is this correct?
	std::uniform_real_distribution<double> udist(0, relativeStop);
  
  if (setSeeds != NULL) {   // Seeds already exist, what do I do?  ** current refuse!
     return;
  } else {
    Seeds newseeds;
    setSeeds = &newseeds;    // Create Seeds object;
    for (int i = 0; i < numseeds; i++) {
      setSeeds->mu_seeds.push_back(PointCov(udist(numgen.mt), 1.));
    }
  }
}

void SeedManager::weightRandomly() {
  for (int i = 0; i < setSeeds->getNumSeeds(); i++) {
    // Weighted by random probability, note doesn't check if already weighted!!
    setSeeds->mu_seeds[i].coverage = numgen.fetchProbability();  
  }
}