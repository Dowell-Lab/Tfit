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


Seeds::Seeds(): mu_seeds() {
  // initialize seeds 
}

std::string Seeds::write_out() {
  return writeSeedsAsBedFields();
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

/**
 * @brief Given the last three fields of a bed12 file (as strings), 
 *   populate the seed data structure.
 * 
 * @param v_numSeeds      field 9 (num exons)
 * @param v_weightsSeeds  field 10 (exon lengths)
 * @param v_relativeSeedStarts  field 11 (exon starts, relative to start)
 */
void Seeds::getSeedsfromBedFields (std::string v_numSeeds, 
        std::string v_weightsSeeds, std::string v_relativeSeedStarts) {
  
  int numseeds = stod(v_numSeeds);
  std::vector<std::string> weightsArray; // Contents of field[11], split on comma (,)
  std::vector<std::string> seedsArray;   // Contents of field[11], split on comma (,)

  weightsArray = string_split(v_weightsSeeds, ',');
  seedsArray = string_split(v_relativeSeedStarts, ',');

  // Note numseeds = size() + 1
  if (weightsArray.size() != numseeds) {
    std::cout << "ERROR: Invalid BED12!" << std::endl;
    std::cout << "Field 11 should contain # elements as specified in field 10!\n"
                << std::endl;
  }
  if (seedsArray.size() != numseeds) {
    std::cout << "ERROR: Invalid BED12, field 12 should contain # elements as specified in field 10!\n"
                << std::endl;
  }

  mu_seeds.reserve(numseeds);
  for (int i = 0; i < numseeds; i++) {
    PointCov singleSeed(stod(seedsArray[i]), stod(weightsArray[i]) / 100);
    mu_seeds.push_back(singleSeed);
  }
}

/**
 * @brief Write out the seeds as the last three fields of a typical bed12.
 * Does NOT lead with a tab.
 * 
 * @return std::string 
 */
std::string Seeds::writeSeedsAsBedFields() {
  std::string output = "";
  output += tfit::prettyDecimal(getNumSeeds(),0);  // field 10: # exons (seeds)
  output += "\t";
  if (getNumSeeds() > 0) {  // field 11: exon lengths (used here for weights) 
    output += tfit::prettyDecimal(mu_seeds[0].coverage*100, 0);
    for (int i = 1; i < mu_seeds.size(); i++) { 
      output += "," + tfit::prettyDecimal(mu_seeds[i].coverage*100, 0); 
    }
  }
  output += "\t";
  if (getNumSeeds() > 0) {  // field 12: seed/exon positions relative to start
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