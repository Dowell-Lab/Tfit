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
  weightsArray = string_split(v_weightsSeeds, ',');

  std::vector<std::string> seedsArray;   // Contents of field[11], split on comma (,)
  seedsArray = string_split(v_relativeSeedStarts, ',');

  // Note numseeds = size() + 1
  if (weightsArray.size() != numseeds) {
    std::cout << "ERROR: Invalid BED12!" << std::endl;
    std::cout << "Field 11 should contain # elements as specified in field 10!\n"
                << std::endl;
    std::cout << to_string(weightsArray.size()) + " : " + to_string(numseeds) << std::endl;
  }
  if (seedsArray.size() != numseeds) {
    std::cout << "ERROR: Invalid BED12, field 12 should contain # elements as specified in field 10!\n"
                << std::endl;
  }

  //mu_seeds.reserve(numseeds+1);
  // mu_seeds.empty();   // Remove anything there already?
  double position, weight;
  for (int i = 0; i < numseeds; i++) {
    position = stod(seedsArray[i]);
    weight = stod(weightsArray[i]) / 100;
    PointCov singleseed(position, weight);
    mu_seeds.push_back(singleseed);
  }
  // cout << "After: " + to_string(mu_seeds.size()) << std::endl;
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

std::vector<double> SeedManager::grabSeedSet(int K) {
  std::vector<double> pickedseeds;
  if (setSeeds == NULL) { // No seeds allocated
     // Should allocate seeds, fill them (but this will need coords!)
     // Or throw an error?  
  }
  for (int k = 0; k < K; k++) {
    /*
   if seeds_still_exist, grab seed (delete it) randomly
     if randomize (add_noise) sample with Normal at seed, r_mu (stddev)
   else grab seed from Normal at midpoint, r_mu (stddev, defaults = 0)

    if (mu_seeds.size()>0  ){
      i 	= sample_centers(mu_seeds ,  p);
      mu 	= mu_seeds[i];
      if (r_mu > 0){   // if randomize (default is FALSE)
        mu = ran_num_generator.fetchNormal(mu,r_mu);
      }
    }else{
      // This appears to seed all to the midpoint
      mu = ran_num_generator.fetchNormal((data->minX+data->maxX)/2., r_mu);
    }
    mus[k] 	= mu;
    if (mu_seeds.size() > 0  ){
      mu_seeds.erase (mu_seeds.begin()+i);
    }

   */
  }
  return pickedseeds;
}

/*
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

*/