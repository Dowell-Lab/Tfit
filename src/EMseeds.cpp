/**
 * @file EMseeds.cpp
 * @author Robin Dowell 
 * @brief Contains the EM algorithm for model inference
 * @version 0.1
 * @date 2022-06-02
 * 
 */
#include "EMseeds.h"

#include <iostream>
#include "split.h"
#include "helper.h" // Random
#include "Intervals.h"   // gInterval

Seeds::Seeds(): numgen() {
   // initalization 
}

std::string Seeds::write_out() {
  return "No seeding yet!";
}

double Seeds::grabSeed() {
  std::discrete_distribution<int> ddistro{weights.begin(), weights.end()};
  return mu_seeds[ddistro(numgen.mt)];
}

void Seeds::setupRandomSeeds(int numseeds, gInterval *region) {
  double relativeStop = (region->stop - region->start)-1;
	std::uniform_real_distribution<double> udist(0, relativeStop);
  for (int i = 0; i < numseeds; i++) {
    mu_seeds.push_back(udist(numgen.mt));
    weights.push_back(1.);
  }
}

void Seeds::weightRandomly() {
  for (int i = 0; i < mu_seeds.size(); i++) {
    weights[i] = numgen.fetchProbability();  // Weighted by random probability
  }
}

void Seeds::grabSeedsfromBed12 (std::vector<std::string> lineArray) {
  if (lineArray.size() < 6) {  // accept BED3 or BED4
    return;   // This is a Bed3, bed4 or bed6 file!
  }
  if (lineArray.size() >= 12) { // Or should this be 11?
    if (weights.size() > 0) weights.clear();
    if (mu_seeds.size() > 0) mu_seeds.clear();
/*
      // Check these are correct for this interval?
      double start = stod(lineArray[6]);
      double stop = stod(lineArray[7]);
    */

    int numseeds = stod(lineArray[9]);
    std::vector<std::string> seedsArray; // Contents of field[12], split on comma (,)

    // Read in weights from field 11:
    seedsArray = string_split(lineArray[10], ',');
    if (seedsArray.size() != numseeds) {
      std::cout << "ERROR: Invalid BED12, field 11 should contain # elements as specified in field 10!\n" << std::endl;
    }
    weights.reserve(numseeds);
    for (int i = 0; i < numseeds; i++) {
      weights[i] = stod(seedsArray[i]);
    }
    seedsArray.clear();

    // Read in seed positions from field 12:
    seedsArray = string_split(lineArray[11], ',');
    if (seedsArray.size() != numseeds) {
      std::cout << "ERROR: Invalid BED12, field 12 should contain # elements as specified in field 10!\n" << std::endl;
    }
    mu_seeds.reserve(numseeds);
    for (int i = 0; i < numseeds; i++) {
      mu_seeds[i] = stod(seedsArray[i]);
    }
  }
}

/**
 * @brief We can write out seeds as a bed12 line
 * 
 * As bed12:  bed6 is the interval: chr start stop identifier score strand
 * field 7  thick_start   typically start of CDS or repeats start
 * field 8  thick_end     typically end of CDS or repeats stop
 * field 9  RGB value     0,0,0 formatted
 * field 10 blockCount (#exons)
 * field 11 blockSizes
 * field 12 blockStarts   positions relative to start, e.g. first is zero
 * 
 * Using blockStarts as seed locations and blockSizes as relative weighting.
*/
std::string Seeds::writeToBed12 (bed12 *region) {
  std::string output = region->write_asBED();   // Need to force bed6!
  output += writeHalfBed12(region->start, region->stop);
  return output;
}  

std::string Seeds::writeHalfBed12(double start, double stop) {
  std::string output;
  // No seeds, no second half, i.e. stick with Bed6!
  if (mu_seeds.size() == 0) {   return output; }
  // Otherwise we're going to output this as 2nd half of a bed12:
  output += "\t" + tfit::prettyDecimal(start, 0);              // field 7
  output += "\t" + tfit::prettyDecimal(stop, 0);            // field 8
  output += "\t00,00,00";                                   // field 9, color = black #000000
  output += "\t" + tfit::prettyDecimal(mu_seeds.size(),0);  // field 10: # seeds
  output += "\t";
  if (mu_seeds.size() > 0) {  // field 11: all size = 1
    output += tfit::prettyDecimal(weights[0]*100, 0);
    for (int i = 1; i < mu_seeds.size(); i++) { 
      output += "," + tfit::prettyDecimal(weights[i]*100, 0); 
    }
  }
  output += "\t";
  if (mu_seeds.size() > 0) {  // field 12: all positions relative to start
    output += tfit::prettyDecimal(mu_seeds[0], 0);
    for (int i = 1; i < mu_seeds.size(); i++) { 
      output += "," + tfit::prettyDecimal(mu_seeds[i],0); 
    }
  }
  return output;
}
