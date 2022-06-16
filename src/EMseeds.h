/**
 * @file EMseeds.h
 * @author Robin Dowell 
 * @version 0.1
 * @date 2022-06-02
 * 
 */
#ifndef EMseeds_H 
#define EMseeds_H 

#include <string>
#include "helper.h"   //Random
#include "Intervals.h"  // gInterval

/**
 * @brief How do we want to generate and manage seeds?
 * 
 * Seeds are associated with a gInterval, use coordinates relative
 * to the gInterval (e.g. gInterval->start is position zero).
 * Seeds are constrained to reside within the gInterval.
 * 
 * Can read/write to files in Bed12 format using the "exon starts" as
 * seed locations (size = 1)
 * 
 */
class Seeds { 
  public:
  // set of seeds (read from file or from alg)
  std::vector<double> mu_seeds;    // Set of preferred seed points for mu
  std::vector<double> weights;    // Weights on the seeds, order all three by weights!

  // These could be private:
  Random numgen;    // Random number generator
  // std::vector<bool> used;    // Indicator, have we used this seed?

  // Constructor
  Seeds();

  //Functions
  std::string write_out();

  double grabSeed(); // Grab a seed (mark it as used?)

  // Build x random choosen seeds, equally weighted across region of interest
  void setupRandomSeeds(int numseeds, gInterval *region);  
  void weightRandomly();    // randomly generate weights for seeds

  // Run a seeding algorithm (template matching or other)

  // User Input (i.e. I/O of seeds, perhaps bed12?)
  void grabSeedsfromBed12 (std::vector<std::string> lineArray);
  bed6 setfromBedLine(std::string line);
  std::string writeToBed12 (bed12 *region);
  std::string writeHalfBed12(double start, double stop);

};

#endif
