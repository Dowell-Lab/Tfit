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
#include "Bed.h"  // gInterval
#include "Data.h"  // gInterval

/**
 * @brief Collection of seeds (may have arisen by several possibilities)
 * 
 * Seeds are associated with a bed4, use coordinates relative
 * to the bed4 (e.g. bed4->start is position zero).
 * Seeds are constrained to reside within the bed4.
 * 
 * Can read/write to files in Bed12 format using the "exon starts" as
 * seed locations (size = 1) and "exon widths" as seed weights.
 */
class Seeds { 
  public:
  // Set of preferred seed points for mu with weights (in coverage position)
  std::vector<PointCov> mu_seeds;    

  // Constructor
  Seeds();

  //Functions
  std::string write_out();

  double getMaxWeight();
  double getMinWeight();
  int getNumSeeds();

  // Sort seeds by position or weight
  void SortByWeights();
  void SortByPositions();

  // User Input (i.e. I/O of seeds, perhaps bed12?)
  void getSeedsfromBedFields (std::string v_numSeeds, std::string v_weightsSeeds, 
                            std::string v_relativeSeedStarts);
  std::string writeSeedsAsBedFields();
};

/**
 * @brief How do we want to generate and manage seeds?
 * 
 */
class SeedManager { 
  public:
  dInterval *belongsTo;
  Seeds *setSeeds;  // convenience pointer

  // These could be private:
  Random numgen;    // Random number generator
  // std::vector<bool> used;    // Indicator, have we used this seed?

  // Constructor
  SeedManager();

  //Functions
  std::string write_out();

  void setupDataLink(dInterval *v_data);

  std::vector<PointCov> setupRandomSeeds(int numseeds, double v_lastposition);
  void weightRandomly(std::vector<PointCov> *seeds);
  void shuffleSeeds();
  double addUncertainty(PointCov *seed, double variance);
  std::vector<PointCov> grabSeedSet(int K);
};

#endif
