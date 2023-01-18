/**
 * @file test_main.cpp
 * @author Robin Dowell 
 * @copyright Copyright 2022 Dowell Lab 
 * @brief This is for testing and debugging only.
 * @version 0.1
 * @date 2022-02-21
 * 
 */

#include "load.h"  // contains segment_fits class
#include "model.h" // contains fit2 implementation (classifier object)

#include "Bedfile.h"

int g_testing {};

/**
 * @brief Program main.  
 * @param argc
 * @param argv
 * @return 
 */
int main(int argc, char* argv[]) {

  /* Read in a segment of data */
  std::string data_file = "../test/examples/typical_region.bg";
  Bedgraph bg; 
  bg.load_file(data_file, 0);

  // Grab first interval for testing
  std::vector<bed4*> roi = bg.setRegions.regions[0];
  bed4 *testregion = roi[0];

  /* Attempt to fit a single model (K=1) */
}
