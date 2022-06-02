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
#include "ModelSets.h" // ModelContainer
#include "Data.h"		// dInterval

/**
 * @brief How do we want to generate and manage seeds?
 * 
 * Seeds should probably be associated with a gInterval.
 * The seed coordinates should be constrained to reside within the gInterval.
 * 
 * As bed12:  bed6 is the interval: chr start stop identifier score strand
 * field 7  thick_start   typically start of CDS or repeats start
 * field 8  thick_end     typically end of CDS or repeats stop
 * field 9  RGB value     0,0,0 formatted
 * field 10 blockCount (#exons)
 * field 11 blockSizes
 * field 12 blockStarts   positions relative to start, e.g. first is zero
 */
class Seeds { 
  public:
  // set of seeds (read from file or from alg)
  // constraints on seeds (protection_zone, others)
  
  // Constructor
  Seeds();

  //Functions
  std::string write_out();

  // Grab a seed (mark it as used or remove)

  // Randomly seed: draw from a range
  // Run a seeding algorithm (template matching or other)
  // User Input (i.e. I/O of seeds, perhaps bed12?)

};

#endif
