/**
 * @file Bedfile.h
 * @author Robin Dowell
 * @brief Header file for Bedfile and Bedgraph classes
 * @version 0.1
 */
#ifndef Bedfile_H 
#define Bedfile_H 

#include <map>
#include <string>
#include <vector>

#include "Intervals.h"      // gInterval, bed6
#include "Data.h"       // RawData, dInterval
#include "Regions.h"    // SetROI, Segment
#include "ITree.h"
#include "helper.h"


/**
 * @brief The complete contents of a bedfile.
 * Entries may be bed3, bed4, or bed6 (or more, but ignored if > 6).
 * 
 * Stores in sets if ITrees, one per chromosome, also provides methods
 * for looking up intervals (i.e. with chromosome identifier).
 */
class Bedfile {
public:
   std::string filename;       // read in file name
   SetROI setRegions;        // Container for set of bed6  regions 
 
   // Constructors
	Bedfile();	// empty constructor

	/* FUNCTIONS: */
  // Debugging functions:
  std::string reportBedfileContents();        // used for debugging

  // Read the file
  void load_file (std::string);
};

/**
 * @brief The complete contents of a bedGraph.
 * We think of a bedgraph as regions of interest (similar to a bed file),
 * but with data.
 * 
 * Thus a large part of the responsibility of this class is to ensure
 * that data exists for each gInterval. 
 * 
 */
class Bedgraph: public Bedfile {
public:
  bool useExistingIntervals;       // has ROI from bed separately or not?

	// Constructors
  Bedgraph();

	/* FUNCTIONS: */
  std::string reportBedGraphContents();        // used for debugging

  // Read the file
  void load_file (std::string, bool);

};

#endif
