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

#include "Regions.h"    // SetROI, Segment

/**
 * @brief The complete contents of a bedfile.
 * Entries may be bed3, bed4, bed6 or bed12 
 */
class Bedfile {
public:
   std::string filename;       //!< file name 
   SetROI setRegions;        //!< Container for set of bed6 regions 
 
   // Constructors
	Bedfile();	// empty constructor

	/* FUNCTIONS: */
  // Debugging functions:
  std::string reportBedfileContents();        // used for debugging

  // Read the file
  void load_file (std::string); //!< Nuts and Bolts of this class, reads bed file.

  void write_file (std::string filename);  // inverse of load
};

/**
 * @brief The complete contents of a bedGraph.
 * We think of a bedgraph as regions of interest (similar to a bed file),
 * but with data.
 */
class Bedgraph: public Bedfile {
public:
  bool useExistingIntervals;       //<! has ROI from bed separately or not?

	// Constructors
  Bedgraph();

	/* FUNCTIONS: */
  std::string reportBedGraphContents();        // used for debugging

  // Read the file
  void load_file (std::string, bool); //!< Nuts and Bolts for reading a bedGraph file
  // void write_file (std::string filename);  // inverse of load

};

#endif
