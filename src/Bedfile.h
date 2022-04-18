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

#include "Intervals.h"
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
   std::map<int, CITree *> intervals;   // index to interval tree 
   Bimap chr_names;     // system for converting names to indexes (& vice versa)

	// Constructors
	Bedfile();	// empty constructor

	/* FUNCTIONS: */
    // Debugging functions:
    std::string print_tree_at_chromosome(std::string);
    std::string reportBedfileContents();        // used for debugging

    // Read the file
    void load_file (std::string);

    // Find overlapping intervals (search)
    std::vector<gInterval *>findOverlapIntervals(gInterval *);

};

/**
 * @brief The complete contents of a bedGraph.
 * 
 */
class Bedgraph{
public:
   std::string filename;       // read in file name
   Bimap chr_names;     // system for converting names to indexes (& vice versa)
   std::map<int, CITree *> intervals;   // index (to segments?!?!)

	// Constructors
	Bedgraph();	// empty constructor

	/* FUNCTIONS: */
    // Debugging functions:
    std::string reportBedGraphContents();        // used for debugging

    // Read the file
    void load_file (std::string); 

    // Find overlapping intervals (search)

};

#endif
