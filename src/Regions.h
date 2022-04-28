/**
 * @file Regions.h
 * @author Robin Dowell
 * @brief Header file for regions classes (Segments and SetROI)
 *
 * @version 0.1
 */
#ifndef Regions_H 
#define Regions_H 

#include <string>
#include <vector>

#include "helper.h"   // Bimap
#include "Intervals.h"   // gInterval, bed6
#include "Data.h"   // RawData, Segment 
#include "ITree.h"    // CITree

/**
 * @brief This is the genomic contents of regions.
 * With no pre-defined ROI, this ends up being one gInterval
 * per chromosome/contig. 
 * 
 * Thus a bed file defines regions of interest (subsets of coords)
 * and bedGraphs (et. al.) have implicit the whole chromosome/contig as a gInterval.
 */
class SetROI {
  public:
   Bimap chr_names;     // system for converting names to indexes (& vice versa)
   std::map<int, std::vector<gInterval *>> regions;  // Unstructured collection

   // The Interval tree is really only needed when searching.  
   std::map<int, CITree *> searchable;   // index to interval tree 
   int treesExist;    // Indicator for Interval tree built 

   // Constructors
   SetROI();

  // Functions
  void addRegionToSet(gInterval *);
  std::string print_tree_at_chromosome(std::string chromo);
  std::string write_out();
   
  // Find overlapping intervals (search)
  void createSearchIndex();   // Builds the interval trees given the regions 
  std::vector<gInterval *>findOverlapIntervals(gInterval *);
  void clearTree();

};

/**
 * @brief A segment is a portion of the genome with data.
 * Cross references a gInterval to a dInterval
 * 
 */
class Segment {
public:
  bed6  *genCoords;
  dInterval  *data;

  // Constructors
  Segment();
  Segment(std::string); // from just a chr name 
  // Segment(bed6);  // from an existing roi (bed6 assumed)

  /* FUCTIONS: */
  // Reporting out 
  std::string write_out();

};


#endif
