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
#include "Bed.h"   // bed12
#include "Data.h"   // RawData, Segment 
#include "ITree.h"    // CITree

/**
 * @brief This is the genomic contents of regions.
 * With no pre-defined ROI, this ends up being one bed4
 * per chromosome/contig. 
 * 
 * If created from a bed file, then ROI are read directly from the file.
 * If created from a bedGraph, then ROI must be "made" and adjusted to
 * include all points in the bedGraph.
 * 
 */
class SetROI {
public:
  Bimap chr_names;     // system for converting names to indexes (& vice versa)
  std::map<int, std::vector<bed12 *>> regions;  // collection, one per "chr"

  // The Interval tree is really only needed when searching.  
  std::map<int, CITree *> searchable;   // index to interval tree 
  bool treesExist;    // Indicator for Interval tree built 

  // Constructors
  SetROI();
  // Destructor
  ~SetROI();

  // Functions
  void addRegionToSet(bed12 *);
  std::string print_tree_at_chromosome(std::string chromo);
  std::string write_out();
   
  // Find overlapping intervals (search)
  void createSearchIndex();   // Builds the interval trees given the regions 
  std::vector<bed12 *>findOverlapIntervals(bed12 *);
  void clearTrees();

  // Adding data to the set of regions (with vs without roi)
  void addDataCreateROI(std::string chr, double start, double stop, double coverage);
  bool addDataToExistingROI(std::string chr, double start, double stop, double coverage);
  void ConditionDataSet(int v_delta, int v_scale);
};


#endif
