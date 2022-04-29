/**
 * @file Regions.cpp
 * @author Robin Dowell 
 * @brief Class code for Sets of ROI and segments 
 * 
 * @version 0.1
 * @date 2022-04-28
 * 
 */
#include "Regions.h"

#include <string>
#include <vector>
#include <iostream>

#include "split.h"
#include "Intervals.h"
#include "Data.h"

/**
 * @brief Construct a new SetROI::SetROI object
 * 
 */
SetROI::SetROI()
    : chr_names() {
  treesExist = 0;
}

/**
 * @brief Given a new region (bed6), add it to this set.
 * 
 * @param nregion 
 */
void SetROI::addRegionToSet(gInterval * nregion) {
  // be sure we have an identifier for this chromosome
  chr_names.addIdentifier(nregion->chromosome);

  // Now add the interval to the correct set.
  int idx = chr_names.lookupIndex(nregion->chromosome);
  regions[idx].push_back(nregion);
}

/**
 * @brief Creates a searchable centered Interval tree for each
 * chromosome
 */
void SetROI::createSearchIndex() {
    // setup interval trees, one per chromosome
    std::map<int, std::vector<gInterval *>>::iterator it;
    for (it = regions.begin(); it != regions.end(); it++) {
      searchable[it->first] = new CITree(it->second);
    }
    treesExist = 1;
}

/**
 * @brief find all the intervals that overlap a given input interval.
 * Since the ITree search work is done by CITree::overlapSearch 
 * the goal here is really just to do that on the right tree (i.e. using
 * the chromosome identifier).
 * 
 */
std::vector<gInterval *> SetROI::findOverlapIntervals(gInterval *input) {
  // std::cout << input->write_out() << std::endl;

  int index = chr_names.lookupIndex(input->chromosome);
  //  std::cout << "Index: " << index << std::endl;
  if (index >= 0) {  // exists in the bimap, search the right tree
    if (!treesExist) {
      createSearchIndex();
    }
    return searchable[index]->overlapSearch(input);
  } else {  // this index isn't found, return empty vector;
    std::vector<gInterval *> empty_vector;
    return empty_vector;
  }
}

/**
 * @brief Destroy the Interval Tree (clear up memory)
 * 
 */
void SetROI::clearTree() {
  searchable.clear(); // Removes Interval Trees to save space
  treesExist = 0;  
}

/**
 * @brief Debugging function for printing ITree for a chromosome
 * 
 * @param chromo 
 * @return std::string 
 */
std::string SetROI::print_tree_at_chromosome(std::string chromo) {
  int idx = chr_names.lookupIndex(chromo);
  if (treesExist) {
    return searchable[idx]->write_Full_Tree();
  } else {
    // Write vector of gIntervals for this chromosome
    return "";
  }
}

/**
 * @brief Debugging function for printing contents of a set.
 * 
 * @return std::string 
 */
std::string SetROI::write_out() {
   std::string report;
   // Then one per line:  chromosome name \t # intervals on that chromosome
   // Currently only chromosome names reported.
   std::map<int, CITree *>::iterator it;
   for (it = searchable.begin(); it != searchable.end(); it++) {
       report += "\nIndex: " + std::to_string(it->first) + " " + chr_names.lookupName(it->first);
          // + "\t" + std::to_string(searchable[it->first]->getSize());
   }
   return report;    
}

/**
 * @brief Add data ONLY to existing intervals in the collection
 * Used when adding data to roi, typically loaded in from a bed file
 * 
 */
void SetROI::addDataToROI(std::string chr, double start, double stop, double coverage) {
   std::vector<gInterval *> setToAdd;
   gInterval dataInterval(chr, start, stop, "temp");
   // Search existing CITree for overlapping intervals
   setToAdd = findOverlapIntervals(&dataInterval);

   // Add this point to those gIntervals
   std::vector<gInterval *>::iterator it;
   for (it=setToAdd.begin(); it != setToAdd.end(); it++) {
     (*it)->addDataPoint(start,stop,coverage,0);
   }
}

/**
 * @brief Add data, building/expanding existing intervals/segments 
 * as necessary.   This is primarily used when adding data from a 
 * bedGraph without pre-existing ROI. 
 * 
 */
void SetROI::addDataToSegments(std::string chr, double start, double stop, double coverage) {
   // Fetch index for this chromosome (adding if necessary)
  chr_names.addIdentifier(chr);
  int idx = chr_names.lookupIndex(chr);

  // If no gInterval exists, add to the regions list
  if (regions[idx].size() <= 0 ) { 
    regions[idx].push_back(new gInterval(chr, start, stop, "temp"));
  }
  // Add this point to the gInterval, expand if needed.
  regions[idx].front()->addDataPoint(start,stop,coverage,1);
}


/************  Segments ****************/

Segment::Segment() {
  genCoords = NULL;
  rawdata = NULL;
  cdata = NULL;
}

Segment::Segment(gInterval *v_genomicCoords) {
  genCoords = v_genomicCoords;
  cdata = NULL;
  rawdata = NULL;
}

std::string Segment::write_out() {
  std::string output = genCoords->write_out();
  output += "\n" + rawdata->write_out(); 
  return output;
}

/**
 * @brief Add data points (a line from bedGraph) to
 * this region. 
 */
void Segment::addDataPoints(double st, double sp, double cov) {
  if (rawdata == NULL) { // First point
    rawdata = new RawData(); 
    rawdata->belongsTo = this->genCoords;   // Link data to gInterval
  } 
  rawdata->addDataPoints(st, sp, cov);
}

void Segment::ConditionData(int v_delta, int v_scale) {
  if (rawdata == NULL) { return; }
  cdata = new dInterval(rawdata, v_delta, v_scale);
}
