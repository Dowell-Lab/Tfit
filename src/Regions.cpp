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
   // Currently only chromosome names reported.
   report = "\n" + chr_names.print_index_names();
   std::map<int, std::vector<gInterval *>>::iterator it;
   std::vector<gInterval *>::iterator intervals;
   for (it = regions.begin(); it != regions.end(); it++)  {
     for (intervals = it->second.begin(); intervals != it->second.end(); intervals++)  {
       report += "\n\t" + (*intervals)->write_out();
     }
   }
   return report;    
}

/**
 * @brief Add data ONLY to existing intervals in the collection
 * Used when adding data to roi, typically loaded in from a bed file
 * 
 */
bool SetROI::addDataToExistingROI(std::string chr, double start, double stop, double coverage) {
   std::vector<gInterval *> setToAdd;
   gInterval dataInterval(chr, start, stop, "temp");
   // Search existing CITree for overlapping intervals
   setToAdd = findOverlapIntervals(&dataInterval);

   if (setToAdd.size() == 0) {    // We found no overlapping intervals
      // std::cout << "NOT FOUND!" + dataInterval.write_out() << std::endl;
      return 0;
   }
   // std::cout << "FOUND!" + dataInterval.write_out() << std::endl;

   // Add this point to the overlapping gIntervals
   std::vector<gInterval *>::iterator it;
   for (it=setToAdd.begin(); it != setToAdd.end(); it++) {
     (*it)->addDataPoint(start,stop,coverage,0);
   }
   return 1;
}

void SetROI::ConditionDataSet(int v_delta, int v_scale) {
  std::map<int, std::vector<gInterval *>>::iterator mit;    // Outer Map iterator
  std::vector<gInterval *>::iterator  it;   //Inner vector iterator
  for (mit = regions.begin(); mit != regions.end(); mit++) {
    // For every index, lets go through the vector of gIntervals.
    for (it = mit->second.begin(); it != mit->second.end(); it++) {
      if ((*it)->data != NULL) {
        (*it)->data->cdata = new dInterval((*it)->data, v_delta, v_scale);
      }
    }
  }
}

/**
 * @brief Add data, building/expanding existing intervals/segments 
 * as necessary.   This is primarily used when adding data from a 
 * bedGraph without pre-existing ROI. 
 * 
 */
void SetROI::addDataCreateROI(std::string chr, double start, double stop, double coverage) {
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

