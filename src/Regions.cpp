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


/************  Segments ****************/

Segment::Segment()
  : genCoords(),
    data() {
  
}

Segment::Segment(std::string v_chr)
  : Segment() {
 genCoords->chromosome = v_chr;
}

std::string Segment::write_out() {
  std::string output = genCoords->write_out();
  output += "\n" + data->write_out(); 
  return output;
}

