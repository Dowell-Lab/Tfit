/**
 * @file Bedfile.cpp
 * @author Robin Dowell 
 * @brief Complete contents of a bedfile 
 * @version 0.1
 * @date 2022-04-05
 * 
 */
#include "Bedfile.h"

#include <stdio.h>

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "Intervals.h"
#include "ITree.h"
#include "split.h"

/**
 * @brief Construct a new Bedfile::Bedfile object
 */
Bedfile::Bedfile() {
  filename = "";
}

/**
 * @brief Debugging function for printing ITree for a chromosome
 * 
 * @param chromo 
 * @return std::string 
 */
std::string Bedfile::print_tree_at_chromosome(std::string chromo) {
    int idx = chr_names.lookupIndex(chromo);
    return intervals[idx]->write_Full_Tree();
}

/**
 * @brief  load a bedfile of intervals
 * @author Robin Dowell 
 * @param fle  name of bedfile containing intervals (example: singleregion.bed)
 */
void Bedfile::load_file(std::string file) {
  filename = file;
  ifstream FH(filename);

  bool EXIT 		= false;  // file error indicator

  // collection of gIntervals from file, one "set" per chromosome ID
  std::map<int, std::vector<gInterval *>> regions;  

  if (FH){
    std::string line;   // We are going to read this file in one line at a time.
    int 	i = 0;    // line counter

    // Reading input file line by line
    while(getline(FH, line)){
      if (line.substr(0,1)!="#") { // ignore comment lines

      bed6 *iregion = new bed6();
      // Note that the bed6 object is doing sanity checking on the line.
      iregion->setfromBedLine(line);  // This interval's info.


      // be sure we have an identifier for this chromosome
      chr_names.addIdentifier(iregion->chromosome);

      // Now add the interval to the correct set.
      int idx = chr_names.lookupIndex(iregion->chromosome);

      regions[idx].push_back(iregion);

      } // for all lines in bedfile that aren't comments 
      i++;    // line counter, could be useful later.
    } // for each line in bedfile
  } else {  // filehandle error
    printf("couldn't open %s for reading\n", filename.c_str() );
    EXIT 	= true;
  }

   if (!EXIT) {
    // setup interval trees, one per chromosome
    std::map<int, std::vector<gInterval *>>::iterator it;
    for (it = regions.begin(); it != regions.end(); it++) {
      intervals[it->first] = new CITree(it->second);
    }
   }
}

/**
 * @brief find all the intervals that overlap a given input interval.
 * Since the ITree search work is done by CITree::overlapSearch 
 * the goal here is really just to do that on the right tree (i.e. using
 * the chromosome identifier).
 * 
 */
std::vector<gInterval *> Bedfile::findOverlapIntervals(gInterval *input) {
  int index = chr_names.lookupIndex(input->chromosome);
  //  std::cout << "Index: " << index << std::endl;
  if (index >= 0) {  // exists in the bimap, search the right tree
    return intervals[index]->overlapSearch(input);
  } else {  // this index isn't found, return empty vector;
    std::vector<gInterval *> empty_vector;
    return empty_vector;
  }
}

/**
 * @brief Provide a basic report on a bedfile.
 * 
 * @return std::string 
 */
std::string Bedfile::reportBedfileContents() {
   // Summary should include: name of file:
   std::string report = filename;
   // Then one per line:  chromosome name \t # intervals on that chromosome
   // Currently only chromosome names reported.
   std::map<int, CITree *>::iterator it;
   for (it = intervals.begin(); it != intervals.end(); it++) {
       report += "\nIndex: " + std::to_string(it->first) + " " + chr_names.lookupName(it->first);
          // + "\t" + std::to_string(intervals[it->first]->getSize());
   }
   return report;    
}

