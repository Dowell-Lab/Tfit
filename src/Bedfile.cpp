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
      iregion->setfromBedLine(line);  // This interval's info.

      // be sure we have an identifier for this chromosome
      chr_names.addIdentifier(iregion->chromosome);

      // Now add the interval to the correct set.
      int idx = chr_names.lookupIndex(iregion->chromosome);

      // std::cout <<  "New interval: " << idx << " " << iregion.write_out() << std::endl;
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

