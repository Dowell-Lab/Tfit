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
  num_chr = 0;
}
Bedfile::Bedfile(std::string FILE) {
  filename = FILE; 
  num_chr = 0;
}

/**
 * @brief Adding a new identifier (chromosome name) to the container
 * 
 * @param chrid  // string: name of current chromosome
 */
void Bedfile::addChromosome(std::string chrid) {
   if (!(chr2index.count(chrid))) {  // a new chromosome
    // Lets add a new chromosome to the index list
    chr2index[chrid] = num_chr; 
    // Reverse index
    IDindex[num_chr] = chrid;
    num_chr++;
   }
}

/**
 * @brief  load a bedfile of intervals
 * @author Robin Dowell 
 * @param FILE name of bedfile containing intervals (example: singleregion.bed)
 * @param spec_chrom  string spec_chrom 	= P->p["-chr"];
 * @param pad   int pad = stoi(P->p["-pad"])+1;
 */
void Bedfile::load_file() {
  ifstream FH(filename);

  bool EXIT 		= false;  // file error indicator
  // collection of gIntervals from file, one "set" per chromosome ID
  // std::map<std::string, std::vector<bed6 *>> regions;  
  map<int, std::vector<gInterval *>> regions;  

  if (FH){
    std::string line;   // We are going to read this file in one line at a time.
    int 	i = 0;    // line counter

    // Reading input file line by line
    while(getline(FH, line)){
      if (line.substr(0,1)!="#") { // ignore comment lines

      bed6 interval;  // This interval's info.
      interval.setfromBedLine(line);

      // Is the current chromosome one we've seen before?
      if (!(chr2index.count(interval.chromosome))) {
        addChromosome(interval.chromosome);
      }
      // Now add the interval to the correct set.
      int idx = chr2index[interval.chromosome];
      regions[idx].push_back(&interval);

      } // not a comment
      i++;    // line counter
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

