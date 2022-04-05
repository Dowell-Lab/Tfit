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
Bedfile::Bedfile(std::string FILE) {
  filename = FILE; 
}

/**
 * @brief  load a bedfile of intervals
 * @author Robin Dowell 
 * @param FILE name of bedfile containing intervals (example: singleregion.bed)
 * @param spec_chrom  string spec_chrom 	= P->p["-chr"];
 * @param pad   int pad = stoi(P->p["-pad"])+1;
 */
void Bedfile::load_file (std::string spec_chrom, int pad, bool center) {
  bool debug = false;    // a debugging indicator
  ifstream FH(filename);

  // if (debug) { printf("\n  spec_chrom: %s, pad: %d\n", spec_chrom.c_str(), pad);}
  bool EXIT 		= false;

  if (FH){
    std::string line;   // We are going to read this file in one line at a time.
    int 	i = 0;
    bool PASSED 	= true; // check on identifier.

    // These are the relevant Bed line fields:
    std::string chrom;  // field[0]
    int start, stop;    // field[1] and field[2]
    // field[3] is the identifier, will make one internally if not provided.
    std::string strand; // field[5] will set to "." if not provided

    // Reading input file line by line
    while(getline(FH, line)){
      // if (debug) { printf("  %s\n", line.c_str()); }

      if (line.substr(0,1)!="#") { // ignore comment lines

      // Should have the objects parse the line.  But do I keep bed4 and bed6 
      // objects distinctly?  Should do check for "proper bed"?

      // Should objects be padded?  Centered?

      } // not a comment
    } // for each line in bedfile
  }else{  // filehandle error
    printf("couldn't open %s for reading\n", filename.c_str() );
    EXIT 	= true;
  }

  // setup interval trees, one per chromosome
  // setup indexes to each tree

}

