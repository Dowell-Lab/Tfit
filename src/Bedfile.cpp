/**
 * @file Bedfile.cpp
 * @author Robin Dowell 
 * @brief Load/store complete contents of a bed and bedgraph files 
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

#include "Intervals.h"    // gInterval, bed6
#include "Data.h"     // RawData, dInterval 
#include "Regions.h"    // SetROI, Segment
#include "ITree.h"
#include "split.h"

/**
 * @brief Construct a new Bedfile::Bedfile object
 */
Bedfile::Bedfile()
  : regions() {
  filename = "";
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

  if (FH){
    std::string line;   // We are going to read this file in one line at a time.
    int 	linenum = 0;    // line counter

    // Reading input file line by line
    while(getline(FH, line)){
      if (line.substr(0,1)!="#") { // ignore comment lines

      bed6 *iregion = new bed6();
      // Note that the bed6 object is doing sanity checking on the line.
      iregion->setfromBedLine(line);  // This interval's info.

      // Add region to collection.
      regions.addRegionToSet(iregion);

      } // for all lines in bedfile that aren't comments 
      linenum++;    // line counter, could be useful later.
    } // for each line in bedfile
  } else {  // filehandle error
    printf("couldn't open %s for reading\n", filename.c_str() );
    EXIT 	= true;
  }

  // Post file parsing setup steps.
  if (!EXIT) {
    regions.createSearchIndex();
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
   report += regions.write_out();
   return report;    
}

/*******************  BEDGRAPH ********************/


Bedgraph::Bedgraph()
  : Bedfile() {
  useExistingIntervals = 0; 
}

std::string Bedgraph::reportBedGraphContents() {
   // Summary should include: name of file:
   std::string report = filename;
   // And the rest of the object!!!
   return report;
}

/**
 * @brief Load a bedgraph file
 * 
 * @param v_filename   assumes a joint bedgraph file
 */
void Bedgraph::load_file(std::string v_filename) {
  filename = v_filename;
  ifstream FH(filename);

  bool EXIT 		= false;  // file error indicator
  // Variables needed temporarily
  vector<string> lineArray; // Contents of file, split on tab (\t) 
	string prevChrom="";	// What chrom was on the previous line?

  if (FH){
    std::string line;   // We are going to read this file in one line at a time.
    int 	linenum = 0;    // line counter

    // Reading input file line by line
    while(getline(FH, line)){
      if (line.substr(0,1)!="#") { // ignore comment lines
        // bedgraphs are always 4 column: chr start stop coverage 
        lineArray = string_split(line, '\t');
        if (lineArray.size() != 4) {
          EXIT = true;
          printf("\nLine number %d  in file %s was not formatted properly\nPlease see manual\n", linenum, filename.c_str());
          break;
        } else {
         if (useExistingIntervals) { // Add points to existing intervals.
           // Add this point to all the relevant segments
           regions.addDataToROI(lineArray[0], std::stod(lineArray[1]), 
                      std::stod(lineArray[2]),std::stod(lineArray[3]));
         } else { // There are no ROI, we are creating ROI as we go ...
           regions.addDataToSegments(lineArray[0], std::stod(lineArray[1]),
                        std::stod(lineArray[2]),std::stod(lineArray[3]));
         }
        }
      }
      linenum++;    // line counter
    } // for each line in bedfile
  } else {  // filehandle error
    printf("couldn't open %s for reading\n", filename.c_str() );
    EXIT 	= true;
  }

  // Post file parsing setup steps:
  if (!EXIT) {
    // Condition all data
  }
}
