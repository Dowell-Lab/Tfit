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
  : setRegions() {
  filename = "";
}

/**
 * @brief  load a bedfile of intervals
 * @author Robin Dowell 
 * @param file  name of bedfile containing intervals (example: singleregion.bed)
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

      bed12 *iregion = new bed12();
      // Note that the bed6 object is doing sanity checking on the line.
      iregion->setfromBedLine(line);  // This interval's info.

      // Add region to collection.
      setRegions.addRegionToSet(iregion);

      } // for all lines in bedfile that aren't comments 
      linenum++;    // line counter, could be useful later.
    } // for each line in bedfile
  } else {  // filehandle error
    printf("couldn't open %s for reading\n", filename.c_str() );
    EXIT 	= true;
  }

  // Post file parsing setup steps.
  if (!EXIT) {
    setRegions.createSearchIndex();
  }
}

/**
 * @brief Write out a bedfile contents.
 * 
 * Q: Do we test if the filename exists first?
 * Q: Should we *sort* the intervals by start so that we get a consistent
 *   output ordering?
 * 
 * @param filename 
 */
void Bedfile::write_file (std::string filename) {
  ofstream FH(filename);    // Test for existence before using?
  
  // For every chromosome in setRegions
  std::map<int, std::vector<bed12 *>>::iterator it;
  std::vector<bed12 *>::iterator interval;
  for (it = setRegions.regions.begin(); it != setRegions.regions.end(); it++)  {
    // For every region with a given chromosome
    for (interval = it->second.begin(); interval != it->second.end(); interval++)  {
      // Output the region in appropriate bed format (3, 4, 6, or 12)
       FH << (*interval)->write_asBEDline() << std::endl;
    }
  }
  FH.close();
}



/**
 * @brief Provide a basic report on a bed file (debugging).
 * 
 * @return std::string  printable string of object contents 
 */
std::string Bedfile::reportBedfileContents() {
   // Summary should include: name of file:
   std::string report = filename;
   report += setRegions.write_out();
   return report;    
}

/*******************  BEDGRAPH ********************/

Bedgraph::Bedgraph()
  : Bedfile() {
  useExistingIntervals = 0; 
}

/**
 * @brief Provide a basic report on a bedGraph file.
 * 
 * @return std::string  printable string of object contents 
 */
std::string Bedgraph::reportBedGraphContents() {
   // Summary should include: name of file:
   std::string report;
   if (useExistingIntervals) {
     report += "Use: TRUE";
   } else {
     report  += "Use: FALSE";
   }
   // And the rest of the object!!!
   report += reportBedfileContents();
   return report;
}

/**
 * @brief Load a bedgraph file
 * 
 * @param v_filename   assumes a joint bedgraph file
 */
void Bedgraph::load_file(std::string v_filename, bool useExisting) {
  filename = v_filename;
  ifstream FH(filename);

  useExistingIntervals = useExisting;
  bool EXIT 		= false;  // file error indicator
  // Variables needed temporarily
  vector<string> lineArray; // Contents of file, split on tab (\t) 
	string prevChrom="";	// What chrom was on the previous line?

  // cout << "In Load File" << std::endl;

  if (FH){
    std::string line;   // We are going to read this file in one line at a time.
    int 	linenum = 0;    // line counter

    // Reading input file line by line
    while(getline(FH, line)){
      // cout << "reading line " + std::to_string(linenum) << std::endl;
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
           setRegions.addDataToExistingROI(lineArray[0], std::stod(lineArray[1]), 
                      std::stod(lineArray[2]),std::stod(lineArray[3]));
         } else { // There are no ROI, we are creating ROI as we go ...
           setRegions.addDataCreateROI(lineArray[0], std::stod(lineArray[1]),
                        std::stod(lineArray[2]),std::stod(lineArray[3]));
         }
        }
      }
      linenum++;    // line counter
        // cout << "lets start " + std::to_string(linenum) << std::endl;
    } // for each line in bedfile
  } else {  // filehandle error
    printf("couldn't open %s for reading\n", filename.c_str() );
    EXIT 	= true;
  }

  // Post file parsing setup steps:
  int delta = 10; int scale = 1;  // These should come from params;
  if (!EXIT) {
    // Condition all data
    setRegions.ConditionDataSet(delta, scale);
  }
}
