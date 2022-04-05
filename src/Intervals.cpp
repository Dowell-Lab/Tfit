/**
 * @file Intervals.cpp
 * @author Robin Dowell 
 * @brief Class code for data intervals
 * Two kinds of intervals:  genomic (gIntervals) and data (dIntervals)
 * 
 * Current design: 
 * gIntervals maintain genomic coordinates and are the fundamental datatype
 * within the Interval trees (see ITree.cpp). Basic gInterval assumes BED4,
 * exists a bed6 class that extends for score and strand.
 * 
 * dIntervals contain data (two strands) per an interval but do so in zero
 * based coordinates.  Can translate back to gInterval if correspondance is setup.
 * 
 * Data intervals must have rapid data access for EM algorithm.
 * @version 0.1
 * @date 2022-01-27
 * 
 */
#include "Intervals.h"

#include <string>
#include <vector>
#include <iostream>

#include "split.h"

/**
 * @brief Construct a new g Interval::g Interval object
 * 
 * Representation is as in BED files:
 *  thus expects 0-based half open coordinates
 * 
 * @param v_chromosome   Chromosome name (field 1 of BED)
 * @param v_start        Start position (field 2 of BED)
 * @param v_stop         Stop position (field 3 of BED)
 * @param v_identifier   Identifier (field 4 of BED)
 */
gInterval::gInterval(std::string v_chromosome, double v_start, double v_stop, std::string v_identifier) {
  chromosome = v_chromosome;
  start = v_start;
  stop = v_stop;
  // Should we check that start < stop?  Throw an error when it doesn't?
  identifier = v_identifier;
}

gInterval::gInterval() {
  chromosome = "NA";
  identifier = "empty";
  start = 0;
  stop = 0;
}

/**
 * @brief This is a standard content dump function.  Primarily
 * used in debugging.  Notice the string does NOT end in a newline.
 * 
 * @return std::string 
 */
std::string gInterval::write_out() {
  std::string text = ("#" + chromosome + ":" + std::to_string((int)start) + "-" 
		+ std::to_string((int)stop) + "," + identifier);
  return text;
}

/**
 * @brief Writes out the interval as a line in a BED file.  
 * Note does NOT include newline.
 * 
 * @return std::string 
 */
std::string gInterval::write_asBED() {
  std::string text = (chromosome + "\t" + std::to_string((int)start) + "\t" 
		+ std::to_string((int)stop) + "\t" + identifier);
  return text;
}

/**
 * @brief Set this gInterval based on a single line from a BED4 file.
 * 
 * @param line expects a single line of a BED file 
 */
void gInterval::setfromBedLine(std::string line) {
  std::vector<std::string> lineArray; // Contents of line, split on tab (\t) 
  lineArray=string_split(line, '\t');
  setBED4fromStrings(lineArray);
}

bool gInterval::Overlap(gInterval *i2) {
   return ((start <= i2->stop) && (i2->start <= stop));
}

bool gInterval::Contains(double point) {
  return ((point <= stop) && (point >= start));
}

/**
 * @brief Helper function that sets the basic BED4 contents from a vector of strings.
 * 
 * @param lineArray expects: chromosome start stop ID in lineArray and ignores rest.
 */
void gInterval::setBED4fromStrings(std::vector<std::string> lineArray) {
  if (lineArray.size() < 3) {  // accept BED3 or BED4
    // Should this throw an exception?
    cout << "ERROR: Invalid BED file\n" << std::endl;
    chromosome = "formaterror";   // hacking an internal indicator of error for now.
  }
  chromosome = lineArray[0];
  start = stod(lineArray[1]);
  stop = stod(lineArray[2]);

  if (lineArray.size() >= 4) {  // At least a BED4, so identifier is present.
    identifier = lineArray[3];
  } else {
    identifier = "";    // identifier not present in BED3 files; should use internal
  }
}

/********************  BED6 ***********************/

/**
 * @brief Construct a new bed6::bed6 object
 * 
 * Representation is as in BED6 files:
 *  thus expects 0-based half open coordinates
 * 
 * @param v_chromosome   Chromosome name (field 1 of BED)
 * @param v_start        Start position (field 2 of BED)
 * @param v_stop         Stop position (field 3 of BED)
 * @param v_identifier   Identifier (field 4 of BED)
 * @param v_score        Score (field 5 of BED6)
 * @param v_strand       Strand info (field 6 of BED6); read as string, stored as char
 */
bed6::bed6(std::string v_chromosome, double v_start, double v_stop, 
    std::string v_identifier, int v_score, std::string v_strand) 
    : gInterval(v_chromosome, v_start, v_stop, v_identifier) {
  score = v_score;
  if (v_strand.compare(0,1,"+") == 0) {
    strand = '+';
  } else if(v_strand.compare(0,1, "-") == 0) {
    strand = '-';
  } else {
    strand = '.';
  }
}

bed6::bed6():gInterval() {
  score = 0;
  strand = '.';
}

/**
 * @brief This is a standard content dump function.  Primarily
 * used in debugging.  Notice the string does NOT end in a newline.
 * 
 * @return std::string 
 */
std::string bed6::write_out() {
  std::string text = gInterval::write_out();
  text += "," + std::to_string(score) + "," + strand;
  return text;
}

/**
 * @brief Writes out the interval as a line in a BED6 file.  
 * Note does NOT include newline.
 * 
 * @return std::string 
 */
std::string bed6::write_asBED() {
  std::string text = gInterval::write_asBED();
  text += "\t" + std::to_string(score) + "\t" + strand;
  return text;
}

/**
 * @brief Set this bed6 based on a single line from a BED file.
 * 
 * @param line expects a single line of a BED file 
 */
void bed6::setfromBedLine(std::string line) {
  std::vector<std::string> lineArray; // Contents of line, split on tab (\t) 
  lineArray=string_split(line, '\t');
  setBED4fromStrings(lineArray);  // setup basic BED3/4 contents.

  if (lineArray.size() >= 6) {  // At least a BED6, so score and strand are present.
    score = stod(lineArray[4]);
    if (lineArray[5].compare(0,1,"+") == 0) {
      strand = '+';
    } else if(lineArray[5].compare(0,1,"-") == 0) {
      strand = '-';
    } else {
      strand = '.';
    }
  } else {
    // score and strand unavailable in BED3 and BED4
    // Not currently throwing error, but should we?
    score = -1;   // using -1 to indicate unavailable (i.e. wrong file type?)
    strand = '.';
  }
}

/****************** dInterval *********************/

/**
 * @brief Constructors: dInterval class
 * @author Robin Dowell
 *
 * Purpose: create/allocate instances of a data Interval 
 *
 * The data Interval class contains both strands of data associated
 * with a particular region.  There is an empty constructor option. 
 *
 * @param identifier 
 * @param st   Start
 * @param sp   Stop
 * @param Integer identifier (opt)
 * @param STR  Strand (as string) (opt)
 *
 * @bug STR should be a char with only '.' '+' and '-' as valid.
 *
 */
dInterval::dInterval(std::string v_identifier, int v_min, int v_max) {
  ID = v_identifier;
  minX = v_min;
  maxX = v_max;
  XN = 1;
  SCALE = 1;

  N = 0;
  fN = 0;
  rN = 0;
}

// empty constructor
dInterval::dInterval() {
  ID = "empty";
  minX = 0;
  maxX = 0;
  XN = 1;
  SCALE = 1;
  N = 0;
  fN = 0;
  rN = 0;
}

/**
 * @brief This is a standard content dump function.  Primarily
 * used in debugging.  Notice the string does NOT end in a newline.
 * 
 * @return std::string 
 */
std::string dInterval::write_out() {
  std::string text = ("#" + ID + ":d:" + std::to_string((int)minX) + "-" 
		+ std::to_string((int)maxX) + "," + std::to_string((int)XN) + ","
    + std::to_string((int)N));
  return text;
}


/**
 * @brief The number of data points in the interval.  
 * This could be nucleotides (smallest unit) to bins (groups of nts)
 * to the number of (possibly random) points along the interval.
 * 
 * @return double  Size of the interval in usable steps.
 */
double dInterval::num_elements() {
  return XN;
}
/**
 * @brief Data from forward strand at xth index.
 * Note a "unit" here could be nucleotides, bins (groups of nts),
 * or points along the interval.
 * 
 * @arg x The index of an element of the forward strand data.
 * Note that this is not a position (genomic or scaled).
 * 
 * @return double 
 */
double dInterval::forward(int x) {
  return X[1][x];
}
/**
 * @brief Data from reverse strand at xth index.
 * Note a "unit" here could be nucleotides, bins (groups of nts),
 * or points along the interval.
 * 
 * @arg x The index an element of the reverse strand data.
 * Note that this is not a position (genomic or scaled).
 * 
 * @return double 
 */
double dInterval::reverse(int x) {
  return X[2][x];
}

/**
 * @brief Position at the xth index.
 * 
 * @arg x The index of a data element.
 * Note that this is not a position (genomic or scaled).
 * 
 * @return double 
 */
double dInterval::position(int x) {
  return X[0][x];
}
  
double dInterval::sum_Region() {
  return N;
}
double dInterval::sum_forward() {
  return fN;
}
double dInterval::sum_reverse() {
  return rN;
}
