/**
 * @file Intervals.cpp
 * @author Robin Dowell 
 * @brief Class code for genomic intervals
 * 
 * Current design: 
 * gIntervals maintain genomic coordinates and are the fundamental datatype
 * within the Interval trees (see ITree.cpp). Basic gInterval assumes BED4,
 * exists a bed6 class that extends for score and strand.
 * 
 * @version 0.1
 * @date 2022-01-27
 * 
 */
#include "Intervals.h"

#include <string>
#include <vector>
#include <iostream>

#include "split.h"
#include "Regions.h"

gInterval::gInterval() {
  chromosome = "NA";
  identifier = "empty";
  start = 0;
  stop = 0;
  data = NULL;
}

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

  data = NULL;
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
  // Because coordinates are half open, these are strickly less than.
  return ((start < i2->stop) && (i2->start < stop));
}

bool gInterval::Contains(double point) {
  // Note that stop is half open (e.g. strickly less than)
  return ((point < stop) && (point >= start));
}

/**
 * @brief Add a set of points (per bedGraph info) to this interval.
 * 
 */
void gInterval::addDataPoint(double st, double sp, double cov, bool expand) {
  // std::cout << "gInterval: " + to_string(st) + "," + to_string(sp) << std::endl;
  // If pointer to segment doesn't exist, create it.
  if (data == NULL) {
    data = new RawData(this);
  }   
  // Adjust for edge cases, expanding the region if permissible.
  double pt_edge = st;
  if (st < start) {
    if (expand) start = st;   // Adjust the gInterval
    pt_edge = start;    // Don't add non-overlapping points (edge case)
  }
  double pt_edge2 = sp;
  if (sp > stop) {
    if (expand) stop = sp;  // Adjust the gInterval
    pt_edge2 = st;    // Don't add non-overlapping points (edge case)
  }
  // std::cout << "ADDing: " + to_string(pt_edge) + "," + to_string(pt_edge2) << std::endl;
  data->addDataPoints(pt_edge, pt_edge2, cov);
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

  // Check that start < stop.  If not, throw an error?
  // Check start >= 0?
  start = stod(lineArray[1]);
  stop = stod(lineArray[2]);

  if (lineArray.size() >= 4) {  // At least a BED4, so identifier is present.
    identifier = lineArray[3];
  } else {  // bed3 files don't have an identifier.
    identifier = "";    // should we make up something internally?
  }
}

/********************  BED6 ***********************/

bed6::bed6():gInterval() {
  score = 0;
  strand = '.';
}


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
  score = v_score;  // Technically constrained to be between 0 and 1000 -- enforce?
  if (v_strand.compare(0,1,"+") == 0) {
    strand = '+';
  } else if(v_strand.compare(0,1, "-") == 0) {
    strand = '-';
  } else {
    strand = '.';
  }
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

