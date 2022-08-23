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
#include "EMseeds.h"

gInterval::gInterval() {
  chromosome = "NA";
  identifier = "empty";
  start = 0;
  stop = 0;
  data = NULL;
}

/**
 * @brief Construct a new gInterval::gInterval object
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
  // Note: prettyDecimal with -1 sigfigs removes the decimal!
  std::string text = ("#" + chromosome + ":" + tfit::prettyDecimal(start,-1) + "-" 
		+ tfit::prettyDecimal(stop,-1) + "," + identifier);
  return text;
}

/**
 * @brief Writes out the interval as a line in a BED file.  
 * Note does NOT include newline.
 * 
 * @return std::string 
 */
std::string gInterval::write_asBEDline() {
  std::string text = (chromosome + "\t" + tfit::prettyDecimal(start,-1) + "\t" 
		+ tfit::prettyDecimal(stop,-1) + "\t" + identifier);
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
  setBEDfromStrings(lineArray);
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
void gInterval::addDataPoint(double v_start, double v_stop, double cov, bool expand) {
  // If pointer to segment doesn't exist, create it.
  if (data == NULL) {
    data = new RawData(this);
  }   
  // Adjust for edge cases, expanding the region if permissible.
  if (expand && (v_start < start)) {
    start = v_start;   // Adjust the gInterval
  }
  if (expand && (v_stop > stop)) {
    stop = v_stop;  // Adjust the gInterval
  }
  // Adding these data points to the set
  data->addDataPoints(v_start, v_stop, cov);
}

/**
 * @brief Helper function that sets the basic BED4 contents from a vector of strings.
 * 
 * @param lineArray expects: chromosome start stop ID in lineArray and ignores rest.
 */
void gInterval::setBEDfromStrings(std::vector<std::string> lineArray) {
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
  text += "," + tfit::prettyDecimal(score,4) + "," + strand;
  return text;
}

/**
 * @brief Writes out the interval as a line in a BED6 file.  
 * Note does NOT include newline.
 * 
 * @return std::string 
 */
std::string bed6::write_asBEDline() {
  std::string text = gInterval::write_asBEDline();
  text += "\t" + tfit::prettyDecimal(score, 4) + "\t" + strand;
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
  setBEDfromStrings(lineArray);
}

void bed6::setBEDfromStrings(std::vector<std::string> lineArray) {
  gInterval::setBEDfromStrings(lineArray);  // setup basic BED3/4 contents.

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


/********************  BED12 ***********************/

bed12::bed12(): bed6() {
  seeds = NULL;
}

bed12::bed12(std::string v_chromosome, double v_start, double v_stop, 
    std::string v_identifier, int v_score, std::string v_strand, 
    std::string v_nextsix) : 
    bed6(v_chromosome, v_start, v_stop, v_identifier,v_score, v_strand) {
      // Allocate a Seeds object to fill here.
    Seeds seedinfo;
    seeds = &seedinfo;
    setfromLastSix(v_nextsix);
}

/**
 * @brief This is a standard content dump function.  Primarily
 * used in debugging.  Notice the string does NOT end in a newline.
 * 
 * @return std::string 
 */
std::string bed12::write_out() {
  std::string text = bed6::write_out();
  if (seeds != NULL) { text += "," + seeds->write_out(); }
  return text;
}

/**
 * @brief Writes out the interval as a line in a BED12 file.  
 * Note does NOT include newline.
 * 
 * As bed12:  bed6 is the interval: chr start stop identifier score strand
 * field 7  thick_start   typically start of CDS or repeats start
 * field 8  thick_end     typically end of CDS or repeats stop
 * field 9  RGB value     0,0,0 formatted
 * field 10 blockCount (#exons)
 * field 11 blockSizes
 * field 12 blockStarts   positions relative to start, e.g. first is zero
 * 
 * Using blockStarts as seed locations and blockSizes as relative weighting.
*/
std::string bed12::write_asBEDline() {
  std::string text = bed6::write_asBEDline();
  // Need fields 7-9!!
  if (seeds != NULL) text += "\t" + seeds->writeSeedsAsBedFields(); // last 3 fields
  return text;
}

/**
 * @brief Set this bed6 based on a single line from a BED file.
 * 
 * @param line expects a single line of a BED file 
 */
void bed12::setfromBedLine(std::string v_line) {
  std::vector<std::string> lineArray; // Contents of line, split on tab (\t) 
  lineArray=string_split(v_line, '\t');
  bed6::setBEDfromStrings(lineArray);
  // Confirm BED12?!?

  std::string numSeeds = lineArray[9];
  std::string seedWeights = lineArray[10];
  std::string seedStarts = lineArray[11];
  seeds->getSeedsfromBedFields(numSeeds, seedWeights, seedStarts);
}

/**
 * @brief  Used for creating seeds from a string (half bed12 line)
 * 
 * @param v_lastsix Assumes does NOT start with \t
 */
void bed12::setfromLastSix(std::string v_lastsix) {
  std::vector<std::string> lineArray; // Contents of line, split on tab (\t) 
  lineArray=string_split(v_lastsix, '\t');

  // Note here we expect half the string so indexs shift.
  std::string numSeeds = lineArray[3];
  std::string seedWeights = lineArray[4];
  std::string seedStarts = lineArray[5];
  seeds->getSeedsfromBedFields(numSeeds, seedWeights, seedStarts);
}

