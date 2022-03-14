/**
 * @file Intervals.cpp
 * @author Robin Dowell 
 * @brief Class code for data intervals
 * @version 0.1
 * @date 2022-01-27
 * 
 */
#include "dInterval.h"

#include <string>
#include <vector>

/**
 * @brief Constructors: dInterval class
 * @author Robin Dowell
 *
 * Purpose: create/allocate instances of a data Interval 
 *
 * The data Interval class contains both strands of data associated
 * with a particular region.  There is an empty constructor option. 
 *
 * @param chr  Chromosome 
 * @param st   Start
 * @param sp   Stop
 * @param Integer identifier (opt)
 * @param STR  Strand (as string) (opt)
 *
 * @bug STR should be a char with only '.' '+' and '-' as valid.
 *
 */
dInterval::dInterval(std::string name, int start, int stop) {
  ID = name;
  minX = start;
  maxX = stop;
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
  std::string text = ("#" + ID + ":" + std::to_string((int)minX) + "-" 
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
