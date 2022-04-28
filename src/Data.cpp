/**
 * @file Data.cpp
 * @author Robin Dowell 
 * @brief Class code for data intervals
 * 
 * Current design: 
 * dIntervals contain data (two strands) per an interval but do so in zero
 * based coordinates.  Can translate back to gInterval if correspondance is setup.
 * 
 * Data intervals must have rapid data access for EM algorithm.
 * @version 0.1
 * @date 2022-04-28
 * 
 */
#include "Data.h"

#include <string>
#include <vector>
#include <iostream>

#include "split.h"
#include "Intervals.h"


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
 * @param v_identifier  A name/identifier for this interval
 */
dInterval::dInterval(std::string v_identifier) {
  ID = v_identifier;
  minX = 0;
  maxX = 0;
  X = NULL;
  XN = 0;
  SCALE = 1;
  N = 0;

  belongsTo = NULL;
}

// empty constructor
dInterval::dInterval() {
  ID = "empty";
  minX = 0;
  maxX = 0;
  X = NULL;
  XN = 0;
  SCALE = 1;
  N = 0;

  belongsTo = NULL;
}

/**
 * @brief This is a standard content dump function.  Primarily
 * used in debugging.  Notice the string does NOT end in a newline.
 * 
 * @return std::string 
 */
std::string dInterval::write_out() {
  std::string text = ("#" + ID + ":min:" + std::to_string((int)minX) + ":max:" 
		+ std::to_string((int)maxX) + ":num_bins:" + std::to_string((int)XN) + ":total:"
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
 * Note that this is not necessarily a genomic position.
 * 
 * @return double 
 */
double dInterval::position(int x) {
  return X[0][x];
}
  
double dInterval::sum_Region() {
  return N;
}

