/**
 * @file gInterval.h
 * @author Robin Dowell
 * @brief Header file for genomic interval container class -- has interval (bed) 
 * and data (dInterval)
 *
 * @version 0.1
 */
#ifndef gInterval_H 
#define gInterval_H 

#include <string>
#include <vector>

#include "Bed.h"
#include "Data.h"

/**
 * @brief gInterval 
 * 
 * Container that links a basic interval (bed) to the data therein (dInterval).
 * 
 * bed format is 0-based, half open coordinate scheme.
 * 
 */
class gInterval {
public:
  bed12  *interval;
  dInterval *tdf;

  // Constructors
  gInterval();  // default

  /* FUCTIONS: */
  // Reporting out 
  std::string write_out();
};

#endif
