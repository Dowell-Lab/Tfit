/**
 * @file Segment.h
 * @author Robin Dowell
 * @brief Header file for segment class (gInterval + dInterval)
 *
 * @version 0.1
 */
#ifndef Segment_H 
#define Segment_H 

#include <string>
#include <vector>

#include "Intervals.h"

/**
 * @brief Segment 
 *
 * A Segment of the genome contains both a gInterval (i.e. genomic coordinates)
 * and a dInterval (data over the region).
 * 
 */
class Segment {

public:
  // Should these be pointers?
  bed6 coords;
  dInterval data;

  // Constructors
  Segment();

  /* FUCTIONS: */
  // Reporting out 
  std::string write_out();
  std::string write_Interval();
  std::string write_data();

protected:
private:
};


#endif