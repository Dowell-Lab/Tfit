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
 * Responsibilities: 
 *  -- convert between genomic and dInterval coordinates
 *  -- given RawData, make a dInterval
 * 
 */
class gInterval {
public:
  bed12  *interval;
  dInterval *tdf;
  RawData *raw;

  // Constructors
  gInterval();  // default

  /* FUCTIONS: */
  // Reporting out 
  std::string write_out();

  void createDataInterval(int v_delta, int v_scale);
  void raw2dInterval();

  //Convert between index, data and genomic coordinates 
  int getIndexfromGenomic(double);   // given genomic coordinate, give index 
  double getGenomeCoordfromIndex(int);  // given an index to the dInterval, what is the genomic Coord
  double getGenomefromData(double); // given data coordinate, get genomic Coord
  double getDataCoordfromGenomeCoord(double); /// given genomic coords, what is data coords?
};

#endif
