/**
 * @file Intervals.h
 * @author Robin Dowell
 * @brief Header file for genomic interval class (gInterval and bed6)
 *
 * gIntervals maintain genomic coordinates and are the fundamental datatype
 * within the Interval trees (see ITree.cpp).  These are all obtained from 
 * bed files (bed3, bed4, bed6, bed9, bed12).  
 * 
 * They have a pointer to data when the interval contains data.  If the interval
 * does not contain data, this pointer is NULL.
 * 
 * @version 0.1
 */
#ifndef Intervals_H 
#define Intervals_H 

#include <string>
#include <vector>

class RawData;    // Forward declaration
class Seeds;    // Forward delcaration

/**
 * @brief gInterval
 * 
 * Basic genomic interval.  Equivalent to basic information stored in bed4 file.
 * This is also the fundamental unit maintained in interval trees.
 * 
 * BED4 is 0-based, half open coordinate scheme.
 * 
 */
class gInterval {

public:
  std::string identifier; // field 4 in BED
  std::string chromosome;  // field 1 in BED
  double start, stop;   // genomic coordinates, fields 2 and 3 in bed

  RawData *data;  // pointer to data if this interval has it.

  // Constructors
  gInterval();  // default
  gInterval(std::string, double, double, std::string);  // input is BED4

  /* FUCTIONS: */
  // Reporting out 
  std::string write_out();

  // Read and write a bed* line
  std::string write_asBEDline();
  void setfromBedLine(std::string);  // converts from a single line from file 
  void setBEDfromStrings(std::vector<std::string> lineArray); // helper for setfromBedLine

  // Comparison of points/intervals
  bool Overlap(gInterval *);    // Does your interval of interest overlap this one?
  bool Contains(double point);  // Is a coordinate in this region?

};

/**
 * @brief a gInterval that also includes score and strand information
 * 
 * Follows the BED6 format.
 *
 * BED6 is 0-based, half-open coordinate scheme. 
 * In BED6 score is between 0 and 1000 and used in heatmap generation.
 * In BED6 strand is "." (no strand), "+" (positive) or "-" (negative).
 */
class bed6: public gInterval {
public:
  char strand;
  int score;

  // Constructors
  bed6();
  bed6(std::string, double, double, std::string, int, std::string);

  // Reporting out 
  std::string write_out();

  // Read and write a bed* line
  std::string write_asBEDline();
  void setfromBedLine(std::string);  // converts from a single line from file 
  void setBEDfromStrings(std::vector<std::string> lineArray); // helper function
};

/**
 * @brief By default, everything is a bed12 because this allows
 * for seeds to the EM.   Therefore, this is where data can be added.
 */
class bed12: public bed6 {
  public: 
  Seeds *seeds;  

  // Constructors
  bed12();
  bed12(std::string, double, double, std::string, int, std::string, std::string);
  bed12(std::string, double, double, std::string);
  ~bed12();

  // Functions
  std::string write_out();

  // Read and write a bed* line
  std::string write_asBEDline(); // std::string writeToBed12 ();
  void setfromBedLine(std::string);  // converts from a single line from file 
  void setfromLastSix(std::string); // converts ONLY the last 6 fields of bedfile

  // add points to data (RawData):
  void addDataPoint(double start,double stop,double cov,bool expand);

};

#endif
