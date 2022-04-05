/**
 * @file Intervals.h
 * @author Robin Dowell
 * @brief Header file for data interval class (dInterval)
 *
 * Two kinds of intervals:  genomic (gIntervals) and data (dIntervals)
 * gIntervals maintain genomic coordinates and are the fundamental datatype
 * within the Interval trees (see ITree.cpp).  These are all obtained from 
 * bed files (bed3, bed4, bed6, bed9, bed12).  
 * 
 * dIntervals contain data (two strands) per an interval but do so in zero
 * based coordinates.  Can translate back to gInterval if correspondance is setup.
 * 
 * @version 0.1
 */
#ifndef Intervals_H 
#define Intervals_H 

#include <string>
#include <vector>

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

  // Constructors
  gInterval(std::string, double, double, std::string);  // input is BED4
  gInterval();

  /* FUCTIONS: */
  // Reporting out 
  std::string write_out();
  std::string write_asBED();

  void setfromBedLine(std::string);  // converts from a single line from file 

  bool Overlap(gInterval *);
  bool Contains(double point);

protected:
  void setBED4fromStrings(std::vector<std::string> lineArray); // helper function

private:
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
  bed6(std::string, double, double, std::string, int, std::string);
  bed6();

  /* FUCTIONS: */
  // Reporting out 
  std::string write_out();
  std::string write_asBED();

  void setfromBedLine(std::string);  // converts from a single line from file 
};

/**
 * @brief dInterval 
 * 
 * Data interval class which contains two strands of data
 * over an interval.  Data is expected to be forward and reverse
 * strand info (as double). 
 * 
 * The data may, or may not, be scaled and binned relative 
 * to genomic coordinates.
 * 
 * Goals: Storage of the data should be hidden from the user.
 * 
 * @author Robin Dowell
 */
class dInterval {

public:
  std::string ID; //!< Does each data interval need a unique name?

  // Ultimately unsure if we need to keep these.
  double minX; //!< This is the minimum value of the interval
  double maxX;   //!< This is the maximum value of the interval 

  // Constructors
  dInterval(std::string, int , int);
  dInterval();

  /* FUCTIONS: */
  // Reporting out 
  std::string write_out();

  /*
   * Wrappers to effectively obtain an iterator across data.
   * This is a bit clunky currently, but a step in the right direction.
   * In all of these cases, you are expected to iterate from 0 to num_elements()
   */
  double num_elements();
  double forward(int);
  double reverse(int);
  double position(int);

  double sum_Region();  // both strands
  double sum_forward(); // forward only
  double sum_reverse(); // reverse only

  void load_forward_strand();
  void load_reverse_strand();
  void scale_and_bin();

private:
	/**
	 * @brief This (X) is the internal representation of the data.
	 * Vector[0] is coordinate (possibly scaled); [1] is forward (summed for bin)
	 * [2] is reverse (summed for bin).  
	 * 
	 * This is the meat and potatoes data representation that is used in the EM (fit2).  
	 * 
	 * These data elements are currently identical to how originally written 
	 * in segment class.
	 */
	double ** X;  //!< Smoothed data inner is [3] dimensions
	double XN; //!< total number of bins
	double SCALE;  //!< scaling factor

	// I think these are for convenience (calculate once)
	double N;	//!< Total sum of values 
	double fN;	//!< Sum of forward values 
	double rN;	//!< Sum of reverse values 
};

#endif
