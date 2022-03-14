/**
 * @file Intervals.h
 * @author Robin Dowell
 * @brief Header file for data interval class (dInterval)
 *
 * @version 0.1
 */
#ifndef Intervals_H 
#define Intervals_H 

#include <string>
#include <vector>

/**
 * @brief Data interval class which contains two strands of data
 * over an interval.  
 * 
 * Data Intervals are expected to represent two strands of data 
 * over a set of positions.  The data may, or may not, be scaled
 * and binned relative to genomic coordinates.
 * 
 * Questions: do we need to keep the genomic coordinates here?
 * Does each interval need a unique identifier?
 * 
 * Goals: Storage of the data should be hidden from the user.
 * But for efficiency, the iterator should be indexable?
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
