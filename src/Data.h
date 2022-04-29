/**
 * @file Data.h
 * @author Robin Dowell
 * @brief Header file for data interval classes (RawData, dInterval, Segment)
 * 
 * RawData holds the raw, unnormalized, unconditioned data.
 *
 * dIntervals contain data (two strands) per an interval but do so in zero
 * based coordinates.  
 * 
 * A segment links a gInterval to RawData and dIntervals.  It can also translate 
 * coordinates among these representations. 
 * 
 * @version 0.1
 */
#ifndef Data_H 
#define Data_H 

#include <string>
#include <vector>

class bed6;   // Forward declaration

/**
 * @brief This is a container for the raw data as read from a bedgraph.  
 * This class is responsible for data sanity checks.  
 * 
 */
class RawData {
  public:
    // These are coordinate coverage vectors that are used temporarily to load in bedgraph
	std::vector< std::vector<double> > forward; 
	std::vector< std::vector<double> > reverse; 

  //Constructors
  RawData();
};

/**
 * @brief dInterval 
 * 
 * Data interval class which contains two strands of data
 * over an interval.  Data is expected to be forward and reverse
 * strand info.
 * 
 * The data is always in zero based coordinates.  
 * 
 * The data may, or may not, be scaled and binned relative 
 * to genomic coordinates.
 * 
 * @author Robin Dowell
 */
class dInterval {
public:
  std::string ID; //!< Does each data interval need a unique name?

  // Ultimately unsure if we need to keep these.
  double minX; //!< This is the minimum value of the interval
  double maxX;   //!< This is the maximum value of the interval 

  bed6 *belongsTo;  // The coordinates associated with this data

  // Constructors
  dInterval(std::string);
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

private:    // Should these be private?
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
	double N;	//!< Total sum of values 
};

#endif
