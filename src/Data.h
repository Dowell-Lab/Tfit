/**
 * @file Data.h
 * @author Robin Dowell
 * @brief Header file for data interval classes (RawData, dInterval)
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

class gInterval;   // Forward declaration
class dInterval;

/**
 * @brief This is a container for the raw data as read from a bedgraph.  
 * This class is responsible for data sanity checks.  
 * 
 */
class RawData {
  public:
    // These are coordinate coverage vectors that are used to load in bedgraph
    // We read the bedgraph contents in directly, without checking for repeated entries.
	std::vector< std::vector<double> > forward; 
	std::vector< std::vector<double> > reverse; 

  double minX; //!< This is the minimum value of the interval
  double maxX;   //!< This is the maximum value of the interval 

  gInterval *belongsTo;  // The coordinates associated with this data
  dInterval *cdata;   // The conditioned data.

  //Constructors
  RawData();
  RawData(gInterval *);

  std::string write_out();  // Debugging

  double Length();
  void ClearData ();    // Deallocates the forward and reverse vectors;
  void addDataPoints(double st, double sp, double cov);

  std::string data_dump();

  void Sort();
  void RemoveDuplicates();
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
	double ** X;  //!< Smoothed data inner is [3] dimensions
	int bins;  //!< total number of bins
  int delta;  // step size (nts per bin)
  int scale;  // all positions are divided by this factor e.g. 1 -> 1/scale

	double N;	//!< Total sum of values 

  RawData *raw;  // The raw data from which this was built.

  // Constructors
  dInterval();
  dInterval(RawData *, int, int); // Convert RawData into binned/scaled data 

  /* FUCTIONS: */
  // Reporting out 
  std::string write_out();  // debugging

  // These can be used as an iterator:
  double num_elements();
  double forward(int);
  double reverse(int);
  double position(int);

  double sum_Region();  // both strands

  // Functions for doing the conditioning.
  void initializeData(int length);  // Sets up the internal matrix
  void BinStrands(RawData *data);
  void ScaleDown(int);    // Convert to zero based coords
  void CompressZeros();   // Remove positions that are zero on both strands

  void ClearX();   // Deallocates X, leaves other variables intact.

  std::string data_dump();

  // Need functions for converting between three coordinate systems: genomic, conditioned, indexed.

	};

#endif
