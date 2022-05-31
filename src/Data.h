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
 * This class is responsible for data sanity checks.  The raw data
 * can be deleted (freeing memory). 
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
  //Destructor
  ~RawData();

  std::string write_out();  // Debugging

  double Length();    
  void freeDataMemory();    //!< Deallocates the forward and reverse vectors;
  void addDataPoints(double st, double sp, double cov);

  std::string data_dump();  //!< a debugging function

  void Sort();    //!< Sorts both forward and reverse by coordinate
  void RemoveDuplicates();  //!< Removes any coordinate with zero reads on both strands
};

/**
 * @brief dInterval 
 * 
 * Data interval class which contains two strands of data
 * over an interval. 
 *
 * The data is always in zero based coordinates.  
 * 
 * The data may, or may not, be scaled and binned relative 
 * to genomic coordinates.
 * 
 * This is the format of data expected by the EM algorithm.
 * 
 * @author Robin Dowell
 */
class dInterval {
public:
	double ** X;  //!< Smoothed data inner is [3] dimensions
	int bins;  //!< total number of bins (e.g. size of X[])
  int delta;  //!< step size (nts per bin)
  int scale;  //!< all positions are divided by this factor e.g. 1 -> 1/scale

	double N;	//!< Total sum of values (convenience variable)

  RawData *raw;  // The raw data from which this was built.

  // Constructors
  dInterval();
  dInterval(RawData *, int, int); // Convert RawData into binned/scaled data 
  ~dInterval();   // Destructor

  /* FUCTIONS: */
  // Reporting out 
  std::string write_out();  // debugging
  std::string data_dump();  // debugging

  // Accessors, can be used to iterate through:
  double num_elements();
  double forward(int);
  double reverse(int);
  double position(int);
  double getLength();

  // Summary information on this data
  double sumForward();
  double sumReverse();
  double sumInterval(int, int, char);  // on single strand
  double sumAlldata();  // both strands

  //Convert between index, data and genomic coordinates 
  int getIndexfromGenomic(double);   // given genomic coordinate, give index 
  int getIndexfromData(double);  // given data coordinate, get closest index
  double getGenomeCoordfromIndex(int);  // given an index to the dInterval, what is the genomic Coord
  double getGenomefromData(double);
  double getDataCoordfromIndex(int);
  double getDataCoordfromGenomeCoord(double);

  // Functions for doing the conditioning
  void initializeData(int length);  // Sets up the internal matrix
  void BinStrands(RawData *data);
  void ScaleDown(int);    // Convert to zero based coords
  void CompressZeros();   // Remove positions that are zero on both strands

  // Other useful functions
  int getWithinRangeofPosition(double position, double dist);
  void DeallocateX();   // Deallocates X, leaves other variables intact.

	};

#endif
