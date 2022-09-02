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
class bed12;   // Forward declaration
class dInterval;

/**
 * @brief Coverage per genomic coordinate
 * [coordinate, coverage]
 * Forces coverage to be positive (i.e. always takes abs()).
 */
class PointCov {
  public:
  double coordinate;
  double coverage;

  //Constructor
  PointCov();
  PointCov(double v_coord, double v_cov);

  std::string write_out();

  static bool sortOnCoordComp(const PointCov pt1, const PointCov pt2);
  static bool sortOnCovComp(const PointCov pt1, const PointCov pt2);

};

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
	std::vector<PointCov> forward;
	std::vector<PointCov> reverse; 

  double minX; //!< This is the minimum value of the interval
  double maxX;   //!< This is the maximum value of the interval 

  bed12 *belongsTo;  // The coordinates associated with this data
  dInterval *cdata;   // The conditioned data.

  //Constructors
  RawData();
  RawData(bed12 *);
  //Destructor
  ~RawData();

  std::string write_out();  // Debugging

  double Length();    // length of region with data (maxX-minX)   
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

  RawData *raw;  // The raw data from which this was built.

  // Constructors
  dInterval();
  dInterval(RawData *, int, int); // Convert RawData into binned/scaled data 
  ~dInterval();   // Destructor

  /* FUCTIONS: */
  // Reporting out 
  std::string write_out();  // debugging
  std::string data_dump();  // debugging

  // Accessors, makes code far more readable and can be used for iteration.
  double num_elements();  // equivalent to bins
  double forward(int);    // equivalent to X[0][i]
  double reverse(int);    // equivalent to X[1][i]
  double position(int);   // equivalent to X[2][i]
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
  // int getWithinRangeofPosition(double position, double dist);
  void DeallocateX();   // Deallocates X, leaves other variables intact.

	};

#endif
