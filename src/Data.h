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
 * A gInterval container links a bed to a dInterval.  
 * 
 * @version 0.1
 */
#ifndef Data_H 
#define Data_H 

#include <string>
#include <vector>

class bed12;   // Forward declaration
class gInterval;

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

  gInterval *parent;

  //Constructors
  RawData();
  RawData(gInterval *);
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
 * Data interval class which contains two strands of data over an interval. 
 * 
 * The data may, or may not, be scaled and binned relative 
 * to genomic coordinates.
 * 
 * X[0] is data coordinate
 * X[1] is forward coverage
 * X[2] is reverse coverage
 * 
 * The coordinate here is zero based X[0][0] = 0, with step size delta
 * (e.g. X[0][i] = i*delta), where i is the current index.
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

  gInterval *seg;

  // Constructors
  dInterval();
  dInterval(int, int);
  // dInterval(RawData *, int, int); // Convert RawData into binned/scaled data 
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
  double getLength();   // Note data matrix is 0 ... getLength

  // Summary information on this data
  double sumForward();
  double sumReverse();
  double sumInterval(int, int, char);  // on single strand
  double sumAlldata();  // both strands

  // Convert data position to index and index to data position
  int getIndexfromData(double);  // given data coordinate, get closest index
  double getDataCoordfromIndex(int);  // given index, what is data coordinate
  double getMinGenomeCoord();
  double getMaxGenomeCoord();

  

  // Functions for doing the conditioning
  void setBinNum(double length);
  void initializeData(int startCoord);  // Sets up the internal matrix
  void BinStrands(RawData *data);
  void ScaleDown(int);    // Convert to zero based coords
  void CompressZeros();   // Remove positions that are zero on both strands

  // Other useful functions
  // int getWithinRangeofPosition(double position, double dist);
  void DeallocateX();   // Deallocates X, leaves other variables intact.

};

#endif
