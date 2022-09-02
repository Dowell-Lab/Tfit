/**
 * @file Data.cpp
 * @author Robin Dowell 
 * @brief Class code for data intervals
 * 
 * Current design: 
 * 
 * RawData keeps points [coord,coverage] for both strands.  It's a container for
 * reading in from bedGraphs.  Once the bedGraph is fully read, then we can
 * convert the RawData into a dInterval.
 * 
 * A dInterval contains points [coord,forward_coverage,reverse_coverage] that 
 * is zero based (e.g. RawData->minX becomes Zero), binned (delta) and scaled.
 * 
 * Data intervals must have rapid data access for EM algorithm.
 * @version 0.1
 * @date 2022-04-28
 * 
 */
#include "Data.h"

#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "split.h"
#include "helper.h"
#include "Intervals.h"

/**********  PointCov ****************/

PointCov::PointCov() {
  coordinate = 0;
  coverage = 0;
}

PointCov::PointCov(double v_coord, double v_cov) {
  coordinate = v_coord;
  coverage = abs(v_cov);
}

std::string PointCov::write_out() {
  std::string output;
  output = "[" + tfit::prettyDecimal(coordinate,0) + "," 
    + tfit::prettyDecimal(coverage, 2) + "]";
  return output;
}

bool PointCov::sortOnCoordComp(const PointCov pt1, const PointCov pt2) {
   return (pt1.coordinate < pt2.coordinate);
}

bool PointCov::sortOnCovComp(const PointCov pt1, const PointCov pt2) {
   return (pt1.coverage < pt2.coverage);
}

/**********************  RawData ********************/
RawData::RawData() {
  minX = maxX = 0;
  belongsTo = NULL; 
  cdata = NULL;
}

RawData::RawData(bed12 *v_ginterval) {
  minX = maxX = 0;
  belongsTo = v_ginterval; 
  cdata = NULL;
}


/**
 * @brief Length of the RawData (max-min)
 * Assumes half open coordinates.
 * 
 * @return double returns length of this segment
 */
double RawData::Length() {
  return (maxX-minX);   // Half open coordinates means no +1
}

/**
 * @brief Release the memory of raw data
 * This will be useful/necesary for min memory version.
 */
void RawData::freeDataMemory () {
  forward.clear();
  reverse.clear();
}

/**
 * @brief Adds RawData.  Assumes the joint bedGraph conventions, 
 * namely zero based, half open coordinates and the sign of the 
 * coverage indicates strand.  A range of points is input as 
 * each individual point within the set.
 * 
 * @bug Should this check to see if the point is in the range of
 * belongsTo?
 * 
 * @param st  start (this point is included)
 * @param sp   stop (with half open, this is strickly less than!)
 * @param cov  signed count/coverage at this point
 */
void RawData::addDataPoints(double start, double stop, double cov) {
  if (maxX == 0) { minX = start; maxX = stop; } // first data point

  // Adjust min/max as gather data points
  if (start < minX) minX = start;
  if (stop > maxX) maxX = stop;

  // Now add per position coverage info, 
  // Half open coordinates make this strickly less than!
  for (int i = start; i < stop; i++) {
    PointCov newpoint((double)i, cov); // PointCov will force abs(cov)!
    if (cov >= 0) {
      forward.push_back(newpoint);
    } else { 
      reverse.push_back(newpoint);
    }
  }
}

/**
 * @brief Once all the data is read in, we can sort the
 * forward and reverse strands.
 */
void RawData::Sort() {
  // Sort vector<double> by first entity (e.g. coordinates)
  std::sort(forward.begin(),forward.end(), PointCov::sortOnCoordComp);
  // std::cout << "After sort: " + data_dump() << std::endl;
  std::sort(reverse.begin(),reverse.end(), PointCov::sortOnCoordComp);
}

/**
 * @brief Remove duplicate entries for any position.
 * Arbitrarily keeps the first one of duplicates
 * Erase is notoriously slow, so this could be costly for big datasets.
 */
void RawData::RemoveDuplicates() {
  Sort();
  for (int i = 1; i < forward.size(); i++) {
    if (forward[i-1].coordinate == forward[i].coordinate)  {
      // This is a duplicate!  Should we throw an error!?!?
      forward.erase(forward.begin()+i);
      // Because this shifts the index, need to compare the
      // new i so we'll shift the index.
      i--;
    }
  }  
  for (int i = 1; i < reverse.size(); i++) {
    if (reverse[i-1].coordinate == reverse[i].coordinate)  {
      // This is a duplicate!  Should we throw an error!?!?
      reverse.erase(reverse.begin()+i);
      i--;
    }
  }  
}

/**
 * @brief Stringify the object for debugging purposes
 * 
 * @return std::string  contents of object as a string.
 */
std::string RawData::write_out() {
  std::string ID;
  if (belongsTo != NULL ) { ID = belongsTo->identifier; }
  else { ID = "noID"; }

  std::string output = ID + ":";
  output += tfit::prettyDecimal(minX,2) + ":" + tfit::prettyDecimal(maxX,2);

  output += "  " + tfit::prettyDecimal(forward.size(), 2);
  output += "  " + tfit::prettyDecimal(reverse.size(), 2);
  return output;
}

/**
 * @brief Outputs full contents of the points, this can be large
 * use with caution!
 * 
 * @return std::string  The complete contents of data as string (will be big!).
 */
std::string RawData::data_dump() {
  std::string output = "Forward: ";
  for (int i=0; i < forward.size(); i++) {
    output += forward[i].write_out();
  }
  output += "\nReverse: ";
  for (int i=0; i < reverse.size(); i++) {
    output += reverse[i].write_out();
  }
  return output;
}

/**
 * @brief Cleanup function.  Necessary to avoid memory leaks.
 */
RawData::~RawData() {
  freeDataMemory();
  // What about belongsTo and cdata?
  if (cdata != NULL) {
    delete(cdata);
  }
}

/****************** dInterval *********************/

/**
 * @brief Constructors: dInterval class
 * @author Robin Dowell
 *
 * Purpose: create/allocate instances of a data Interval 
 */
dInterval::dInterval() {
  X = NULL;
  delta = 1;
  scale = 1;
  bins = -1;  // illegal value since we dont know this yet.
  raw = NULL;
}

/**
 * @brief Construct a new dInterval::dInterval object
 * 
 * @param data      The rawdata
 * @param v_delta   The binning size (nts per bin)
 * @param v_scale   The scaling factor
 */
dInterval::dInterval(RawData *data, int v_delta, int v_scale) {
  raw = data;
  raw->RemoveDuplicates();      // Is this necessary?  It's potentially time consuming!
  delta = v_delta;  scale = v_scale; 
  if (scale <= 0) { scale = 1;}

  // Establish num of bins
  bins =  raw->Length()/delta;
  if (bins < 1) { // Min value
    bins = 1; 
  } else {
    if (delta*bins < raw->Length()) {
      bins += 1;    // The end stuff which isn't of length delta!
    }
  } 

  initializeData(raw->minX);
  BinStrands(raw);
  ScaleDown(raw->minX);
  // std::cout << data_dump() << std::endl;
  CompressZeros();
  //std::cout << "AFTER: " + to_string(bins) << std::endl;
}

/**
 * @brief Cleanup function to avoid memory leaks!
 */
dInterval::~dInterval() {
  if (X != NULL) { DeallocateX();}
}

/**
 * @brief Stringify object for debugging!
 * 
 * @return std::string contents of object as a string.
 */
std::string dInterval::write_out() {
  std::string output;
  if (raw != NULL) {
    output = raw->write_out();
  }
  output += "\tbins: " + tfit::prettyDecimal(bins,0) + "," 
      + tfit::prettyDecimal(delta,4) + "," + tfit::prettyDecimal(scale,3);
  return output;
}

/**
 * @brief Outputs full contents of the points, this can be large
 * use with caution!
 * 
 * @return std::string Contents of X[][] array as string.
 */
std::string dInterval::data_dump() {
  std::string output = "Points: ";
  for (int i = 0; i < bins; i ++ ){
    output += "[" + tfit::prettyDecimal(position(i),2) + ":" + 
          tfit::prettyDecimal(forward(i),2) + "," + tfit::prettyDecimal(reverse(i),2) + "]";
  }
  return output;
  
}

/**
 * @brief The number of data points in the interval.  
 * This could be nucleotides (smallest unit) to bins (groups of nts)
 * to the number of (possibly random) points along the interval.
 * 
 * @return double  Size of the interval in usable steps.
 */
double dInterval::num_elements() {
  return bins;
}

/**
 * @brief Data from forward strand at xth index.
 * Note a "unit" here could be nucleotides, bins (groups of nts),
 * or points along the interval.
 * 
 * @arg x The index of an element of the forward strand data.
 * Note that this is not a position (genomic or scaled).
 * 
 * @return double 
 */
double dInterval::forward(int x) {
  return X[1][x];
}

/**
 * @brief Data from reverse strand at xth index.
 * Note a "unit" here could be nucleotides, bins (groups of nts),
 * or points along the interval.
 * 
 * @arg x The index an element of the reverse strand data.
 * Note that this is not a position (genomic or scaled).
 * 
 * @return double 
 */
double dInterval::reverse(int x) {
  return X[2][x];
}

/**
 * @brief The position (in either genomic coords or data coords) at xth index.
 * 
 * @arg x The index an element 
 * 
 * @return double 
 */
double dInterval::position(int x) {
  return X[0][x];
}

/**
 * @brief Return the length of this dInterval.
 * Will use the raw object's length if available (genome coords?),
 * otherwise assumes the last element of X[][] is the length (i.e. zero based).
 * 
 * @return double  Length of this region/segment.
 */
double dInterval::getLength() {
  if (raw == NULL) {
    if (X == NULL) {
      return 0;     // No data
    }
    return position(bins-1);    // conditioned data, no raw
  }
  return (raw->Length());   //Genome coords length
}

/**
 * @brief Helper function: sum the data on forward strand
 * 
 * @return double SUM (forward strand coverage)
 */
double dInterval::sumForward() {
  double sum = 0;
  for(int i = 0; i < num_elements(); i++) {
     sum += forward(i);
  }
  return sum;
}

/**
 * @brief Helper function: sum the data on reverse strand
 * 
 * @return double SUM (reverse strand coverage)
 */
double dInterval::sumReverse() {
  double sum = 0;
  for(int i = 0; i < num_elements(); i++) {
     sum += reverse(i);
  }
  return sum;
}

/**
 * @brief Total sum of all the data. 
 * Currently a stored value, not computed.
 * 
 * @return double Sum of all reads on both strands in interval.
 */
double dInterval::sumAlldata() {
  double sum = 0;
  sum += sumForward();
  sum += sumReverse();
  return sum;
}

/**
 * @brief Sum the data on one strand between indicies.
 * Assumes start and stop are indexes into dInterval (e.g. between 0 and X[0][bins-1]).
 * 
 * @return double Sum of data on one strand between indicies.
 */
double dInterval::sumInterval(int start, int stop, char strand) {
  int strandidx = 1;
  if (strand == '-') { strandidx = 2;}

  double sum = 0;
  for (int i= start; i< stop; i++) {
    sum += X[strandidx][i];
  }
  return sum;
}

/*****Convert between coordinate systems *****/

/**
 * @brief Convert genomic coordinate into the correct index
 * 
 * @param genomicCoord   genomic coordinate 
 * @return int    correct index into X[]
 */
int dInterval::getIndexfromGenomic(double genomicCoord) {
  return int((genomicCoord - raw->minX)/scale);
}

/**
 * @brief Get the Index from a data coordinate 
 * Uses a binary search approach.  Returns the index with 
 * the closest data coordinate.
 * 
 * @param dataCoord  Scaled/binned data coordinate
 * @return * int  correct index into X[]
 */
int dInterval::getIndexfromData(double dataCoord) {
  int indexMin = 0;
  int indexMax = bins-1;
  int indexCenter;
  bool found = false;

  while (!found) {
    if ((indexMax - indexMin) < 2) {
      found = true;
      // Take closest
      int distmin = abs(dataCoord - position(indexMin));
      int distmax = abs(dataCoord - position(indexMax));
      if(distmin < distmax) {
        return indexMin;         
      } else {
        return indexMax;         
      }
    } else if (dataCoord == position(indexMin)) {
      found = true;
      return indexMin;
    } else if (dataCoord == position(indexMax)) {
      found = true;
      return indexMax;
    }
    indexCenter = (int)(indexMax -indexMin)/2;

    if ((dataCoord > position(indexMin)) && (dataCoord < position(indexCenter))) {
      indexMax = indexCenter;
    } else { 
      // Edge case
      if (indexCenter == indexMin) { indexMin += 1;}
      else {indexMin = indexCenter;}
    }
    if (indexMax < indexMin) { found = true;}
  }
  std::cerr << "Error!! Unreachable in getIndexfromData." << std::endl;
  return -1;
}

/**
 * @brief Convert index to genomic coordinate
 * Note that because of binning these will always be on 
 * minX + delta intervals.
 * 
 * @param index  position index into X[]
 * @return double   genomic Coordinate cooresponding
 */
double dInterval::getGenomeCoordfromIndex(int index) {
  double coord = ((index * scale) + raw->minX);
  // Necessary to account for the last bin being uneven sized.
  if (coord > raw->maxX) { return raw->maxX; }
  else { return coord; }
}

/**
 * @brief Convert data coordinates to genomic coordinates
 * 
 * @param dataCoord  Scaled/binned data coordinate
 * @return double genomic coordinates corresponding
 */
double dInterval::getGenomefromData(double dataCoord) {
   int index = getIndexfromData(dataCoord);
   return getGenomeCoordfromIndex(index);
}

/**
 * @brief given an index, what is the data coordinate?
 * 
 * @param index  position index into X[]
 * @return double Scaled/binned data coordinate
 */
double dInterval::getDataCoordfromIndex(int index) {
  return X[0][index];
}

/**
 * @brief Given Genomic coordinates, what are data coordinates?
 * 
 * @param genomicCoord   genomic coordinate 
 * @return double Scaled/binned data coordinate
 */
double dInterval::getDataCoordfromGenomeCoord(double Gcoord) {
  int index = getIndexfromGenomic(Gcoord);
  return getDataCoordfromIndex(index);
}

/***** Conditioning / Scaling / Binning Data *****/
/**
 * @brief Allocate and setup the X matrix given the size of the RawData
 * 
 * @param startCoord the left edge of this interval
 */
void dInterval::initializeData(int startCoord) {
  // if (X != NULL) {
    // std::cout << "Clearing old data!" << std::endl;
    // ClearX(); 
  // }
  X = new double*[3];
  for (int j = 0 ; j < 3;j++){
    X[j] 		= new double[bins];
  }
  X[0][0] 		= double(startCoord);
  // std::cout << "MIN: " + to_string(minX) << std::endl;
  X[1][0]=0,X[2][0]=0;
	
  for (int i = 1; i < bins; i++){
    X[0][i] 	= X[0][i-1] + delta;
    X[1][i] 	= 0;
    X[2][i] 	= 0;
  }
}

/**
 * @brief   Converts rawData into binned data. 
 * We have three arrays to walk through: 
 *   X[][binnum] : conditioned data  i = 0 to bins-1
 *   forward[fi][2]  : raw data for forward strand  fi = 0; fi < forward.size()
 *   reverse[ri][2]  : raw data for reverse strand  ri = 0; ri < reverse.size()
 * 
 * @param data  the RawData container for the strand coverage info
 */
void dInterval::BinStrands(RawData *data) {
  int fi = 0; int ri = 0;   // forward and reverse indicies
  double binStart, binEnd;    // genomic coordinate bounds for a given bin
  for (int i = 0; i < bins; i++) {  // Each bin in X
    // Genomic coordinate range for this bin:
    binStart = position(i); 
    binEnd = binStart + delta;
    while ((fi < data->forward.size()) &&     // Still more forward data
           (data->forward[fi].coordinate >= binStart) &&  // This coordinate is above start
           (data->forward[fi].coordinate < binEnd)) {   // And below end of range for this bin
      X[1][i] += data->forward[fi].coverage;    // Add to the bin
      fi++;       // Increment the forward index
    }
    // As above, but now for the reverse strand:
    while ((ri < data->reverse.size()) && 
          (data->reverse[ri].coordinate >= binStart) && 
          (data->reverse[ri].coordinate < binEnd)) {
      X[2][i] += data->reverse[ri].coverage;
      ri++;
    }
  }
}

/**
 * @brief Scales down the coordinates to Zero based.
 * 
 * @param startCoord the left edge of this region (genomic coords) 
 */
void dInterval::ScaleDown(int startCoord) {
  for (int i = 0; i < bins; i ++ ){
    X[0][i] 	= (position(i)-startCoord)/scale;
  }
}

/**
 * @brief Get rid of data points where both strands are zero.
 * 
 */
void dInterval::CompressZeros() {
  int realN 		= 0;	// number of non-zero bins
  for (int i = 0; i < bins;i++){
    if ((forward(i) > 0) or (reverse(i) > 0)) {
      realN++;
    }
  }
  // going to remove the zero bins
  double **newX = new double *[3];
  for (int j = 0; j < 3; j++)
  {
    newX[j] = new double[realN];
  }
  int j = 0;
  for (int i = 0; i < bins; i++) {
    if ((forward(i) > 0) or (reverse(i) > 0)) {
      newX[0][j] = position(i);
      newX[1][j] = forward(i);
      newX[2][j] = reverse(i);
      j++;
    }
  }
  if (realN != j) {
    printf("WHAT? %d,%d\n", j, realN);
  }
  this->DeallocateX();
  X = newX;
  bins = realN;
}

/**
 * @brief  This is currently a rewrite of Joey's get_nearest_position
 * into this data structure.
 * 
 * @param position  Typically a mu, in transformed coordinates
 *  -- should this instead be in genomic coordinates?  Or is that an
 *  alternative function?
 * @param dist    Typically s(1/lambda)
 * @return int    index of the position closest to position + dist (signed dist)
int dInterval::getWithinRangeofPosition(double qspot, double dist) {
	int i;

	if (dist < 0 ){   // negative strand
		i=0;
		while (i < (bins-1) and (position(i)-qspot) < dist){
			i++;
		}
	}else{
		i=bins-1;
		while (i >0 and (position(i)-qspot) > dist){
			i--;
		}
	}
	return i;
}
 */

/**
 * @brief Cleans up existing X matrix (memory clearing)
 * 
 */
void dInterval::DeallocateX() {
  if (X == NULL) return;
  for (int i = 0; i < 3; i++) {
    delete [] X[i];
  }
  delete [] X;
  X = NULL;
}
