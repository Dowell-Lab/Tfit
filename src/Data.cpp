/**
 * @file Data.cpp
 * @author Robin Dowell 
 * @brief Class code for data intervals
 * 
 * Current design: 
 * dIntervals contain data (two strands) per an interval but do so in zero
 * based coordinates.  Can translate back to gInterval if correspondance is setup.
 * 
 * Data intervals must have rapid data access for EM algorithm.
 * @version 0.1
 * @date 2022-04-28
 * 
 */
#include "Data.h"

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

#include "split.h"
#include "helper.h"
#include "Intervals.h"

/**********************  RawData ********************/
RawData::RawData() {
  minX = maxX = 0;
  belongsTo = NULL; 
  cdata = NULL;
}

RawData::RawData(gInterval *v_ginterval) {
  minX = maxX = 0;
  belongsTo = v_ginterval; 
  cdata = NULL;
}

/**
 * @brief Length of the RawData (max-min)
 * Assumes half open coordinates
 * 
 * @return double 
 */
double RawData::Length() {
  return (maxX-minX);   // Half open coordinates means no +1
}

/**
 * @brief Release the memory of raw data
 * This will be useful/necesary for min memory version.
 */
void RawData::ClearData () {
  forward.clear();
  reverse.clear();
}

/**
 * @brief Adds RawData.  Assumes the joint bedGraph conventions, 
 * namely zero based, half open coordinates and the sign of the 
 * coverage indicates strand.  A range of points is input as 
 * each individual point within the set.
 * 
 * @param st  start (this point is included)
 * @param sp   stop (with half open, this is strickly less than!)
 * @param cov  signed count/coverage at this point
 */
void RawData::addDataPoints(double st, double sp, double cov) {
  if (maxX == 0) { minX = st; maxX = sp; } // first data point

  // Adjust min/max as gather data points
  if (st < minX) minX = st;
  if (sp > maxX) maxX = sp;

  // Now add per position coverage info, 
  // Half open coordinates make this strickly less than!
  for (int i = st; i < sp; i++) {
    double c = abs(cov);
    std::vector<double> point {(double)i,c}; 
    // std::cout << "Point: " + to_string(point[0]) + "," + to_string(point[1]) << std::endl;
    if (cov >= 0) {
      forward.push_back(point);
    } else { 
      reverse.push_back(point);
    }
  }
}

/**
 * @brief Once all the data is read in, we can sort the
 * forward and reverse strands and remove duplicate entries.
 * 
 */
void RawData::Sort() {
  // Sort vector<double> by  first entity (e.g. coordinates)
  std::sort(forward.begin(),forward.end(), 
      [](const std::vector<double>& a, const std::vector<double>& b) { 
        return a[0] < b[0]; });
  // std::cout << "After sort: " + data_dump() << std::endl;
  std::sort(reverse.begin(),reverse.end(), 
      [](const std::vector<double>& a, const std::vector<double>& b) { 
        return a[0] < b[0]; });

}

/**
 * @brief Remove duplicate entries for any position.
 * Arbitrarily take one of the entries or max or avg?
 * Erase is notoriously slow, so this could be costly for big datasets.
 */
void RawData::RemoveDuplicates() {
  Sort();
  vector<vector<double>>::iterator it;
  for (int i = 1; i < forward.size(); i++) {
    if (forward[i-1][0] == forward[i][0])  {
      it = forward.begin() + i;
      // This is a duplicate!  Should we throw an error!?!?
      forward.erase(it);
      // Because this shifts the index, need to compare the
      // new i so we'll shift the index.
      i--;
    }
  }  
  for (int i = 1; i < reverse.size(); i++) {
    if (reverse[i-1][0] == reverse[i][0])  {
      it = reverse.begin() + i;
      // This is a duplicate!  Should we throw an error!?!?
      reverse.erase(it);
      i--;
    }
  }  
}

/**
 * @brief  Debuggging output of object
 * 
 * @return std::string 
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
 * @return std::string 
 */
std::string RawData::data_dump() {
  std::string output = "Forward: ";
  for (int i=0; i < forward.size(); i++) {
    output += "[" + tfit::prettyDecimal(forward[i][0],2) + "," 
          + tfit::prettyDecimal(forward[i][1],2) + "]";
  }
  output += "\nReverse: ";
  for (int i=0; i < reverse.size(); i++) {
    output += "[" + tfit::prettyDecimal(reverse[i][0],2) + "," 
      + tfit::prettyDecimal(reverse[i][1], 2) + "]";
  }
  return output;
}


/****************** dInterval *********************/

/**
 * @brief Constructors: dInterval class
 * @author Robin Dowell
 *
 * Purpose: create/allocate instances of a data Interval 
 *
 * The data Interval class contains both strands of data associated
 * with a particular region.  There is an empty constructor option. 
 */
dInterval::dInterval() {
  X = NULL;
  delta = 1;
  scale = 1;
  bins = -1;  // illegal value since we dont know this yet.
  N = 0;

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
  N = 0;  // initial value
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
 * @return std::string 
 */
std::string dInterval::data_dump() {
  std::string output = "Points: ";
  for (int i = 0; i < bins; i ++ ){
    output += "[" + tfit::prettyDecimal(X[0][i],2) + ":" + 
          tfit::prettyDecimal(X[1][i],2) + "," + tfit::prettyDecimal(X[2][i],2) + "]";
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

  
double dInterval::sumAlldata() {
  return N;
}

/**
 * @brief Sum the data on one strand between indicies
 * 
 * @return std::string 
 */
double dInterval::sumInterval(int start, int stop, char strand) {
  double sum = 0;
  int strandidx = 1;
  if (strand == '-') { strandidx = 2;}

  for (int i= start; i< stop; i++) {
    sum += X[strandidx][i];
  }
  return sum;
}

/*****Convert between coordinate systems *****/

/**
 * @brief Convert genomic coordinate into the correct index
 * 
 * @param genomicCoord   genomic Index
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
 * @param dataCoord 
 * @return * int 
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
      int distmin = abs(dataCoord - X[0][indexMin]);
      int distmax = abs(dataCoord - X[0][indexMax]);
      if(distmin < distmax) {
        return indexMin;         
      } else {
        return indexMax;         
      }
    } else if (dataCoord == X[0][indexMin]) {
      found = true;
      return indexMin;
    } else if (dataCoord == X[0][indexMax]) {
      found = true;
      return indexMax;
    }
    indexCenter = (int)(indexMax -indexMin)/2;

    if ((dataCoord > X[0][indexMin]) && (dataCoord < X[0][indexCenter])) {
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
 * @return double 
 */
double dInterval::getGenomefromData(double dataCoord) {
   int index = getIndexfromData(dataCoord);
   return getGenomeCoordfromIndex(index);
}

/**
 * @brief given an index, what is the data coordinate?
 * 
 * @return double 
 */
double dInterval::getDataCoordfromIndex(int index) {
  return X[0][index];
}

/**
 * @brief Given Genomic coordinates, what are data coordinates?
 * 
 * @return double 
 */
double dInterval::getDataCoordfromGenomeCoord(double Gcoord) {
  int index = getIndexfromGenomic(Gcoord);
  return getDataCoordfromIndex(index);
}


/***** Conditioning / Scaling / Binning Data *****/

/**
 * @brief Allocate and setup the X matrix given the size of the RawData
 * 
 * @param length 
 */
void dInterval::initializeData(int minX) {
  // if (X != NULL) {
    // std::cout << "Clearing old data!" << std::endl;
    // ClearX(); 
  // }
  X = new double*[3];
  for (int j = 0 ; j < 3;j++){
    X[j] 		= new double[bins];
  }
  X[0][0] 		= double(minX);
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
    binStart = X[0][i];  binEnd = binStart + delta;
    while ((fi < data->forward.size()) &&     // Still more forward data
           (data->forward[fi][0] >= binStart) &&  // This coordinate is above start
           (data->forward[fi][0] < binEnd)) {   // And below end of range for this bin
      X[1][i] += data->forward[fi][1];    // Add to the bin
      N += data->forward[fi][1];      // Add to the full total (both strands)
      fi++;       // Increment the forward index
    }
    // As above, but now for the reverse strand:
    while ((ri < data->reverse.size()) && 
          (data->reverse[ri][0] >= binStart) && 
          (data->reverse[ri][0] < binEnd)) {
      X[2][i] += data->reverse[ri][1];
      N += data->reverse[ri][1];
      ri++;
    }
  }
}

/**
 * @brief Scales down the coordinates to Zero based.
 * 
 * @param minX 
 */
void dInterval::ScaleDown(int minX) {
  for (int i = 0; i < bins; i ++ ){
    X[0][i] 	= (X[0][i]-minX)/scale;
  }
}

/**
 * @brief Get rid of data points where both strands are zero.
 * 
 */
void dInterval::CompressZeros() {
  int realN 		= 0;	// number of non-zero bins
  for (int i = 0; i < bins;i++){
    if (X[1][i]>0 or X[2][i]>0){
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
    if (X[1][i] > 0 or X[2][i] > 0) {
      newX[0][j] = X[0][i];
      newX[1][j] = X[1][i];
      newX[2][j] = X[2][i];
      j++;
    }
  }
  if (realN != j) {
    printf("WHAT? %d,%d\n", j, realN);
  }
  this->ClearX();
  X = newX;
  bins = realN;
}

/**
 * @brief  This is currently a rewrite of Joey's get_nearest_position
 * into this data structure.
 * 
 * @param position  Typically a mu
 * @param dist    Typically s(1/lambda)
 * @return int    index of the position closest to position + dist (signed dist)
 */
int dInterval::getWithinRangeofPosition(double position, double dist) {
	int i;

	if (dist < 0 ){   // negative strand
		i=0;
		while (i < (bins-1) and (X[0][i] -position) < dist){
			i++;
		}
	}else{
		i=bins-1;
		while (i >0 and (X[0][i] - position) > dist){
			i--;
		}
	}
	return i;
}

/**
 * @brief Cleans up existing X matrix (memory clearing)
 * 
 */
void dInterval::ClearX() {
  for (int i = 0; i < 3; i++) {
    delete X[i];
  }
  delete X;
  X = NULL;
}

