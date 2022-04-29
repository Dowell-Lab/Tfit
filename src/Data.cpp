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

#include "split.h"
#include "Intervals.h"

/**********************  RawData ********************/
RawData::RawData() {
  minX = maxX = 0;
  belongsTo = NULL; 
}

double RawData::Length() {
  return (maxX-minX+1); 
}

void RawData::ClearData () {
  forward.clear();
  reverse.clear();
}

void RawData::addDataPoints(double st, double sp, double cov) {
  if (maxX == 0) { minX = st; maxX = sp; } // first data point

  // Adjust min/max as gather data points
  if (st < minX) minX = st;
  if (sp > maxX) maxX = sp;

  // Now add per position coverage info.
  for (int i = st; i <= sp; i++) {
    double c = abs(cov);
    std::vector<double> point {(double)i,c}; 
    if (cov >= 0) {
      forward.push_back(point);
    } else { 
      reverse.push_back(point);
    }
  }
}

std::string RawData::write_out() {
  std::string ID;
  if (belongsTo != NULL ) ID = belongsTo->identifier;
  else ID = "noID";

  std::string output = ID + ":";
  output += to_string(minX) + "-" + to_string(maxX);

  output += "  " + forward.size();
  output += "  " + reverse.size();
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
  delta = v_delta;  scale = v_scale;  bins =  raw->Length()/delta;
  initializeData(raw->Length());
  BinOneStrand(1,raw->forward);
  BinOneStrand(2,raw->reverse);
  ScaleDown(raw->minX);
  CompressZeros();
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

/**
 * @brief Allocate and setup the X matrix given the size of the RawData
 * 
 * @param length 
 */
void dInterval::initializeData(int minX) {
  if (X != NULL) ClearX(); 
  X = new double*[3];
  for (int j = 0 ; j < 3;j++){
    X[j] 		= new double[bins];
  }
  X[0][0] 		= double(minX);
  X[1][0]=0,X[2][0]=0;
	
  for (int i = 1; i < bins; i++){
    X[0][i] 	= X[0][i-1] + delta;
    X[1][i] 	= 0;
    X[2][i] 	= 0;
  }
}

/**
 * @brief   Converts rawData into binned data. 
 * 
 * @param strand  Expects:  1 for forward; 2 for reverse
 */
void dInterval::BinOneStrand(int strand, std::vector<std::vector<double>>sdata) {
  if ((strand >= 1) && (strand <= 2)) { // Only valid strand input
    int j = 0;
    for (int i = 0; i < sdata.size(); i++) {
      while (j < delta and X[0][j] <= sdata[i][0]) {
        j++;
      }
      if (j < delta and sdata[i][0] <= X[0][j]) {
        X[strand][j - 1] += sdata[i][1];
        N += sdata[i][1];
      }
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
 * @brief Position at the xth index.
 * 
 * @arg x The index of a data element.
 * Note that this is not necessarily a genomic position.
 * 
 * @return double 
 */
double dInterval::position(int x) {
  return X[0][x];
}
  
double dInterval::sum_Region() {
  return N;
}

std::string dInterval::write_out() {
  std::string output;
  if (raw != NULL) {
    output = raw->write_out();
  }
  output += "\t" + to_string(bins) + "," + to_string(delta) + "," + to_string(scale);
  return output;
}