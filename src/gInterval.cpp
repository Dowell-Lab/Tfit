/**
 * @file gInterval.cpp
 * @author Robin Dowell 
 * @brief Container holding coordinates (bed) and data (dInterval). 
 * 
 * @version 0.1
 * @date 2023-01-17
 * 
 */
#include "gInterval.h"

#include <string>
#include <vector>
#include <iostream>

#include "Bed.h"
#include "Data.h"

gInterval::gInterval() {
    interval = NULL;
    tdf = NULL;
    raw = NULL;
}

std::string gInterval::write_out()
{
   std::string output;
   if (interval != NULL) {
      output = "Interval: \n";
      output += interval->write_out();
   }
   if (tdf != NULL) {
      output = "\nData: \n";
      output += tdf->write_out();
   }
   return output;
}

/**
 * @brief Construct a new dInterval from RawData
 * 
 * @param v_delta   The binning size (nts per bin)
 * @param v_scale   The scaling factor
 */
void gInterval::createDataInterval(int v_delta, int v_scale) {
  if (tdf == NULL) {    // No existing dInterval
    tdf = new dInterval(v_delta, v_scale);
    tdf->seg = (this);
  } else {  // Oh shit, an existing interval, what do we do?

  }
}

void gInterval::raw2dInterval() {
   if (raw == NULL) {
      return; // No raw data, should throw error
   }
   if (tdf == NULL) {
      return; // No dInterval set, should throw error or allocate?
   }
   raw->RemoveDuplicates();      // Is this necessary?  It's potentially time consuming!
   if (tdf->scale <= 0) { tdf->scale = 1;}

   // Establish num of bins
   tdf->setBinNum(raw->Length());
   tdf->initializeData(raw->minX);
   tdf->BinStrands(raw);
   tdf->ScaleDown(raw->minX);
   // std::cout << data_dump() << std::endl;
   tdf->CompressZeros();
   //std::cout << "AFTER: " + to_string(bins) << std::endl;
}

/*****Convert between coordinate systems *****/

/**
 * @brief Convert genomic coordinate into the correct index
 * If no interval exists, assume they meant data coordinate.
 * 
 * @param genomicCoord   genomic coordinate 
 * @return int    correct index into X[]
 */
int gInterval::getIndexfromGenomic(double genomicCoord) {
  double dataCoord;
  if (interval != NULL) {
    dataCoord = genomicCoord - interval->start; // gets zero based coord (e.g. data)
  } else {     // Problem.  Assume they provided data coord.  Should test this.
    dataCoord = genomicCoord; 
  }
  return (tdf->getIndexfromData(dataCoord));
}

/**
 * @brief Convert index to genomic coordinate
 * 
 * Note that because of binning these will always be on 
 * minX + delta intervals.
 * 
 * @param index  position index into X[]
 * @return double   genomic Coordinate cooresponding
 */
double gInterval::getGenomeCoordfromIndex(int index) {
  if ((tdf == NULL) || (raw ==NULL)) { return -1; }   // Error!

  double coord = ((index * tdf->scale) + raw->minX);
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
double gInterval::getGenomefromData(double dataCoord) {
   if (tdf == NULL) { return -1; }  // Error case!
   int index = tdf->getIndexfromData(dataCoord);
   return getGenomeCoordfromIndex(index);
}


/**
 * @brief Given Genomic coordinates, what are data coordinates?
 * 
 * @param genomicCoord   genomic coordinate 
 * @return double Scaled/binned data coordinate
 */
double gInterval::getDataCoordfromGenomeCoord(double Gcoord) {
   if (tdf == NULL) { return -1;} // Error case!
  int index = getIndexfromGenomic(Gcoord);
  return tdf->getDataCoordfromIndex(index);
}