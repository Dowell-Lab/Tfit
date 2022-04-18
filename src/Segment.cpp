/**
 * @file Segment.cpp
 * @author Robin Dowell 
 * @brief Class code for segments
 * Segments contain both a genomic interval (BED) and a data 
 * interval (bedGraph).
 *
 * This will replace the segment class in load.cpp eventually.
 * 
 * @version 0.1
 * @date 2022-04-18
 * 
 */
#include "Segment.h"

#include <string>
#include <vector>
#include <iostream>

Segment::Segment() {

}

std::string Segment::write_out() {
   std::string output;
   output = coords.write_out();  
   output += "\n" + data.write_out();

   return output;
}

std::string Segment::write_Interval() {
   return coords.write_out(); 
}

std::string Segment::write_data() {
   return data.write_out(); 
}

