/**
 * @file GenerateReads.h
 * @author Robin Dowell 
 * @version 0.1
 * @date 2022-11-01
 * 
 */
#ifndef GenerateReads_h 
#define GenerateReads_h 

#include <string>
#include <vector>

/**
 * @brief How do we want to generate and summarize reads?
 * 
 */
class ReadGenerator { 
  public:

  // Constructor
  ReadGenerator();

  //Functions
  std::string write_out();

};

#endif
