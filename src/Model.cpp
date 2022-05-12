/**
 * @file Model.cpp
 * @author Robin Dowell 
 * @brief Contains the model objects
 * @version 0.1
 * @date 2022-05-22
 * 
 */
#include "Model.h"

NoiseModel::NoiseModel()
	: noise() {
	
}

std::string NoiseModel::write_out() {
   std::string output;
   output = noise.write_out();	
   return(output);
}

FullModel::FullModel()
	: F_bidir(),
	  R_bidir(),
	  forward(),
	  reverse() {
	
}

std::string FullModel::write_out() {
   std::string output;
   output = F_bidir.write_out();	
   output += " " + forward.write_out();	
   output += " " + R_bidir.write_out();	
   output += " " + reverse.write_out();	
   output += " " + to_string(pi);
   return(output);
}






