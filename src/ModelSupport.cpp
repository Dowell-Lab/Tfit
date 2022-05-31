/**
 * @file ModelSupport.cpp
 * @author Robin Dowell 
 * @brief Contains the responsibilities, hyper parameters and other supporting objects 
 * @version 0.1
 * @date 2022-05-31
 * 
 */
#include "ModelSupport.h"

#include <iostream>
#include "helper.h"
#include "Distro.h"

/********* Sufficiency Statistics *****/

Responsibilities::Responsibilities() {
   r_forward = 0;
	r_reverse = 0;
   ri_forward = 0;
	ri_reverse = 0;
}

void Responsibilities::reset() {
   ri_forward = 0;
	ri_reverse = 0;
}

std::string Responsibilities::write_out() {
   std::string output;
   output = "Ri: f: " + tfit::prettyDecimal(ri_forward,3);
   output += " r: " + tfit::prettyDecimal(ri_reverse,3);
   output = "Rk: f: " + tfit::prettyDecimal(r_forward,3);
   output += " r: " + tfit::prettyDecimal(r_reverse,3);
   return output;
}

double Responsibilities::getResponsibility() {
   return (r_forward + r_reverse);
}

/************** HyperParmaeters *******************/

HyperParameters::HyperParameters() {
  ALPHA_0 = 1;
  ALPHA_1 = 1;
  ALPHA_2 = 1;
  ALPHA_3 = 1;
  BETA_0 = 1;
  BETA_1 = 1;
}

std::string HyperParameters::write_out() {
  std::string output;
  output = "For sigma: Gamma(" + tfit::prettyDecimal(ALPHA_0, 4) + "," 
                  + tfit::prettyDecimal(BETA_0,4) + ")";
  output = "\nFor lambda: Gamma(" + tfit::prettyDecimal(ALPHA_1, 4) + "," 
                  + tfit::prettyDecimal(BETA_1,4) + ")";
  output = "\nFor weight: Dirichlet(" + tfit::prettyDecimal(ALPHA_2, 4) + ")" ;
  output = "\nFor pi: Beta(" + tfit::prettyDecimal(ALPHA_3, 4) + ")";
  return output;
}

