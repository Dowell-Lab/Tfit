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
/****************** sumExpected ********************/

sumExpected::sumExpected() {
   reset(); 
}

void sumExpected::reset() {
   sumRExpX = 0;
   sumRExpY = 0;
   sumRExpX2 = 0;
   sumRXExpY = 0;
}

void sumExpected::addToSumRExpY (double i_term) {
   sumRExpY += i_term;
}

void sumExpected::addToSumRExpX (double i_term) {
  sumRExpX += i_term; 
}

void sumExpected::addToSumRExpX2 (double i_term) {
  sumRExpX2 += i_term; 
}

void sumExpected::addToSumRXExpY (double i_term) {
  sumRXExpY += i_term; 
}

std::string sumExpected::write_out () {
   std::string output;
   output = "Sum Exp[X] * r_i^k: " + tfit::prettyDecimal(sumRExpX, 4);
   output += "Sum Exp[Y] * r_i^k: " + tfit::prettyDecimal(sumRExpY, 4);
   output += "Sum Exp[X2] * r_i^k: " + tfit::prettyDecimal(sumRExpX2, 4);
   output += "Sum (s*(z-mu) - Exp[Y]) * r_i^k: " + tfit::prettyDecimal(sumRXExpY, 4);
   return output;
}

