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

/********* Generic Priors *****/
Priors::Priors() {
   alpha = 1;
   beta = 1; 
}

std::string Priors::write_out() {
   std::string output;
   output = "Prior(" + tfit::prettyDecimal(alpha,2) 
         + "," + tfit::prettyDecimal(beta,2) + ")";
   return output;
}

/********* Container for per Strand Info *****/

perStrandInfo::perStrandInfo() {
   forward = 0;
   reverse = 0; 
}

std::string perStrandInfo::write_out() {
   std::string output;
   output = "Forward: " + tfit::prettyDecimal(forward,4) 
         + " Reverse: " + tfit::prettyDecimal(reverse,4);
   return output;
}

double perStrandInfo::sumBothStrands() {
   return (forward + reverse);
}

/********* Sufficiency Statistics *****/

Responsibilities::Responsibilities(): Ri(), Rk() {
}

void Responsibilities::reset() {
   Ri.forward = 0.;
   Ri.reverse = 0.;
}

std::string Responsibilities::write_out() {
   std::string output;
   output = "Ri: f: " + tfit::prettyDecimal(Ri.forward,3);
   output += " r: " + tfit::prettyDecimal(Ri.forward,3);
   output = "Rk: f: " + tfit::prettyDecimal(Rk.forward,3);
   output += " r: " + tfit::prettyDecimal(Rk.reverse,3);
   return output;
}

double Responsibilities::getResponsibility() {
   return Rk.sumBothStrands();
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

/***********  Constraints on parameters ***************/

bidirConstraints::bidirConstraints() {
  lambdaMin = 0.05;
  lambdaMax = 5; // Constrain lambda so it doesn't eat the elongation
  fp_Min = 0;
  fp_Max = 2.5;  // Constrain footprint so that we aren't looking at adjacent ones.
}

std::string bidirConstraints::write_out() {
   std::string output;
   output = "Lambda: Min: " + tfit::prettyDecimal(lambdaMin,4) 
         + " Max: " + tfit::prettyDecimal(lambdaMax, 4);
   output += " Footprint: Min: " + tfit::prettyDecimal(fp_Min, 3)
      + " Max: " + tfit::prettyDecimal(fp_Max, 3);
   return output;
}
