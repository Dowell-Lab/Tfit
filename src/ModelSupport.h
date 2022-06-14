/**
 * @file ModelSupport.h
 * @author Robin Dowell 
 * @version 0.1
 * @date 2022-05-31
 * 
 */
#ifndef ModelSupport_h
#define ModelSupport_h

#include <string>
#include <vector>

#include "Distro.h"

/**
 * @brief Every model component must be able to report on its 
 * responsiblities, e.g. Prob(k|data point, params) in the EM 
 * algorithm.  
 * 
 */
class Responsibilities {
  public:
  double ri_forward, ri_reverse; //responsibilities per strand
  double r_forward, r_reverse;  // Sum r_i^k 

  //Constructor
  Responsibilities();

  //Functions
  std::string write_out();
  void reset();

  double getResponsibility();
};

/**
 * @brief In the updates to the bidirectional, there are numerous 
 * sum from i to N over all the reads an Exp value * r_i^k.   
 * These are those running sums.
 * 
 */
class sumExpected {
 public:
 double sumRExpY;  // Sum Exp[Y] * r_i^k
 double sumRExpX;  // Sum Exp[X] * r_i^k
 double sumRExpX2;  // Sum Exp[X^2] * r_i^k
 double sumRXExpY;    // Sum (s_i(z_i-mu) - E[Y]) * r_i^k;  Formerly Joey's 'C'

 // Constructor
 sumExpected();

 //Functions
 std::string write_out();
 void reset();

 double getSumRExpY() const { return sumRExpY; }
 double getSumRExpX() const { return sumRExpX; }
 double getSumRExpX2() const { return sumRExpX2; }
 double getSumRXExpY() const { return sumRXExpY; }

 void addToSumRExpY (double);
 void addToSumRExpX (double);
 void addToSumRExpX2 (double);
 void addToSumRXExpY (double);

};

#endif
