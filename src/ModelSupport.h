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
  double r_forward, r_reverse;  // Running totals (sum over i)

  //Constructor
  Responsibilities();

  //Functions
  std::string write_out();
  void reset();

  double getResponsibility();
};

#endif
