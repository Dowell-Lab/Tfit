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

class HyperParameters {
  //FOR SIGMA ; variance in loading, gamma
	double ALPHA_0, BETA_0;
	//FOR LAMBA ; rate of initiation, gamma
	double ALPHA_1, BETA_1;
	//FOR W ; weight , Dirichlet
	double ALPHA_2;
	//FOR PI ; strand prob. , beta
	double ALPHA_3;

  // Other priors?
	// double p;   // What is this?
	// double foot_print;	//  ?? 

  //Constructor
  HyperParameters();

  //Functions
  std::string write_out();

};

#endif
