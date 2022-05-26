/**
 * @file EMalg.h
 * @author Robin Dowell 
 * @version 0.1
 * @date 2022-05-12
 * 
 */
#ifndef EMalg_H 
#define EMalg_H 

#include <string>
#include "Models.h"

class AlgorithmControl { 
  public:
	double convergence_threshold; //convergence check
	int max_iterations; //stop after this many iterations
	bool seed; //seed with a gross peak finder
	double noise_max; //fit a uniform noise component, never let it get above this weight
	bool move_l;	// indicator
	double r_mu;

  AlgorithmControl();
};

class EMalg {
	public:
	bool converged;	// indicator
	// AlgorithmControl *control;  // Do we want to include this as a link?

	// Constructor
	EMalg();

	//Functions
	std::string write_out();

	int fit (); // A data interval, a set of models, algorithm control, priors

};


 /* 
class classifier{
public:
	//===================================================================================
	//Bayesian Priors
	double p;
	double foot_print;
	
	//===================================================================================
	//final important parameters
	double last_diff;
	double r_mu;

	// These are priors again.
	double ALPHA_0, BETA_0, ALPHA_1, BETA_1, ALPHA_2, ALPHA_3;

	component * components;		// These are all the models we fit
	vector<vector<double>> init_parameters;  // set but never used?

 */


#endif
