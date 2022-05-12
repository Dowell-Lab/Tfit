/**
 * @file Model.h
 * @author Robin Dowell 
 * @version 0.1
 * @date 2022-05-12
 * 
 */
#ifndef Model_H
#define Model_H

#include <string>

/**
 * @brief Class that handles a single instance of the model.
 *  Also contains the priors.   
class component{
public:
	EMG bidir;
	UNI forward;
	UNI reverse;
	NOISE noise;

	component * forward_neighbor;
	component * reverse_neighbor;

	// Indicator variables	
	bool type;		// TRUE then have bidir, forward and reverse; FALSE noise
	// bool termination;  // doesn't appear used.
	bool EXIT;

	//=====================================
	//parameters to simulate from
	double alpha_0, alpha_1, alpha_2, beta_0, beta_1, beta_2;
	//=====================================
	//bayesian priors for parameters, MAP
	//FOR SIGMA ; variance in loading, gamma
	double ALPHA_0, BETA_0;
	//FOR LAMBA ; rate of initiation, gamma
	double ALPHA_1, BETA_1;
	//FOR W ; weight , Dirichlet
	double ALPHA_2;
	//FOR PI ; strand prob. , beta
	double ALPHA_3;

	double foot_print;		// set but never used?
	double w_thresh=0;		// set but never used?

	void initialize_bounds(double,  segment *, int , double , double, double, double, double, double);
	double evaluate(double, int);
	void add_stats(double, double , int, double);
	void update_parameters(double,int);
	void set_priors(double,double,double,double,double,double,double, int);
	double get_all_repo();
	void reset();
};
 */

class NoiseModel {
  Uniform noise;

  // Constructor
  NoiseModel();

  // Functions
  std::string write_out();
};

class FullModel {
  public:
	EMG F_bidir;
	EMG R_bidir;
	double pi; 	// strand bias
	Uniform forward;
	Uniform reverse;

  // Constructor
	FullModel();

  // Functions
  std::string write_out();
};


#endif
