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


/**
 * @brief Wrapper class around the EM
 * 
class classifier{
public:
	int K; //number of components, if K=0 only UNI otherwise K+NOISE components
	double convergence_threshold; //convergence check
	int max_iterations; //stop after this many iterations
	bool seed; //seed with a gross peak finder
	// If noise_max > 0, include a noise component (e.g. it's actually K+1)
	double noise_max; //fit a uniform noise component, never let it get above this weight
	double move; //variance in moving the uniform supports
	//===================================================================================
	//Bayesian Priors
	double p;
	double foot_print;
	
	//===================================================================================
	//final important parameters
	double ll,pi; // log likelihood, 
	double last_diff;
	bool converged;	// indicator
	double r_mu;
	bool move_l;	// indicator
	// These are priors again.
	double ALPHA_0, BETA_0, ALPHA_1, BETA_1, ALPHA_2, ALPHA_3;

	component * components;		// These are all the models we fit
	vector<vector<double>> init_parameters;  // set but never used?

	// Constructor
	classifier(int, double, int, double, double, double, double
		, double, double, double, double, double);
	classifier(int, double, int, double, double, double, double
		, double, double, double, double, bool, double);
    classifier(double , int  , double ,
		double , double , double ,
		double , double , double ,double , vector<vector<double>>, double );
	classifier();

	// Functions
	int fit2(segment *,vector<double>, int, int);  // This is the core EM algorithm
	string write_classifier_setup();  // The invariant parts of the object
	string write_classifier_status();  // The parts that track model quality
	string write_components();
    string write_params();

    void computeUniform (segment * data);	// Used when K = 0
};
 */

class EMalg {
	public:
	// Constructor
	EMalg();

	//Functions
	std::string write_out();

};

#endif
