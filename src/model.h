/**
 * @file model.h
 * @author Joey Azofeifa
 * @brief Header file of \ref UNI (uniform), \ref EMG (full model), \ref NOISE,
 *  \ref component (full set of distros, priors), and \ref classifier (EM convergence
 * info and k number of components) classes.  
 * @version 0.1
 * @date 2016-05-20
 * 
 */
#ifndef model_H
#define model_H

#include <string>

#include "load.h"

using namespace std;

/**
 * @brief Uniform distribution class
 * @author Joey Azofeifa  
 * @bug Uses a distinct strand representation from everywhere else in Tfit code.
 */
class UNI{
public:
	double a,b; // Bounds of the uniform interval
	double w;
	double pi;	// strand (=1 if + strand; =0 if - strand)
	int j,k,l; //so j and k are the bounds the uniform can move through, l is the current one
	// double left_SUM, right_SUM;
	int st;		// Strand -- how is this differnet from pi?
	int pos;	

	//sufficient stats
	double ri_forward, ri_reverse; //current responsibility
	double r_forward, r_reverse; //running total
	double delta_a, delta_b;

	// Constructors
	UNI();
	UNI(double, double, double, int, int, double);

	// Functions
	double pdf(double,int);	
	string write();

};

double sum(double * , int );
double LOG(double );

/**
 * @brief Class for a single instance of the model.  Contains mostly 
 * the EMG but also l from the Uniform. 
 * @author Joey Azofeifa  
 * 
 */
class EMG {
public:
	double mu, sigma, lambda, pi, w;
	double foot_print;

	//sufficient stats
	double ri_forward, ri_reverse; //current responsibility
	double r_forward, r_reverse;//running total
	double ey, ex, ex2; 
	double ex_r;

	// Internals?
	double C;
	bool move_fp;
	double prev_mu;

	// Constructors
	EMG();
	EMG(double, double, double, double, double);
	// Functions
	double pdf(double,int);
	double EY(double ,int);
	double EY2(double ,int);
	string write();
	string write_sufficient_stats();
	string write_internals();
};

/**
 * @brief NOISE class is another Uniform.
 * 
 */
class NOISE{
public:
	double a, b, w, pi;
	double ri_forward, ri_reverse; //current responsibility
	double r_forward, r_reverse; //running total
	double pdf(double, int);

	// Constructors
	NOISE();
	NOISE(double, double, double, double);

	string write();
};

/**
 * @brief Class that handles a single instance of the model.
 *  Also contains the priors.   
 */
class component{
public:
	EMG bidir;
	UNI forward;
	UNI reverse;
	NOISE noise;

	/* double linked list */
	component * forward_neighbor;
	component * reverse_neighbor;

	// Indicator variables	
	bool type;		// TRUE then have bidir, forward and reverse; FALSE noise
	// bool termination;  // doesn't appear used.
	bool EXIT;

	//=====================================
	//parameters to simulate from
	// double alpha_0, alpha_1, alpha_2, beta_0, beta_1, beta_2;
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

	// Constructor
	component();
	// Functions
	void initialize_bounds(double,  segment *, int , double , double, double, double, double, double);
	double evaluate(double, int);
	void add_stats(double, double , int, double);
	void update_parameters(double,int);
	void set_priors(double,double,double,double,double,double,double, int);
	double get_all_repo();
	string write();
	void reset();
	string write_priors();
};

/**
 * @brief Wrapper class around the EM
 * 
 */
class classifier{
public:
	int K; //number of components, if K=0 only UNI otherwise K+NOISE components
	double convergence_threshold; //convergence check
	int max_iterations; //stop after this many iterations
	bool seed; //seed with a gross peak finder
	// If noise_max > 0, include a noise component (e.g. it's actually K+1)
	double noise_max; //fit a uniform noise component, never let it get above this weight
	// double move; 
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



#endif
