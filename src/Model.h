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
#include <vector>

class Bidirectional {
  public:
	double mu, sigma, lambda;	// These are "tied" between the two strands
	double pi;		// strand bias
	double footprint;		// the ad hoc footprint parameter 

  // Constructor
  Bidirectional();
  Bidirectional(double, double, double, double, double);

  // Functions
  std::string write_out();

  double pdf(double z, char s);
  double ExpY(double z, char s);
  double ExpY2(double z, char s);

  std::vector<double> generate_data(int n, char s);

// private:
  double normalPDF(double);
  double millsRatio(double);
  double normalCDF(double x);
  int indicatorStrand(char s);
  double applyFootprint (double z, char s);
};


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
	//FOR LAMBDA ; rate of initiation, gamma
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
	double a, b;

  // Constructor
  NoiseModel();
  NoiseModel(double v_a, double v_b);

  // Expectation: pG(b-a)/S where S is length of genome, G is # reads mapped, p is prob noise

  // Functions
  std::string write_out();

  double pdf(double x, char s);

};

class FullModel {
  public:
    Bidirectional bidir;
	// Uniform forward;
	// Uniform reverse;

  // Constructor
	FullModel();

  // Functions
  std::string write_out();
};


#endif
