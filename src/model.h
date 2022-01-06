/**
 * @file model.h
 * @author Joey Azofeifa
 * @brief 
 * @version 0.1
 * @date 2016-05-20
 * 
 */
#ifndef model_H
#define model_H
#include <string>
#include "load.h"
using namespace std;
class UNI{
public:
	double a,b,w,pi;
	int j,k,l; //so j and k are the bounds the uniform can move through, l is the current one
	double left_SUM, right_SUM;
	int st;
	int pos;
	//sufficient stats
	double ri_forward, ri_reverse; //current responsibility
	double r_forward, r_reverse; //running total
	double delta_a, delta_b;
	UNI();
	UNI(double, double, double, int, int, double);
	double pdf(double,int);	
	string print();

};
double sum(double * , int );
double LOG(double );


class EMG{
public:
	double mu, si, l, pi, w;
	//sufficient stats
	double ri_forward, ri_reverse; //current responsibility
	double ey, ex, ex2, r_forward, r_reverse;//running total
	double ex_r;
	double C;
	double foot_print;
	bool move_fp;
	double prev_mu;
	EMG();
	EMG(double, double, double, double, double);
	double pdf(double,int);
	double EY(double ,int);
	double EY2(double ,int);
	string print();
};

/***
class NORMAL{
public:
	double mu,si,pi,w;
	double ex,ex2,r,ri;
	NORMAL();
	NORMAL(double, double, double, double);
	double pdf(double, int);
};
**/

class NOISE{
public:
	double a, b, w, pi;
	double ri_forward, ri_reverse; //current responsibility
	double r_forward, r_reverse; //running total
	double pdf(double, int);
	NOISE();
	NOISE(double, double, double, double);
};


class component{
public:
	EMG bidir;
	UNI forward;
	UNI reverse;
	NOISE noise;
	component * forward_neighbor;
	component * reverse_neighbor;
	
	bool type;
	bool termination;
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

	double foot_print;

	double w_thresh=0;


	component();
	void initialize(double, segment *, int , double , double, double, double);
	void initialize_with_parameters(vector<double>, segment *, int, double, double, double);
	void initialize_with_parameters2(vector<double>, segment *, int, double, double, double);
	void initialize_bounds(double,  segment *, int , double , double, double, double, double, double);
	double evaluate(double, int);
	void add_stats(double, double , int, double);
	double pdf(double , int);
	void update_parameters(double,int);
	void set_priors(double,double,double,double,double,double,double, int);
	double get_all_repo();
	bool check_elongation_support();
	void print();
	void reset();
	string write_out();
};


class classifier{
public:
	int K; //number of components
	double convergence_threshold; //convergence check
	int max_iterations; //stop after this many iterations
	bool seed; //seed with a gross peak finder
	double noise_max; //fit a uniform noise component, never let it get above this weight
	double move; //variance in moving the uniform supports
	//===================================================================================
	//Bayesian Priors
	double p;
	double foot_print;
	int fit(segment *,vector<double>);
	int fit2(segment *,vector<double>, int, int);
	classifier(int, double, int, double, double, double, double
		, double, double, double, double, double);
	classifier(int, double, int, double, double, double, double
		, double, double, double, double, bool, double);
	classifier(double , int  , double ,
		double , double , double ,
		double , double , double ,double , vector<vector<double>>, double );

	classifier();
	void free_classifier();
	string print_out_components();
	int fit_uniform_only(segment * );
	int fit_uniform_only2(segment * );
	
	//===================================================================================
	//final important parameters
	double ll,pi;
	double last_diff;
	component * components;
	bool converged;
	double r_mu;
	bool move_l;
	double ALPHA_0, BETA_0, ALPHA_1, BETA_1, ALPHA_2, ALPHA_3;
	vector<vector<double>> init_parameters;
};







#endif
