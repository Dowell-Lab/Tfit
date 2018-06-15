#ifndef model_single_H
#define model_single_H
#include "load.h"
#include "model.h"

using namespace std;

/** Represents an individual nonlinear regression model over a given segment.
 * Much of the functionality within this class is presently unimplemented, though
 * it does appear capable of computing a probability density sample at a given point.
 */ 
class NLR{
public:
	double mu, si, l, pi, wn , wl, wr, fp;
	double EX, EX2, EY, EPI, EPIN,   WE, WF, WR;
	double alpha_1, beta_1, alpha_2, beta_2, alpha_3;
	EMG bidir;
	UNI forward;
	UNI reverse;
	NLR();
	void init( segment *, double, int, double,
	 double, double, double, double, double);
	double addSS(double, double, double);
	double pdf(double );
	double get_all();
	void resetSS();
	void set_new_parameters(double);
	void print();

};

/** Represents a binary classifier that makes judgements based on a single segment of input data.
 * It would appear based on prior documentation that this class was created in order to provide support for the new
 * parameters introduced with the latest development revision of Tfit. Despite this, the classes present in 
 * model_single.h and model_single.cpp are unused in Tfit.
 */
class classifier_single{
public:
	double ll, covergence_threshold, max_iterations;
	int K, type;
	double scale;
	NLR * components;
	double alpha_1, beta_1, alpha_2, beta_2, alpha_3;
	double fit(segment *, double);
	classifier_single();
	classifier_single(double, int, int, int, 
		double, double, double, double, double, double);

};


#endif
