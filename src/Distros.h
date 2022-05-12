/**
 * @file Distros.h
 * @author Robin Dowell 
 * @version 0.1
 * @date 2022-05-12
 * 
 */
#ifndef Distros_H
#define Distros_H

#include <string>

/**
 * @brief Uniform distribution class
 * @author Joey Azofeifa  
 * @bug Uses a distinct strand representation from everywhere else in Tfit code.
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

	// Functions
	double pdf(double,int);	

};
*/
class Uniform {
  public:
    double a,b; // lower and upper bounds of uniform, respectively
                
	// Constructors
  Uniform();

	// Functions
  std::string write_out();
};



/**
 * @brief Class for a single instance of the model.  Contains mostly 
 * the EMG but also l from the Uniform. 
 * @author Joey Azofeifa  
class EMG {
public:
	double mu, sigma, lambda, pi, w;
	double foot_print;

	//sufficient stats
  double ri_forward, ri_reverse; //current responsibility
	double ey, ex, ex2, r_forward, r_reverse;//running total
	double ex_r;

	// Internals?
	double C;
	bool move_fp;
	double prev_mu;

	double pdf(double,int);
	double EY(double ,int);
	double EY2(double ,int);
};
 */

class EMG2 {  // Note "2" is not not class with Joey's code.  Will resolve later.
  public:
	  double mu, sigma, lambda, pi;
	  double foot_print;

	// Constructors
  EMG2();
  
  // Functions
  std::string write_out();
}

/* Do we need a normal here too? */

#endif
