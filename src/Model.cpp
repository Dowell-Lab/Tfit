/**
 * @file Model.cpp
 * @author Robin Dowell 
 * @brief Contains the model objects
 * @version 0.1
 * @date 2022-05-22
 * 
 */
#include "Model.h"

#include <iostream>
#include "helper.h"
#include "Distro.h"

// Empty Constructor

Bidirectional::Bidirectional()
	: loading(), initiation() {
   pi = 0.5;
   footprint = 40;	
}

Bidirectional::Bidirectional(double v_mu, double v_sigma, double v_lambda, double v_pi, double v_footprint) {
   loading.mu = v_mu;
   loading.sigma = v_sigma;
   initiation.lambda = v_lambda;
	pi 	= v_pi;     // Check that it's a probability?
   footprint = v_footprint;
}

/**
 * @brief Mill's Ratio (formerly just "R")
 * 
 * @param x 
 * @return double 
 */
double Bidirectional::millsRatio(double x){ 
   // Mill's Ratio asymptotic behavior as x->inf
	if (x > 10){ return 1.0 / x; }

   Normal standard;

   // for smaller x: 
	double N = standard.cdf(x); // normalCDF(x);
	double D = standard.pdf(x); // normalPDF(x);
   // This needs to be adjusted to avoid that pow() call.
	if (D < pow(10,-15)){ //machine epsilon
		return 1.0 / pow(10,-15);
	}
	return exp(log(1. - N)-log(D));
}

/**
 * @brief This is the indicator function on strand
 * 
 * @param s 
 * @return int 
 */
int Bidirectional::indicatorStrand(char s) {
   if (s == '+') {
      return 1;
   } 
   return -1;
}

double Bidirectional::applyFootprint (double z, char s) {
	if (s == '+'){ 
      z-=footprint; 

   }else{
		z+=footprint;
	}
   return z;
}

/**
 * @brief EMG probability density function
 * h(z,s; mu, sigma, lambda, pi) in Azofeifa 2018
 * 
 * @param z current position/read
 * @param s strand 
 * @return double 
 */
double Bidirectional::pdf(double z, char s){
   // if (w==0){ return 0.0; }
   double w = 1;     // dummy placeholder -- no w yet.

   // Offset the position by the footprint 
   z = applyFootprint(z,s);

	double vl 		= (initiation.lambda/2.0)*(indicatorStrand(s)*2*(loading.mu-z) + initiation.lambda*pow(loading.sigma,2));

	double h;   // Called h(z,s;mu, sigma, lambda, pi) in the paper.
	if (vl > 100){ //potential for overflow, inaccuracies
      Normal standard;
		h 			= initiation.lambda*standard.pdf((z-loading.mu)/loading.sigma)*millsRatio(initiation.lambda*loading.sigma - indicatorStrand(s)*((z-loading.mu)/loading.sigma));
	}else{
		h 			= (initiation.lambda/2)*exp(vl)*erfc((indicatorStrand(s)*(loading.mu-z) + initiation.lambda*pow(loading.sigma ,2) )/(sqrt(2)*loading.sigma));
	}
	// Doesn't this line negate the whole "if" statement above??
	h     = (initiation.lambda/2)*exp(vl)*erfc((indicatorStrand(s)*(loading.mu-z) + initiation.lambda*pow(loading.sigma ,2) )/(sqrt(2)*loading.sigma));
	h     = h*w*pow(pi, std::max(0, indicatorStrand(s)) )*pow(1-pi, std::max(0, -1*indicatorStrand(s)));

	if (h < pow(10,7) and not std::isnan(float(h)) ){
	  return h; 
	}
	return 0.0;
}

/**
 * @brief conditional expectation of Y given z_i 
 * Eq 9 E[Y|params] in Azofeifa 2018
 * 
 * @param z 
 * @param s 
 * @return double 
 */
double Bidirectional::ExpY(double z, char strand){
   z = applyFootprint(z,strand);
   double s = indicatorStrand(strand);
	return std::max(0. , s*(z-loading.mu) 
      - initiation.lambda*pow(loading.sigma, 2) 
      + (loading.sigma / millsRatio(initiation.lambda*loading.sigma 
            - s*((z-loading.mu)/loading.sigma))));
}

/**
 * @brief conditional expectation of Y given z_i 
 * Eq 9 E[X|params] in Azofeifa 2018
 * 
 * @param z 
 * @param s 
 * @return double 
 */
double Bidirectional::ExpX(double z, char strand){
   z = applyFootprint(z,strand);
   double s = indicatorStrand(strand);
   return (z - s*(ExpY(z,strand)-footprint));
}

/**
 * @brief conditional expectation of Y^2 given z_i
 * Eq 9 E[Y^2|params] in Azofeifa 2018
 * 
 * @param z 
 * @param s 
 * @return double 
 */
double Bidirectional::ExpY2(double z, char strand){
   z = applyFootprint(z,strand);
   double s = indicatorStrand(strand);
	return pow(initiation.lambda,2)*pow(loading.sigma,4) 
     + pow(loading.sigma, 2)*(2*initiation.lambda*s*(loading.mu-z)+1 ) 
     + pow(loading.mu-z,2) 
     - ((loading.sigma*(initiation.lambda*pow(loading.sigma,2) 
        + s*(loading.mu-z)))/millsRatio(initiation.lambda*loading.sigma 
        - s*((z-loading.mu)/loading.sigma))); 
}

/**
 * @brief conditional expectation of Y^2 given z_i
 * Eq 9 E[X^2|params] in Azofeifa 2018
 * 
 * @param z 
 * @param s 
 * @return double 
 */
double Bidirectional::ExpX2(double z, char strand){
   z = applyFootprint(z,strand);
	return (pow(ExpX(z,strand),2)+ExpY2(z,strand)-ExpY(z,strand));
}

std::string Bidirectional::write_out() {
   std::string output = "Bidir(" + loading.write_out() + ";" + initiation.write_out()
            + ";" + tfit::prettyDecimal(pi,3) + "," + tfit::prettyDecimal(footprint,2) + ")";
   return output;
}

/**
 * @brief Generate n samples from this bidirectional model
 * 
 * @param n          Number of samples to create
 * @return std::vector<double> 
 */
std::vector<double> Bidirectional::generate_data(int n) {
   Random num_gen;

   std::vector<double> results;
   double rnum;  int signStrand = 1;
   for (int i = 0; i < n; i++) {
    // Flip the coin on strand, based on pi (strand bias)
    rnum = num_gen.fetchProbability();
    if (rnum <= pi) {
      signStrand = 1;
    } else {
      signStrand = -1;
    }
    // Given that strand, generate a read from EMG
    double norm = num_gen.fetchNormal(loading.mu,loading.sigma);
    double expon = num_gen.fetchExponential(initiation.lambda);
    results.push_back(signStrand *(norm + signStrand * expon));
   }
   return results;
}


/******************** Noise Model ************************/

NoiseModel::NoiseModel()
	: noise() {
}

NoiseModel::NoiseModel(double v_a, double v_b) {
   noise.lower = v_a;
   noise.upper = v_b;
}

std::string NoiseModel::write_out() {
   return ("Noise: " + noise.write_out());
}

/**
 * @brief NOISE density function
 * Note that NOISE is just a uniform.
 * 
 * @param x          NOT used!
 * @param strand 
 * @return double 
 */
double NoiseModel::pdf(double x, char s){
   double w = 1;     // Place holder!
   double pi = 1;    // Place hodler!
	if (s == '+'){
		return (w*pi) / abs(noise.upper-noise.lower);
	}
	return (w*(1-pi)) / abs(noise.upper-noise.lower);
}

/**********************  Full model (with Elongationg) *************/


FullModel::FullModel()
	: bidir(), forwardElongation(), reverseElongation() {
   w_forward = 0.5;
   w_reverse = 0.5;	
}

std::string FullModel::write_out() {
   std::string output = bidir.write_out();
   return(output);
}




