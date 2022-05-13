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

// Empty Constructor
Bidirectional::Bidirectional() { }

Bidirectional::Bidirectional(double v_mu, double v_sigma, double v_lambda, double v_pi, double v_footprint) {
	mu 	= v_mu;
	sigma 	= v_sigma;
	lambda  	= v_lambda;
	pi 	= v_pi;     // Check that it's a probability?
   footprint = v_footprint;
}

/** 
 * @brief the standard Normal PDF
 * @param x  position
 * @return double 
 */
double Bidirectional::normalPDF(double x) { 
	return exp(-pow(x,2)*0.5)/sqrt(2*M_PI);
}
/**
 * @brief the standard Normal CDF
 * @param x 
 * @return double 
 */
double Bidirectional::normalCDF(double x){ 
	return 0.5*(1+erf(x/sqrt(2)));
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

   // for smaller x: 
	double N = normalCDF(x);
	double D = normalPDF(x);
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

	double vl 		= (lambda/2.0)*(indicatorStrand(s)*2*(mu-z) + lambda*pow(sigma,2));

	double h;   // Called h(z,s;mu, sigma, lambda, pi) in the paper.
	if (vl > 100){ //potential for overflow, inaccuracies
		h 			= lambda*normalPDF((z-mu)/sigma)*millsRatio(lambda*sigma - indicatorStrand(s)*((z-mu)/sigma));
	}else{
		h 			= (lambda/2)*exp(vl)*erfc((indicatorStrand(s)*(mu-z) + lambda*pow(sigma ,2) )/(sqrt(2)*sigma));
	}
	// Doesn't this line negate the whole "if" statement above??
	h     = (lambda/2)*exp(vl)*erfc((indicatorStrand(s)*(mu-z) + lambda*pow(sigma ,2) )/(sqrt(2)*sigma));
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
	return std::max(0. , s*(z-mu) - lambda*pow(sigma, 2) 
      + (sigma / millsRatio(lambda*sigma - s*((z-mu)/sigma))));
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
	return pow(lambda,2)*pow(sigma,4) + pow(sigma, 2)*(2*lambda*s*(mu-z)+1 ) 
     + pow(mu-z,2) 
     - ((sigma*(lambda*pow(sigma,2) + s*(mu-z)))/millsRatio(lambda*sigma - s*((z-mu)/sigma))); 
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
   std::string output = "Bidir(" + tfit::prettyDecimal(mu,2) + "," + tfit::prettyDecimal(sigma,2)
            + "," + tfit::prettyDecimal(lambda,4) + "," + tfit::prettyDecimal(pi,3) + "," 
            + tfit::prettyDecimal(footprint,2) + ")";
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
    double norm = num_gen.fetchNormal(mu,sigma);
    double expon = num_gen.fetchExponential(lambda);
    results.push_back(signStrand *(norm + signStrand * expon));
   }
   return results;
}


/******************** Noise Model ************************/

NoiseModel::NoiseModel(){} //empty constructor

NoiseModel::NoiseModel(double v_a, double v_b) {
	a=v_a;
	b=v_b;
}

std::string NoiseModel::write_out() {
   std::string output = "Noise: U(" + tfit::prettyDecimal(a,4) + "," 
      + tfit::prettyDecimal(b,4) + ")";
   return(output);
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
		return (w*pi) / abs(b-a);
	}
	return (w*(1-pi)) / abs(b-a);
}



/**********************  Full model (with Elongationg) *************/

FullModel::FullModel() {
	
}

std::string FullModel::write_out() {
   std::string output = bidir.write_out();
   return(output);
}




