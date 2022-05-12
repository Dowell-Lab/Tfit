/**
 * @file Model.cpp
 * @author Robin Dowell 
 * @brief Contains the model objects
 * @version 0.1
 * @date 2022-05-22
 * 
 */
#include "Model.h"

// Empty Constructor
Bidirectional::Bidirectional() { }

Bidirectional::Bidirectional(double v_mu, double v_sigma, double v_lambda, double v_pi, double v_footprint) {
	mu 	= v_mu;
	sigma 	= v_sigma;
	lambda  	= v_lambda;
	pi 	= v_pi;
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
 * @brief Mill's Ratio (formerly just "R")
 * 
 * @param x 
 * @return double 
 */
double Bidirectional::millsRatio(double x){ 
	if (x > 4){
		return 1.0 / x;
	}
	double N = NormalCDF(x);
	double D = NormalPDF(x);
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
      z-=foot_print; 

   }else{
		z+=foot_print;
	}
   return z;
}

/**
 * @brief EMG density function
 * 
 * @param z current position/read
 * @param s strand 
 * @return double 
 */
double Bidirectional::pdf(double z, char s){
   // if (w==0){ return 0.0; }

   // Offset the position by the footprint 
   z = applyFootprint(z,s);

	double vl 		= (lambda/2.0)*(indicatorStrand(s)*2*(mu-z) + lambda*pow(sigma,2));
	double p;
	if (vl > 100){ //potential for overflow, inaccuracies
		p 			= lambda*NormalPDF((z-mu)/sigma)*MillsRatio(lambda*sigma - indicatorStrand(s)*((z-mu)/sigma));
	}else{
		p 			= (lambda/2)*exp(vl)*erfc((indicatorStrand(s)*(mu-z) + lambda*pow(sigma ,2) )/(sqrt(2)*sigma));
	}
	// Doesn't this line negate the whole "if" statement above??
	p     = (lambda/2)*exp(vl)*erfc((indicatorStrand(s)*(mu-z) + lambda*pow(sigma ,2) )/(sqrt(2)*sigma));
	p     = p*w*pow(pi, max(0, indicatorStrand(s)) )*pow(1-pi, max(0, -1*indicatorStrand(s)));

	if (p < pow(10,7) and not isnan(float(p)) ){
	  return p; 
	}
	return 0.0;
}

/**
 * @brief conditional expectation of Y given z_i
 * 
 * @param z 
 * @param s 
 * @return double 
 */
double Bidirectional::ExpY(double z, char s){
   z = applyFootprint(z,s);
	return max(0. , s*(z-mu) - lambda*pow(sigma, 2) + (sigma / MillsRatio(lambda*sigma - s*((z-mu)/sigma))));
}

/**
 * @brief conditional expectation of Y^2 given z_i
 * 
 * @param z 
 * @param s 
 * @return double 
 */
double Bidirectional::ExpY2(double z, char s){
   z = applyFootprint(z,s);
	return pow(lambda,2)*pow(sigma,4) + pow(sigma, 2)*(2*lambda*s*(mu-z)+1 ) + pow(mu-z,2) - ((sigma*(lambda*pow(sigma,2) + s*(mu-z)))/MillsRatio(lambda*sigma - s*((z-mu)/sigma) )); 
}


/******************** Noise Model ************************/

NoiseModel::NoiseModel()
	: noise() {
	
}

std::string NoiseModel::write_out() {
   std::string output;
   output = noise.write_out();	
   return(output);
}

/**********************  Full model (with Elongationg) *************/

FullModel::FullModel()
	: F_bidir(),
	  R_bidir(),
	  forward(),
	  reverse() {
	
}

std::string FullModel::write_out() {
   std::string output;
   output = F_bidir.write_out();	
   output += " " + forward.write_out();	
   output += " " + R_bidir.write_out();	
   output += " " + reverse.write_out();	
   output += " " + to_string(pi);
   return(output);
}




