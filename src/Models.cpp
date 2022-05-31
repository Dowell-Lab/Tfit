/**
 * @file Models.cpp
 * @author Robin Dowell 
 * @brief Contains the model objects
 * @version 0.1
 * @date 2022-05-22
 * 
 */
#include "Models.h"

#include <iostream>
#include "helper.h"
#include "Distro.h"

/************** Basic functionality required of all models ************/

BasicModel::BasicModel()
  : sufficiencyStats() {
  weight = 0.; 
}

void BasicModel::resetSufficiency() {
   sufficiencyStats.reset();
}

/***** BIDIRECTIONAL MODEL *****/
// Empty Constructor

Bidirectional::Bidirectional()
  : BasicModel(), loading(), initiation() {
   pi = 0.5;
   footprint = 40;	
}

Bidirectional::Bidirectional(double v_mu, double v_sigma, double v_lambda, 
   double v_pi, double v_footprint) : BasicModel(), loading(), initiation() {
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
 * In Azofeifa 2018, \f$h(z,s;\mu,\sigma,\lambda,\pi)\f$ is defined as:
 * \f$\lambda \phi(\frac{z-\mu}{\sigma}) R(\lambda\sigma - s \frac{z-\mu}{\sigma}) \f$
 * where \f$\phi\f$ is the standard normal distribution.
 * 
 * Wikipedia states the PDF as given in Grushka 1972 as:
 * \f$\frac{\lambda}{2} e^{\frac{\lambda}{2}(2\mu + \lambda\sigma^2 - 2z)}
 *       erfc (\frac{\mu + \lambda\sigma^2 - z}{\sqrt{2}\sigma)})\f$
 *
 * In both cases, the mean difference is multipled by -1 for negative strand
 * and the entire pdf is multipled by an indicator function \f$I = \pi\f$ 
 * if pos strand and \f$I = (1-\pi)\f$ if neg strand.
 * 
 * The Wikipedia version has fewer numerical issues and is used here.
 * @param z current position/read
 * @param s strand 
 * @return double 
 */
double Bidirectional::pdf(double z, char s){
   // Offset the position by the footprint 
   z = applyFootprint(z,s);
	double h;   // Called h(z,s;mu, sigma, lambda, pi) in the paper.

   // Reused intermediates in calculation:
   double pointMeandiff = indicatorStrand(s) * (loading.mu - z);
   double lambdaSigmaSQ = initiation.lambda * pow(loading.sigma,2);
   double halfLambda = initiation.lambda/2.0;
	double exponentvalue = (halfLambda)* (2*(pointMeandiff) + lambdaSigmaSQ);

   // This has s and Indicator added:
	h     = (halfLambda)*exp(exponentvalue)*
      erfc((pointMeandiff + lambdaSigmaSQ)/(sqrt(2)*loading.sigma));
   if (indicatorStrand(s) > 0) {
      h = h*pi; // pow(pi, std::max(0, indicatorStrand(s)) = 0 if s = -1; 1 otherwise 
   } else {
      h = h*(1-pi);
   }

   // Check for outofbounds issues?
	if ((h < 1e7) && (! std::isnan(float(h)))){
	  return h; 
	}
	return 0.0;
}

/**
 * @brief This is the "alternative formulation for computation" of the 
 * exponentially modified gaussian probability density function:
 * \f$f(x;\mu,\sigma,\lambda)\f$ as described in Kalambet et. al. 2011.  
 * Briefly:
 * NOTE: \f$\pi\f$ here is the numerical constant (3.14...) NOT the strand bias.
 * 
 * EQ#1: 
 * \f$\frac{h\sigma}{\tau} \sqrt{\frac{\pi}{2}} e^{0.5(\frac{\sigma}{\tau})^2 - \frac{x-\mu}{\tau}}
 *    erfc(\frac{1}{\sqrt{2}}(\frac{\sigma}{\tau} - \frac{x-\mu}{\sigma})) \f$
 * where:
 * \f$ * h = \frac{1}{\sigma \sqrt{2*\pi}} \f$
 * \f$\tau = \frac{1}{\lambda} \f$
 * 
 * EQ#2:
 * \f$ h e^{-0.5 (\frac{x-\mu}{\sigma})^2} \frac{\sigma}{\tau} \sqrt{\frac{\pi}{2}}
 *    erfcx (\frac{1}{\sqrt{2}} ( \frac{\sigma}{\tau} - \frac{x-\mu}{\sigma}))\f$
 * where: \f$erfcx(t) = exp(t^2) * erfc(t)\f$
 * 
 * EQ#3: 
 * \f$  (h e^{-0.5(\frac{x-\mu}{\sigma})^2}) / (1 + \frac{(x-\mu)\tau}{\sigma^2})\f$
 * 
 * Where DECISION on formula is based on 
 * \f$z = \frac{1}{\sqrt{2}}(\frac{\sigma}{\tau} - \frac{x-\mu}{\sigma})\f$
 * if z < 0 use EQ #1
 * if 0 <= z <= 6.71 x 10^7 use EQ #2
 * if z > 6.71 x 10^7 use EQ#3 
 * 
 * In all cases, the mean difference is multipled by -1 for negative strand
 * and the entire pdf is multipled by an indicator function \f$I = \pi\f$ 
 * if pos strand and \f$I = (1-\pi)\f$ if neg strand, i.e. this is where we 
 * incorporate strand bias (the parameter \f$\pi\f$).
 *     
 * @param z 
 * @param s 
 * @return double 
 */
double Bidirectional::pdf_alt(double z, char s){
  // Offset the position by the footprint 
   z = applyFootprint(z,s);

   double pointMeandiff = indicatorStrand(s) * (z - loading.mu);
   double tau = 1.0/initiation.lambda;

   double decision = 1.0/sqrt(2) * 
      (loading.sigma*initiation.lambda - (pointMeandiff/loading.sigma));
   double h = 1.0/(loading.sigma * sqrt(2*M_PI));

   double f = 0;     // As described in Kalambet et. al. 2011
   if (decision > 6.71e7) {     // Equation #3
      f = h * exp(-0.5 * pow((pointMeandiff/loading.sigma),2));
      f = f / (1.0 + (pointMeandiff*tau) / pow(loading.sigma,2));
   } else {
      double sigmaLambda = loading.sigma*initiation.lambda; // sigma/tau
      double CON = sqrt(M_PI/2); // constant
      double secondTerm = 1.0/sqrt(2) * (sigmaLambda - (pointMeandiff/loading.sigma));
      if (decision < 0) {  // Equation #1
         double firstTerm = exp(0.5 * pow(sigmaLambda,2) - (pointMeandiff/tau));
         f = h * sigmaLambda * CON * firstTerm * erfc(secondTerm);
      } else {  // Equation #2
         double fTerm = exp(-0.5 * pow((pointMeandiff/loading.sigma),2));
         f = h * fTerm * sigmaLambda * CON * erfc(secondTerm) * exp(pow(secondTerm,2));
      }
   }

   // Factor in strand bias
   if (indicatorStrand(s) > 0) {
      f = f*pi; // pow(pi, std::max(0, indicatorStrand(s)) = 0 if s = -1; 1 otherwise 
   } else {
      f = f*(1-pi);
   }
   return f;
}

/**
 * @brief conditional expectation of Y given z_i 
 * Eq 9 E[Y|params] in Azofeifa 2018:
 * \f$ E[Y|g_i,\theta^t] = s_i(z-\mu) - \lambda\sigma^2
 * + \frac{\sigma}{R(\lambda\sigma - s_i(z_i-mu)/\sigma)}
 * \f$
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
 * \f$ E[X|g_i,\theta^t = z_i - s_i E[Y|g_i,\theta^t] \f$
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
 * \f$ E[Y^2|g_i,\theta^t] = \lambda^2\sigma^4 +
 * \sigma^2(2\lambda(\mu-z)s_i+1) + (z_i - \mu)^2 
 * \frac{\sigma(\lambda\sigma^2 + s_i(\mu - z_i))}{R(\lambda\sigma - s_i(z_i-\mu)/\sigma}\f$
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
 * Eq 9 E[X^2|params] in Azofeifa 2018:
 * \f$ E[X^2|g_i,\theta^t] = E[X|g_i,\theta^t] +
 * E[Y^2|g_i,\theta^t] - E[Y|g_i,\theta^t] \f$
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

/*************** Uniform Model (Elongation and Noise) ************************/
UniformModel::UniformModel(): BasicModel(), uni() {
  // Just calls the parents currently.
}

UniformModel::UniformModel(double v_a, double v_b): BasicModel(), uni() {
   uni.lower = v_a;
   uni.upper = v_b;
}

std::string UniformModel::write_out() {
   return ("Uni: " + uni.write_out());
}

/**
 * @brief Uniform Model density function
 * 
 * @param x      
 * @param strand     currently NOT used!
 * @return double 
 */
double UniformModel::pdf(double x, char s){
   // Do we need to do *anything* with strand info?
   return (weight * uni.pdf(x));
}

double UniformModel::getResponsibility() {
   // Strand??
   return (sufficiencyStats.r_forward);
}

/**********************  Full model (with Elongationg) *************/

FullModel::FullModel()
	: bidir(), forwardElongation(), reverseElongation() {
      forwardElongation.weight = 0.5;
      reverseElongation.weight = 0.5;
}

std::string FullModel::write_out() {
   std::string output = bidir.write_out();
   output += " " + forwardElongation.write_out();
   output += " " + reverseElongation.write_out();
   return(output);
}

double FullModel::pdf(double z, char s) {
   // This is weighted by component weights, is this correct?
   return (bidir.weight*bidir.pdf(z,s) 
   + reverseElongation.weight*reverseElongation.pdf(z,s)
   + reverseElongation.weight*reverseElongation.pdf(z,s));
}

void FullModel::resetSufficiencyStats() {
   bidir.resetSufficiency();
   forwardElongation.resetSufficiency();
   reverseElongation.resetSufficiency();
}

double FullModel::getResponsibility() {  // formerly called get_all_repo
   return (bidir.sufficiencyStats.r_forward + bidir.sufficiencyStats.r_reverse
         + forwardElongation.sufficiencyStats.r_forward 
         + reverseElongation.sufficiencyStats.r_reverse);
}
