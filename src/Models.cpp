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
#include <algorithm>  // min and max
#include "helper.h"     // tfit::prettyDecimal tift::StrandAsInt
#include "Distro.h"     //Normal, Uniform, Exponential
#include "Data.h"    // dInterval

/************** Basic functionality required of all models ************/

BasicModel::BasicModel()
  : sufficiencyStats() {
  weight = 0.; 
  pi = 0.5;

  alpha_pi = 1;
  alpha_w = 1;
}

void BasicModel::setPriorPi(double v_alpha) {
   alpha_pi = v_alpha;
}

void BasicModel::setPriorWeight(double v_alpha) {
   alpha_w = v_alpha;
}

void BasicModel::resetSufficiency() {
   sufficiencyStats.reset();
}

/***** BIDIRECTIONAL MODEL *****/
// Empty Constructor

Bidirectional::Bidirectional()
  : BasicModel(), loading(), initiation() {
   footprint = 40;	
   // Priors:
   alpha_sigma = 1;
   beta_sigma = 1;
   alpha_lambda = 1;
   beta_lambda = 1;
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

double Bidirectional::applyFootprint (double z, char s) {
	if (s == '+'){ 
      z-=footprint; 

   }else{
		z+=footprint;
	}
   return z;
}

double Bidirectional::getMu() {
  return (loading.mu);
}

double Bidirectional::getSigma() {
  return (loading.sigma);
}

double Bidirectional::getLambda() {
  return (initiation.lambda);
}

void Bidirectional::setMu(double newmu) {
  loading.mu = newmu; 
}

void Bidirectional::setSigma(double newsigma) {
   loading.sigma = newsigma;
}

void Bidirectional::setLambda(double newlambda) {
  initiation.lambda = newlambda;
}

void Bidirectional::setPriorSigma(double v_alpha, double v_beta) {
   alpha_sigma = v_alpha;
   beta_sigma = v_beta; 
}

void Bidirectional::setPriorLambda(double v_alpha, double v_beta) {
   alpha_lambda = v_alpha;
   beta_lambda = v_beta; 
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
   double pointMeandiff = tfit::StrandAsInt(s) * (loading.mu - z);
   double lambdaSigmaSQ = initiation.lambda * pow(loading.sigma,2);
   double halfLambda = initiation.lambda/2.0;
	double exponentvalue = (halfLambda)* (2*(pointMeandiff) + lambdaSigmaSQ);

   // This has s and Indicator added:
	h     = (halfLambda)*exp(exponentvalue)*
      erfc((pointMeandiff + lambdaSigmaSQ)/(sqrt(2)*loading.sigma));
   if (tfit::StrandAsInt(s) > 0) {
      h = h*pi; // pow(pi, std::max(0, tfit::StrandAsInt(s)) = 0 if s = -1; 1 otherwise 
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

   double pointMeandiff = tfit::StrandAsInt(s) * (z - loading.mu);
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
   if (tfit::StrandAsInt(s) > 0) {
      f = f*pi; // pow(pi, std::max(0, tfit::StrandAsInt(s)) = 0 if s = -1; 1 otherwise 
   } else {
      f = f*(1-pi);
   }
   return f;
}

/**
 * @brief conditional expectation of Y given z_i 
 * Eq 9 E[Y|params] in Azofeifa 2018:
 * \f$ E[Y|g_i,\theta^t] = s_i(z-\mu) - \lambda\sigma^2
 * + \frac{\sigma}{R(\lambda\sigma - s_i(z_i-\mu)/\sigma)}
 * \f$
 * 
 * @param z       position
 * @param strand  strand ('+' '-' or '.')
 * @return double  Expected value of Y
 */
double Bidirectional::ExpY(double z, char strand){
   z = applyFootprint(z,strand);
   double s = tfit::StrandAsInt(strand);
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
 * @param z       position
 * @param strand  strand ('+' '-' or '.')
 * @return double  Expected value of X
 */
double Bidirectional::ExpX(double z, char strand){
   z = applyFootprint(z,strand);
   double s = tfit::StrandAsInt(strand);
   return (z - s*(ExpY(z,strand)-footprint));
}

/**
 * @brief conditional expectation of Y^2 given z_i
 * Eq 9 E[Y^2|params] in Azofeifa 2018
 * \f$ E[Y^2|g_i,\theta^t] = \lambda^2\sigma^4 +
 * \sigma^2(2\lambda(\mu-z)s_i+1) + (z_i - \mu)^2 
 * \frac{\sigma(\lambda\sigma^2 + s_i(\mu - z_i))}{R(\lambda\sigma - s_i(z_i-\mu)/\sigma}\f$
 * 
 * @param z       position
 * @param strand  strand ('+' '-' or '.')
 * @return double Expected value of \f$ Y^2 \f$
 */
double Bidirectional::ExpY2(double z, char strand){
   z = applyFootprint(z,strand);
   double s = tfit::StrandAsInt(strand);
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
 * @param z       position
 * @param strand  strand ('+' '-' or '.')
 * @return double Expected value of \f$ X^2 \f$
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

double Bidirectional::getResponsibility() {  // formerly called get_all_repo
   return (sufficiencyStats.r_forward + sufficiencyStats.r_reverse);
}

void Bidirectional::updateParameters(double N, double K) {
   double r = getResponsibility();
   double prevmu = getMu();

   pi = (sufficiencyStats.r_forward + alpha_pi) / (r + alpha_pi * 2);
   // Note that this assumes the noise component is equivalently weighted
   // to each model.  Instead perhaps the noise should be Beta weighted 
   // (i.e. restricted to something low) and this renormalized accordingly.
   // Note that 3K is ~ the number of component weights.
   weight = (r + alpha_w) / (N + alpha_w * K * 3 + K * 3);
/* setMu(bidir.ex / (r + 0.001)); 
   // Note that Joey uses bidir.mu which was set above to be the t+1 instance.
   // Yet the updates should be using the previous step (mu_t).  I believe
   // this is a bug in the original Tfit code.
   double tempSigma = (pow(abs((1. / (r + 3 + alpha_sigma)) * (bidir.ex2 - 2 * bidir.mu * bidir.ex +
                                                     r * pow(bidir.mu, 2) + 2 * beta_sigma)), 0.5);
   setSigma(tempSigma);
   setLambda(min((r + alpha_lambda) / (bidir.ey + beta_lambda), 5.));
*/

   // Note that setting the max at 0.05 is equivalent to setting tau (in bps; 1/lambda)
   // to a upper limit of ~20 bases.  This should be strongly biasing the EMG
   // to look Gaussian!   Ugh.
   setLambda(std::max(getLambda(),0.05));
   /* This is not exactly Joey's logic.  See update_parameters() for the
   interplay of bidir.move_fp, bidir_prev_mu!!! */
   if (abs(getMu() - prevmu) < 0.01) {
    // footprint = min(max(bidir.C / (r + 0.1), 0.0), 2.5);
   }
}

double Bidirectional::calculateRi(double z, char strand) {
   sufficiencyStats.ri_forward = pdf(z,strand);
   sufficiencyStats.ri_reverse = pdf(z,strand);

   return (sufficiencyStats.ri_forward + sufficiencyStats.ri_reverse);
}

/*************** Uniform Model (Elongation and Noise) ************************/
UniformModel::UniformModel(): BasicModel(), uni() {
  pi = 0.;
}

UniformModel::UniformModel(double v_a, double v_b): BasicModel(), uni() {
   uni.lower = v_a;
   uni.upper = v_b;
   pi = 0.;
}

std::string UniformModel::write_out() {
   return ("Uni: " + uni.write_out() + " pi: " + tfit::prettyDecimal(pi,4));
}

/**
 * @brief Uniform Model density function
 * 
 * @param x      
 * @param strand     currently NOT used!
 * @return double 
 */
double UniformModel::pdf(double x, char s){
   // Should the strand use pi?
   return (weight * uni.pdf(x));
}

void UniformModel::setPi(double v_pi) {
   pi = v_pi;
}

double UniformModel::calculateLikelihood(dInterval *data) {
   double ll = 0;
   // Calculate MLE = -n log (pi/l) where pi is strand (1 -> +; -1 -> -)
   for (int i = 0; i < data->num_elements(); i++) {
      if (pi > 0) {
            ll += log(pi / data->getLength()) * data->forward(i);
      }
      if (pi < 1) {
            ll += log((1 - pi) / data->getLength()) * data->reverse(i);
      }
   }
   return ll;
}

void UniformModel::setBounds(double v_lower, double v_upper) {
   uni.lower = v_lower;
   uni.upper = v_upper;
}

double UniformModel::getResponsibility() {
   return (sufficiencyStats.r_forward + sufficiencyStats.r_reverse);
}

void UniformModel::updateParameters(double N, double K) {
   weight = (sufficiencyStats.r_forward + alpha_w) / (N + alpha_w*K*3 + K*3);
   pi = (sufficiencyStats.r_forward + 1) / 
      (sufficiencyStats.r_forward + sufficiencyStats.r_reverse + 2);
}

double UniformModel::calculateRi(double z, char strand) {
   sufficiencyStats.ri_forward = pdf(z,strand);
   sufficiencyStats.ri_reverse = pdf(z,strand);

   return (sufficiencyStats.ri_forward + sufficiencyStats.ri_reverse);
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

void FullModel::updateParameters(double N, double K) {
   bidir.updateParameters(N,K);
   forwardElongation.updateParameters(N,K);
   reverseElongation.updateParameters(N,K);

   // Updates necessary to keep this properly tied together:
   forwardElongation.uni.lower = bidir.getMu();
   reverseElongation.uni.upper = bidir.getMu();

   // Do these get set to zero and one?  If not why?
   if (bidir.weight == 0) {
      forwardElongation.weight = 0;
      reverseElongation.weight = 0;
   }
}

double FullModel::calculateRi(double z, char strand) {
   bidir.calculateRi(z,strand);
   forwardElongation.calculateRi(z,strand);
   reverseElongation.calculateRi(z,strand);

	return (bidir.sufficiencyStats.ri_forward 
      + forwardElongation.sufficiencyStats.ri_forward
      + reverseElongation.sufficiencyStats.ri_forward);
   /*
	return bidir.ri_reverse + reverse.ri_reverse + forward.ri_reverse;
   */
}


