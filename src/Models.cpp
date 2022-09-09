/**
 * @file Models.cpp
 * @author Robin Dowell 
 * @brief Contains the model objects
 * @version 0.1
 * @date 2022-05-22
 * 
 */
#include "Models.h"
#include <math.h>
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

  betaPi.alpha = 1;
  DirichletW.alpha = 1;

  // These parameter *should* never be used!!!
  betaPi.beta = -1;    // Not used!
  DirichletW.beta = -1;    // Not used!
}

void BasicModel::setPriorPi(double v_alpha) {
   betaPi.alpha = v_alpha;
}

void BasicModel::setPriorWeight(double v_alpha) {
   DirichletW.alpha = v_alpha;
}

void BasicModel::setWeight(double v_weight) {
  weight = v_weight; 
}

void BasicModel::setPi(double v_pi) {
  pi = v_pi; 
}

double BasicModel::getWeight() {
  return weight;
}

double BasicModel::getPi() {
  return pi;
}

double BasicModel::pdf(double z, char st) {
   return 1.;     // Note this function *must* be overwritten by models! 
}

double BasicModel::getResponsibility() {  // formerly called get_all_repo
   return sufficiencyStats.Rk.sumBothStrands();
}

void BasicModel::resetSufficiency() {
   sufficiencyStats.resetRi();
}

/**
 * @brief  Use the responsibility stats to update the parameters (pi and weight)
 *
 * Both parameters are needed for the weight update which has a Dirichlet prior. 
 * @param N    total number of *something*
 * @param K    Number of "components"
 */
void BasicModel::updateParameters(double N, double K) {
   double r = getResponsibility(); // sum Rk 
   pi = (sufficiencyStats.Rk.forward + betaPi.getAlpha()) / (r + betaPi.getAlpha() * 2);

   // Note that this assumes the noise component is equivalently weighted
   // to each model.  Instead perhaps the noise should be Beta weighted 
   // (i.e. restricted to something low) and this renormalized accordingly.
   // Note that 3K is ~ the number of component weights.
   weight = (r + DirichletW.getAlpha()) / (N + DirichletW.getAlpha() * K * 3 + K * 3);
   // std::cout << tfit::prettyDecimal(N+DirichletW.getAlpha()*K*3+K*3,0) << std::endl;
}

/**
 * @brief  Given a position (z), update the position based sufficiency stats.
 * 
 * @param z          position in tranformed coords 
 * @param coverage   coverage per strand
 * @return  updated Ri values (a perStrandInfo)
 */
perStrandInfo BasicModel::calculateRi(double z, perStrandInfo coverage) {
   if (coverage.forward) {
      // std::cout << "Have forward coverage" << std::endl;
      // std::cout << tfit::prettyDecimal(pdf(z,'+'), 2) << std::endl;
      sufficiencyStats.Ri.forward = pdf(z, '+');
   } else {
      // Should this be reset?
   }
   if (coverage.reverse) {
      // std::cout << "Have reverse coverage" << std::endl;
      // std::cout << tfit::prettyDecimal(pdf(z,'-'), 2) << std::endl;
      sufficiencyStats.Ri.reverse = pdf(z, '-');
   } else {
      // Should this be reset?
   }

   return sufficiencyStats.Ri;
}


void BasicModel::updateExpectations(perStrandInfo coverage, perStrandInfo normalizedRi) {
   if (normalizedRi.forward) {
      // std::cout << "Have forward coverage" << std::endl;

      // DO we need to worry about division by ZERO?
     sufficiencyStats.Rk.forward += coverage.forward*sufficiencyStats.Ri.forward/normalizedRi.forward;
      // std::cout << tfit::prettyDecimal(sufficiencyStats.Rk.forward, 2) << std::endl;
   } else {
      // Should Rk be reset?
   }
   if (normalizedRi.reverse) {
      // std::cout << "Have reverse coverage" << std::endl;
     sufficiencyStats.Rk.reverse += coverage.reverse*sufficiencyStats.Ri.reverse/normalizedRi.reverse;
   } else {
      // Should Rk be reset?
   }
   sufficiencyStats.resetRi();
}

/***** BIDIRECTIONAL MODEL *****/
Bidirectional::Bidirectional()
  : BasicModel(), loading(), initiation() {
   // Priors: these are defaults
   footprint = 40;	
   setPriorLambda(1,1);
   setPriorSigma(1,1);
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

double Bidirectional::applyFootprint(double z, char s) {
   // Is z allowed to exceed the bounds?
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

void Bidirectional::setMu(double v_mu) {
  loading.mu = v_mu; 
}

void Bidirectional::setSigma(double v_sigma) {
   loading.sigma = v_sigma;
}

void Bidirectional::setLambda(double v_lambda) {
  initiation.lambda = v_lambda;
}

void Bidirectional::setPriorSigma(double v_alpha, double v_beta) {
   GammaSigma.setAlpha(v_alpha);
   GammaSigma.setBeta(v_beta);
}

void Bidirectional::setPriorLambda(double v_alpha, double v_beta) {
   GammaLambda.setAlpha(v_alpha);
   GammaLambda.setBeta(v_beta);
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
double Bidirectional::pdf(double z, char s) {
   // Offset the position by the footprint 
   z = applyFootprint(z,s);
	double h;   // Called h(z,s;mu, sigma, lambda, pi) in the Azofeifa 2018 paper.

   // Reused intermediates in calculation:
   double pointMeandiff = tfit::StrandAsInt(s) * (loading.mu - z);
   double lambdaSigmaSQ = initiation.lambda * loading.sigma * loading.sigma;
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
      f = h * exp(-0.5 * 
         (pointMeandiff/loading.sigma) * (pointMeandiff/loading.sigma));
      f = f / (1.0 + (pointMeandiff*tau) / (loading.sigma * loading.sigma));
   } else {
      double sigmaLambda = loading.sigma*initiation.lambda; // sigma/tau
      double CON = sqrt(M_PI/2); // constant
      double secondTerm = 1.0/sqrt(2) * (sigmaLambda - (pointMeandiff/loading.sigma));
      if (decision < 0) {  // Equation #1
         double firstTerm = exp(0.5 * (sigmaLambda * sigmaLambda) - (pointMeandiff/tau));
         f = h * sigmaLambda * CON * firstTerm * erfc(secondTerm);
      } else {  // Equation #2
         double meanDivsigma = pointMeandiff/loading.sigma;
         double fTerm = exp(-0.5 * meanDivsigma * meanDivsigma);
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
      - initiation.lambda* (loading.sigma * loading.sigma) 
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
   double sigmaSQ = loading.sigma * loading.sigma;
   double lambdaSQ = initiation.lambda * initiation.lambda;
	return lambdaSQ * sigmaSQ * sigmaSQ 
     + sigmaSQ * (2*initiation.lambda*s*(loading.mu-z)+1) 
     + (loading.mu-z) * (loading.mu-z)
     - ((loading.sigma*(initiation.lambda*sigmaSQ
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
   double ExpXSQ = ExpX(z,strand) * ExpX(z,strand);
	return (ExpXSQ + ExpY2(z,strand)-ExpY(z,strand));
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

void Bidirectional::updateExpectations(double z, perStrandInfo coverage, 
                                          perStrandInfo normalizeRi) {
   BasicModel::updateExpectations(coverage,normalizeRi);

	//now adding all the conditional expectations for the convolution, per strand
   double rterm = sufficiencyStats.Ri.forward / normalizeRi.forward;
   if (rterm > 0 and coverage.forward > 0) {
      calcExpectedVals(z, '+', rterm*coverage.forward);
   }
   rterm = sufficiencyStats.Ri.reverse / normalizeRi.reverse;
   if (rterm > 0 and coverage.reverse > 0) {
      calcExpectedVals(z, '-', rterm*coverage.reverse);
   }
}

void Bidirectional::calcExpectedVals (double position, char strand, double rtimescov) {
   double current_EY = ExpY(position, strand);
   sumOverN.addToSumRExpY(current_EY * rtimescov); 

   double current_EX = ExpX(position,strand);
   sumOverN.addToSumRExpX(current_EX * rtimescov);

   double current_EY2 = ExpY2(position, strand);
   sumOverN.addToSumRExpX2((current_EX*current_EX + current_EY2 - current_EY*current_EY)*rtimescov);

   sumOverN.addToSumRXExpY(std::max((strand * (position - getMu()) - current_EY) *rtimescov, 0.0)); 
}

void Bidirectional::updateParameters(double N, double K) {
   BasicModel::updateParameters(N,K);  // Updates pi and weight
   double r = getResponsibility();
   double prevmu = getMu();

  setMu(sumOverN.sumRExpX/ (r + 0.001)); 
   // Note that Joey uses bidir.mu which was set above to be the t+1 instance.
   // Yet the updates should be using the previous step (mu_t).  I believe
   // this is a bug in the original Tfit code.
   // Also notice here he is taking the 0.5 power on this (sqrt), which is NOT in the 
   // original paper.
   double tempSigma = sqrt(abs((1. / (r + 3 + GammaSigma.getAlpha())) 
                  * (sumOverN.sumRExpX2 - 2 * prevmu * sumOverN.sumRExpX +
                  r * prevmu * prevmu + 2 * GammaSigma.getBeta())));
   setSigma(tempSigma);

   // Note that setting the max at 0.05 is equivalent to setting tau (in bps; 1/lambda)
   // to a upper limit of ~20 bases.  This should be strongly biasing the EMG
   // to look Gaussian!   Ugh.
   setLambda(std::min((r + GammaLambda.getAlpha()) / 
            (sumOverN.sumRExpY + GammaLambda.getBeta()), 5.));
   setLambda(std::max(getLambda(),0.05));

	if (abs(getMu() - prevmu) < 0.01) {
		//bidir.move_fp 	= true;
	} else{
		// bidir.prev_mu 	= bidir.mu;
	}
   /*
	if (bidir.move_fp){
		footprint 	= std::min( std::max(sumOverN.sumRXExpY / (r+0.1),0.0) , 2.5);
	}
   */
}

void Bidirectional::setParametersModel(double v_mu, double v_sigma,
                        double v_lambda, double v_weight) {
  
   setMu(v_mu);  
   setSigma(v_sigma);
   setLambda(v_lambda);
   setWeight(v_weight);
   setPi(0.5);
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

void UniformModel::setBounds(double v_lower, double v_upper) {
   uni.lower = v_lower;
   uni.upper = v_upper;
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
double UniformModel::pdf(double x, char s) {
   // Should the strand use pi?
   return (weight * uni.pdf(x));
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

void UniformModel::updateParameters(double N, double K) {
   BasicModel::updateParameters(N,K);
}

void UniformModel::initalizeBounds(double v_minX, double v_maxX, double v_weight, double v_pi) {
    setBounds(v_minX,v_maxX);
    setPi(v_pi);
    setWeight(v_weight);
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
   return (bidir.getResponsibility() + forwardElongation.sufficiencyStats.Rk.forward 
         + reverseElongation.sufficiencyStats.Rk.reverse);
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

perStrandInfo FullModel::calculateRi(double z, perStrandInfo coverage) {
   perStrandInfo fromBidir = bidir.calculateRi(z,coverage);
   perStrandInfo fromForward = forwardElongation.calculateRi(z,coverage);
   perStrandInfo fromReverse = reverseElongation.calculateRi(z,coverage);

   perStrandInfo cummulativeRi;
   cummulativeRi.forward = fromBidir.forward + fromForward.forward + fromReverse.forward;
   cummulativeRi.reverse = fromBidir.reverse + fromForward.reverse + fromReverse.reverse;

   return cummulativeRi;
}

void FullModel::updateExpectations(double z, perStrandInfo coverage, perStrandInfo normalizeRi) {
  bidir.updateExpectations(z, coverage, normalizeRi);
  forwardElongation.updateExpectations(coverage, normalizeRi);
  reverseElongation.updateExpectations(coverage, normalizeRi);
}

void FullModel::initBounds(double v_mu, double v_sigma, double v_lambda, 
               double v_weight, double v_minX, double v_maxX) {
   double tau = 1/ v_lambda;
   forwardElongation.initalizeBounds(v_mu + tau, v_maxX, v_weight, 1.0);
   bidir.setParametersModel(v_mu, v_sigma, v_lambda, v_weight);
   reverseElongation.initalizeBounds(v_minX, v_mu - tau, v_weight, -1.0);
}
