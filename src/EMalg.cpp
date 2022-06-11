/**
 * @file EMalg.cpp
 * @author Robin Dowell 
 * @brief Contains the EM algorithm for model inference
 * @version 0.1
 * @date 2022-05-22
 * 
 */
#include "EMalg.h"

#include <algorithm>
#include <cmath>  // isfinite()

#include "helper.h" // nINF and LOG

/**
 * @brief These are all the control variables (knobs) that 
 * can alter the behavior of the EM algorithm.  Reasonable
 * defaults are applied here.
 */
AlgorithmControl::AlgorithmControl() {
	convergence_threshold = 0.0001; //convergence check
	max_iterations = 2000; //stop after this many iterations
	noise_max = 0.05; //fit a uniform noise component, never let it get above this weight
	seed = true; //seed with a gross peak finder
   elon_move = false;
	//move_l = false;	// indicator
   r_mu = 0;	
   maxUniformIter = 200;
}

/**
 * @return std::string contents of object as string
 */
std::string AlgorithmControl::write_out() {
   std::string output;
   output = "Convergence Threshold: " + tfit::prettyDecimal(convergence_threshold,5);
   output += "\nMax Interations: " + tfit::prettyDecimal(max_iterations, 0);
   output += "\nMaximum Background Weight: " + tfit::prettyDecimal(noise_max,2);
   output += "\nMysterious r_mu thingy: " + tfit::prettyDecimal(r_mu,2);
   if (seed) {
     output += "\nSeed the EM " ;
   } else {
     output += "\ndo NOT Seed the EM " ;
   }
   output += "\nMaximum Uniform Iterations: " + tfit::prettyDecimal(maxUniformIter,0);
   return output;
}


/**********************  The actual algorithm class *************/
EMalg::EMalg(): control(), models() {
   converged = false;	
}

/**
 * @return std::string contents of object as string
 */
std::string EMalg::write_out() {
   std::string output;
   if (converged) { output = "Converged!\n"; }
   else { output = "Not Converged!\n";}
   output += control.write_out();
   output += models.write_out();
   return output;
}

/**
 * @brief This is the core EM fitting algorithm!
 * 
 * @return int an indicator of success
 */
int EMalg::fit (dInterval *data) {
   if (models.K == 0) {    // We only have the noise model.
     models.ll = computeBackgroundModel(data);
     return 1;
   } 
   models.initializeComponents(data);

	//======================================================
	int t 			= 0; //EM loop ticker
	double prevll 	= nINF; //previous iterations log likelihood
	int u 			= 0; //elongation movement ticker
	converged 		= false; //has the EM converged?
	while (t < control.max_iterations && not converged){
      models.resetSufficiencyStats();
      Estep(data);
      Mstep(data);

      // Checks if algorithm completed
		if (abs(models.ll-prevll)< control.convergence_threshold){
			converged=true;
		}
		if (not std::isfinite(models.ll)){
			models.ll 	= nINF;
			return 0;	
		}
      if (u > control.maxUniformIter) { 
         adjustBounds();
         u = 0;
      }
		u++;
		t++;
		prevll=models.ll;
   }
   // Cleanup
   return 1;
}

/**
 * @brief computes the likelihood when the NOISE is the only model.
 * 
 * @param data 
 */
double EMalg::computeBackgroundModel(dInterval *data) {
   double ll = 0;
   double pos = data->sumForward();
   double neg = data->sumReverse();
   double pi = pos / (pos + neg); // calculate strand bias from data 

   // setup NOISE component
   models.noise.setPi(pi);

   ll += models.noise.calculateLikelihood(data);
   return ll;
}

/**
 * @brief E-step, grab all the stats and responsibilities
 * Equation 9 from Azofeifa 2017:
 * 
 * \f$ E[Y|g_i,\theta^t] = s_i(z-\mu) - \lambda\sigma^2
 * + \frac{\sigma}{R(\lambda\sigma - s_i(z_i-\mu)/\sigma)} \f$
 * 
 * \f$ E[X|g_i,\theta^t = z_i - s_i E[Y|g_i,\theta^t] \f$
 * 
 * \f$ E[Y^2|g_i,\theta^t] = \lambda^2\sigma^4 +
 * \sigma^2(2\lambda(\mu-z)s_i+1) + (z_i - \mu)^2 
 * \frac{\sigma(\lambda\sigma^2 + s_i(\mu - z_i))}{R(\lambda\sigma - s_i(z_i-\mu)/\sigma}\f$
 * 
 * \f$ E[X^2|g_i,\theta^t] = E[X|g_i,\theta^t] +
 * E[Y^2|g_i,\theta^t] - E[Y|g_i,\theta^t] \f$
 * 
 */
void EMalg::Estep(dInterval *data) {
   double norm_forward, norm_reverse, N; // helper variables
   models.ll = 0;
   // i -> |D| (Azofeifa 2017 pseudocode)
   for (int i = 0; i < data->num_elements(); i++) { // For every data point
      norm_forward = 0;
      norm_reverse = 0;

      // Equation 7 in Azofeifa 2017: calculate r_i^k
      for (int k = 0; k < models.K; k++) { // computing the responsibility terms per model
         if (data->forward(i)) {
           norm_forward += models.setModels[k]->calculateRi(data->getDataCoordfromIndex(i), '+');
         }
         if (data->reverse(i)) {
           norm_reverse += models.setModels[k]->calculateRi(data->getDataCoordfromIndex(i), '-');
         }
      }  // plus Noise!
      norm_forward += models.noise.calculateRi(data->getDataCoordfromIndex(i), '+');
      norm_reverse += models.noise.calculateRi(data->getDataCoordfromIndex(i), '-');

      if (norm_forward > 0) {
         models.ll += tfit::LOG(norm_forward) * data->forward(i);
      }
      if (norm_reverse > 0) {
         models.ll += tfit::LOG(norm_reverse) * data->reverse(i);
      }

      // now we need to add the sufficient statistics 
      for (int k = 0; k < models.K ; k++) {
         /*
         if (norm_forward) {
            components[k].add_stats(data->Coordinate(i), data->ForwardCoverage(i), 1, norm_forward);
         }
         if (norm_reverse) {
            components[k].add_stats(data->Coordinate(i), data->ReverseCoverage(i), -1, norm_reverse);
         }
         */
      } // plus noise!
   }
}

/**
 * @brief M-step, Equation 10 in Azofeifa 2017, Theta_k^(t+1)
 * 
 * Equation 7 from Azofeifa 2017:
 * 
 * \f$ r_i^k=p(k|g_i;\theta_k^g)=\frac{w_k p(g_i;\theta_k^g)}{\sum_{k\in K} 
 * w_k p(g_i;\theta_k^g)} \f$
 * 
 * Equation 10 from Azofeifa 2017:
 *
 * \f$ w^{t+1}_k =\frac{r_k}{r} \f$
 * 
 * \f$ \pi^{t+1}_k =\frac{\sum_{i=1} r_i^k  I(s_i=1) }{r_k} \f$
 * 
 * \f$ \mu^{t+1}_k =\frac{1}{r_k} \sum_{i=1} E[X|g_i;\theta^t ]  r_i^k \f$
 * 
 * \f$ \frac{1}{\lambda^{t+1}_k} =\frac{1}{r_k}  \sum_{i=1} E[Y|g_i; \theta^t] r_i^k \f$
 * 
 * \f$ \sigma^{t+1}_k =\frac{1}{r_k} [ (\sum_{i=1} E[X^2|g_i; \theta^t]r_i^k 
 * -2\mu_k\sum_{i=1} E[X|g_i; \theta^t] r_i^k] +  \mu_k^2   \f$ 
 * 
 * @param data 
 */
void EMalg::Mstep(dInterval *data) {
	double N=0; //get normalizing constant
   N = models.getAllResponsibilities();
	for (int k = 0; k < models.K; k++){
      models.setModels[k]->updateParameters(N,models.K);
   } // plus noise!
   models.noise.updateParameters(N,models.K);   // Joey's may ignore background?
}

/**
 * @brief Adjust the bounds of the elongation components
 * 
 */
void EMalg::adjustBounds () {
   models.SortByMu(); // Joey's: sort_components(components, K);
   if (control.elon_move) { // We are going to move the bounds on the Uniforms
      // update_j_k(components, data, K, N);
      // update_l(components, data, K);
   }
}
