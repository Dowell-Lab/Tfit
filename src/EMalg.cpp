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
	converged 		= false; //has the EM converged?
	int u 			= 0; //elongation movement ticker
	while (t < control.max_iterations && not converged){  // Iterate on EM
      models.resetAllSufficiencyStats();
      /*
      if (components[k].EXIT)
      {
         converged = false, ll = nINF;
         return 0;
      }
      */
      Estep(data);
      Mstep(data);

      // Checks if algorithm completed
      // if ((r / N) < pow(10, -5)) EXIT;
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

   return models.noise.calculateLikelihood(data);
}

/**
 * @brief E-step, grab all the stats and responsibilities
 * 1. Calculate Sum_Ri over k (partition function/normalizinng factor) 
 * 2. Use those sums to update log liklihoods
 * 3. Then calculate Expectations (Eqn 9 in Azofeifa) and 
 *    add these to the sum of expectations that is kept (i=1..N)
 */
void EMalg::Estep(dInterval *data) {
   perStrandInfo normalizeRi;
   models.ll = 0;
   // i -> |D| (Azofeifa 2017 pseudocode)
   for (int i = 0; i < data->num_elements(); i++) { // For every data point
      // Calculate Ri and normalize over k. 
      normalizeRi = models.calculateAllRi(i, data);
      // Update log likelihood based on normalizing constants and depth
      if (normalizeRi.forward > 0) {
        models.ll += tfit::LOG(normalizeRi.forward) * data->forward(i);
      }
      if (normalizeRi.reverse > 0) {
        models.ll += tfit::LOG(normalizeRi.reverse) * data->reverse(i);
      }
      // Calculate Expected Values
      models.updateExpectations(i, data, normalizeRi);
  }
}

/**
 * @brief M-step,
 * Update parameters given the expected values and responsibilities
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
