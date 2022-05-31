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


AlgorithmControl::AlgorithmControl() {
	convergence_threshold = 0.0001; //convergence check
	max_iterations = 2000; //stop after this many iterations
	noise_max = 0.05; //fit a uniform noise component, never let it get above this weight
	seed = true; //seed with a gross peak finder
	move_l = false;	// indicator
   r_mu = 0;	
   maxUniformIter = 200;
}


/**********************  The actual algorithm class *************/
EMalg::EMalg(): control(), models() {
   converged = false;	
}

std::string EMalg::write_out() {
   return "No contents!";
}

/**
 * @brief This is the core EM fitting algorithm!
 * 
 * @return int 
 */
int EMalg::fit (dInterval *data) {
   if (models.K == 0) {
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
      /*
         if (data->ForwardCoverage(i)) { // if there is actually data point here...
            norm_forward += components[k].calculateRi(data->Coordinate(i), 1);
         }
         if (data->ReverseCoverage(i)) { // if there is actually data point here...
            norm_reverse += components[k].calculateRi(data->Coordinate(i), -1);
         }
         */
      }  // plus Noise!
      if (norm_forward > 0) {
         models.ll += tfit::LOG(norm_forward) * data->forward(i);
      }
      if (norm_reverse > 0) {
         models.ll += tfit::LOG(norm_reverse) * data->reverse(i);
      }

      // now we need to add the sufficient statistics, need to compute expectations
      //  Equation 9 in Azofeifa 2017
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
 * @param data 
 */
void EMalg::Mstep(dInterval *data) {
		double N=0; //get normalizing constant
      N = models.getAllResponsibilities();
		for (int k = 0; k < models.K; k++){
			// components[k].update_parameters(N, K);
		} // plus noise!
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
