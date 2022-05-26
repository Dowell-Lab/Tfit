/**
 * @file EMalg.cpp
 * @author Robin Dowell 
 * @brief Contains the EM algorithm for model inference
 * @version 0.1
 * @date 2022-05-22
 * 
 */
#include "EMalg.h"
#include "helper.h" // nINF


AlgorithmControl::AlgorithmControl() {
	convergence_threshold = 0.0001; //convergence check
	max_iterations = 2000; //stop after this many iterations
	noise_max = 0.05; //fit a uniform noise component, never let it get above this weight
	seed = true; //seed with a gross peak finder
	move_l = false;	// indicator
   r_mu = 0;	
}


/**********************  The actual algorithm class *************/
EMalg::EMalg() {
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
int EMalg::fit () {
   // Compute Uniform Model
   // Initialize Components
   // Seeding (where to start)
   while (1) {
      // Reset Sufficiency Stats
      // E-Step
      // M-Step
   }
   // Cleanup
}
