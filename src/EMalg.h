/**
 * @file EMalg.h
 * @author Robin Dowell 
 * @version 0.1
 * @date 2022-05-12
 * 
 */
#ifndef EMalg_H 
#define EMalg_H 

#include <string>
#include "ModelSets.h" // ModelContainer
#include "Data.h"		// dInterval
#include "EMseeds.h"

/**
 * @brief Parameters for altering the EM algorithm.
 * Defaults are set in the constructor to this object. 
 * 
 * Expect that these will be overwritten by command line parameters.
 * Expect to output these to configuration files for maximum reproducibility.
 * 
 */
class AlgorithmControl { 
  public:
	double convergence_threshold; //!< convergence check
	int max_iterations; //!< stop after this many iterations
	bool seed; //!< seed with a gross peak finder
	double noise_max; //!< fit a uniform noise component, never let it get above this weight
	// bool move_l;	// indicator
	double r_mu;
	bool elon_move;	//!< should we move the uniform bounds?
	int maxUniformIter;	//!< How many iterations before move uniform?
  
  // Constructor
  AlgorithmControl();

  //Functions
  std::string write_out();

  // Setters
  void setConvergenceThreshold(double v_convergenceThreshold) { convergence_threshold = v_convergenceThreshold; }
  void setMaxIterations(int v_maxIterations) { max_iterations = v_maxIterations; }
  void setNoiseMax(double v_noiseMax) { noise_max = v_noiseMax; }
  void setMaxUniformIter(int v_maxUniformIter) { maxUniformIter = v_maxUniformIter; }
  void setRMu(double v_rMu) { r_mu = v_rMu; }
  // Currently only toggling away from default:
  void TurnOnElongationMove() { elon_move = true;}
  void TurnOffSeeding() { seed = false; }

  // Needs a way to write and read from file.
};

/**
 * @brief The main EM algorithm class.  This is essentially one 
 * function (fit) which fits a set of models (ModelContainer)
 * to a dInterval (data).
 * 
 */
class EMalg {
	public:
	bool converged;	//!< indicator
	AlgorithmControl control;  //!< parameters that can alter the way the EM works
	ModelContainer models;		//!< This is the set of models we're trying to infer
	dInterval *data;	// Data on which to fit this model
	SeedManager seeds;	// Handling seeding the algorithm

	// Constructor
	EMalg();

	//Functions
	std::string write_out();	// debugging

	// Setup Functions
	void setDataAndSeeds(dInterval *v_data);	// Links in data and init seeds

    // The main algorithm: Fitting via EM
	bool fit ();

	// Breaking the EM algorithm into smaller manageable chunks
	// Should these be private?
	double computeBackgroundModel();
	void Estep();
	void Mstep();
	void adjustBounds();
};

 /* 	UNUSED/PORTED parts of Joey's original classifier: 
class classifier{
public:
	//===================================================================================
	//final important parameters
	double last_diff;
	double r_mu;

	// These are priors again.
	double ALPHA_0, BETA_0, ALPHA_1, BETA_1, ALPHA_2, ALPHA_3;

	vector<vector<double>> init_parameters;  // set but never used?

 */


#endif
