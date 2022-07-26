/**
 * @file ModelSets.cpp
 * @author Robin Dowell 
 * @brief Contains the model objects
 * @version 0.1
 * @date 2022-05-26
 * 
 */
#include "ModelSets.h"

#include <iostream>
#include <algorithm> // for sort
#include "helper.h"
#include "Data.h"

/**********************  Model Wrapper (for EM) *************/

ModelWrapper::ModelWrapper() {
   type = FMOD;         
   gene = NULL;
   bidir = NULL;
}

ModelWrapper::ModelWrapper(FullModel *v_gene) {
   type = FMOD;
   gene  = v_gene;
   bidir = NULL;
   gene->resetSufficiencyStats();
}

ModelWrapper::ModelWrapper(Bidirectional *v_bidir) {
   type = BIDIR;
   gene  = NULL;
   bidir = v_bidir;
   bidir->resetSufficiency();
}

ModelWrapper::~ModelWrapper() {
   if (gene != NULL) delete gene;
   if (bidir != NULL) delete bidir;
}

double ModelWrapper::getMu() {
  if (type == FMOD) {
     return gene->bidir.getMu();
  } else if (type == BIDIR) {
     return bidir->getMu();
  } else {  // Undefined type
     return 0.;
  }
}

std::string ModelWrapper::write_out() {
   std::string output;
   if (type == FMOD) {
     output = "Full Model:\n";
     output += gene->write_out();
   } else if (type == BIDIR) {
     output = "Bidirectional :\n";
     output += bidir->write_out();
   } else {
     output = "No Model Here!\n";
   }
   return output;
}

void ModelWrapper::resetSufficiencyStats() {
   if (type == FMOD) {
      gene->resetSufficiencyStats();
   } else if (type == BIDIR) {
      bidir->resetSufficiency();
   } else { // Model undefined.
     // Should never be here!
   }
}

double ModelWrapper::getResponsibility() {
   if (type == FMOD) {
      return gene->getResponsibility();
   } else if (type == BIDIR) {
      return bidir->getResponsibility();
   } else {
     // Should never be here!
     return 0.;
   }
}

// component::set_priors
void ModelWrapper::setPriors() {
	/*
	//============================
	//for sigma
	alpha_0 	= 20.46;
	beta_0 		= 10.6;
	//============================
	//for lambda
	alpha_1 	= 20.823;
	beta_1 		= 0.5;
	//==============================
	//for initial length of Uniforms
	alpha_2 	= 1.297;
	beta_2 		= 8260;
	*/

	//*****************************************************
	/* bidir.w 	= (r + ALPHA_2) / (N + ALPHA_2*K*3 + K*3) ; */
		
	//w_thresh= ( ALPHA_2 ) / (N + ALPHA_2*K*3 + K*3 );
}

// component::initialize_bounds
void ModelWrapper::initializeBounds(double mu_seed) {
   /*
   footprint = fp;
   EXIT = false;

   int complexity = 1;
   if (data->strand == ".") {
      complexity = 3;
   } else { 
      complexity = 2;
   }

   //====================================
   double sigma, lambda, pi_EMG, w_EMG;
   double b_forward, w_forward;
   double a_reverse, w_reverse;

   //====================================
   // start sampling
   // for the bidirectional/EMG component
   Random ran_num_gen;

   sigma = ran_num_gen.fetchUniform(1, 250) / scale;
   lambda = scale / ran_num_gen.fetchUniform(1, 250);
   double dist = (1.0 / lambda);
   int j = get_nearest_position(data, mu, dist);
   int k = get_nearest_position(data, mu, forward_bound - mu);
   double pi = 0.5;
   if (data->strand == "+")
   {
      pi = 1.0;
   }
   else if (data->strand == "-")
   {
      pi = 0.;
   }
   forward = UNI(mu + (1.0 / lambda), data->maxX, 1.0 / (complexity * K), 1, j, pi);
   forward.j = j, forward.k = k;

   bidir = EMG(mu, sigma, lambda, 1.0 / (complexity * K), 0.5);
   bidir.foot_print = fp;

   dist = -(1.0 / lambda);
   j = get_nearest_position(data, mu, dist);
   k = get_nearest_position(data, mu, reverse_bound - mu);
   reverse = UNI(data->minX, mu - (1.0 / lambda), 1.0 / (complexity * K), -1, j, 1 - pi);
   reverse.j = k, reverse.k = j;
   termination = (termination > 1);

   type = 1;
   */
}

/**
 * @brief 
 * 
 * @param N 
 * @param K 
 */
void ModelWrapper::updateParameters(double N, int K) {
   if (type == FMOD) {
      gene->updateParameters(N,K);
   } else if (type == BIDIR) {
      bidir->updateParameters(N,K);
   } else {
     // Should never be here!
   }
}

perStrandInfo ModelWrapper::calculateRi(double z, perStrandInfo coverage) {
   if (type == FMOD) {
      return gene->calculateRi(z,coverage);
   } else if (type == BIDIR) {
      return bidir->calculateRi(z,coverage);
   } else {
      perStrandInfo empty;
      return empty;
   }
}

void ModelWrapper::updateExpectations(double z, perStrandInfo coverage, perStrandInfo normalizeRi) {
   if (type == FMOD) {
      gene->updateExpectations(z,coverage, normalizeRi);
   } else if (type == BIDIR) {
      bidir->updateExpectations(z,coverage,normalizeRi);
   } else {
      // Not a valid model!
   }

}

/***************** sets of models *************/

ModelContainer::ModelContainer()
	: noise() {
  K = 0;
  ll = nINF;
}

ModelContainer::ModelContainer(int v_k, ModTypes type) {
  K = v_k;
  ModelWrapper *aModel;
  for (int i = 0; i < K; i++) {
     if (type == FMOD) {   // Full Model
      FullModel *thisModel;
      thisModel = new FullModel();
      aModel = new ModelWrapper(thisModel);
     } else {  // Bidirectionals only
      Bidirectional *thisModel;
      thisModel = new Bidirectional();
      aModel = new ModelWrapper(thisModel);
     }
     setModels.push_back(aModel);
  }
  ll = nINF;
}

ModelContainer::~ModelContainer() {
  for (int i = 0; i < K; i++) {
     delete setModels[i];
  } 
}

std::string ModelContainer::write_out() {
   std::string output;
   output = "Noise: " + noise.write_out();
   output += "\n" + std::to_string(K) + " models:\n";
   for(int i = 0; i < K; i++) {
     output += setModels[i]->write_out();
     output += "\n";
   }
   output += "ll: " + tfit::prettyDecimal(ll,4);
   return output; 
}

/**
 * @brief Pre-EM algorithm setup.  Includes:
 * (1) initialize models with hyperparameters
 * (2) Seed the positions of mu
 * (3) Based on seeds, set the bounds of each model, including Noise
 */
void ModelContainer::initializeComponents(dInterval *data) {
	//initialize(1) components with user defined hyperparameters
	for (int k = 0; k < K; k++){
      setModels[k]->setPriors();
		//components[k].set_priors(ALPHA_0, BETA_0, ALPHA_1, BETA_1, ALPHA_2, ALPHA_3,data->N, K);
	}
	//random seeding, set bounds 
   Seed_SetBounds();    
   SortByMu(); // Joey's: sort_components(components, K);
   //  *** Joey's code resets the links for forward_neightbor and reverse_neighbor here.

   // Initialize the Noise model
   // Needs a pi (from classifier in Joey) and weight, noise_w in Joey
   // noise.setBounds(data->minX, data->maxX); 
}

void ModelContainer::Seed_SetBounds() {
	std::vector<double> mus;
   mus = RandomSeeds();	// randomly seed
   std::sort(mus.begin(), mus.end()); // sort random seed; replaces Joey's sort_vector
   for (int k = 0; k < K; k++) { 
      // Joey: initialize_bounds();
      setModels[k]->initializeBounds(mus[k]);   // use random seeds to initialize
      /** (mus[k],
                                      data, K, data->SCALE, 0., topology, foot_print, data->maxX, data->maxX);
                                      */
   }
}

/**
 * @brief Fetch a random seed.
 * Joey's code passes around a set of random seeds (mu_seeds) obtained by
 * template matching.  It then samples from this set (removing as used). 
 * There is also r_mu which is set at 0.05 by parameters.   It uses this
 * as the variance to a normal centered at the midpoint of the region.  It
 * then samples from THIS Normal for a random mu.
 */
std::vector<double> ModelContainer::RandomSeeds() {
  std::vector<double> mus; mus.reserve(K);
  double mu = 0;  int tempIndex = 0;
  Random ran_num_generator;
  for (int k = 0; k < K; k++){ // Each model
  /*
     if (mu_seeds.size() > 0) {  // mu_seeds are a set of guesses at mu
        tempIndex = sample_centers(mu_seeds, p);
        mu = mu_seeds[i];
        if (r_mu > 0) {
           mu = ran_num_generator.fetchNormal(mu, r_mu);
        }
     } else {
        mu = ran_num_generator.fetchNormal((data->minX + data->maxX) / 2., r_mu);
     }

     mus[k] = mu;
     if (mu_seeds.size() > 0) {
        mu_seeds.erase(mu_seeds.begin() + tempIndex);
     }
     */
  }
  return mus;
}

void ModelContainer::SortByMu() {
   std::sort(setModels.begin(), setModels.end(), tfit::compareMu); 
}

void ModelContainer::resetAllSufficiencyStats() {
   for (int k = 0; k < K; k++) {
      setModels[k]->resetSufficiencyStats();
   }
   noise.resetSufficiency();
}

double ModelContainer::getAllResponsibilities() {
   double N = 0;
	for (int k = 0; k < K; k++){
	   N += setModels[k]->getResponsibility();
	}
   // N +=  noise.getResponsibility(); // Or should this just be zero?!?
   return N;
}

/**
 * @brief 
 * 
 * @param z       data position coordinates
 * @param coverage   Coverage per strand at this position
 * @return perStrandInfo 
 */
perStrandInfo ModelContainer::calculateAllRi(double z, perStrandInfo coverage) {
   perStrandInfo perModelRi;
   perStrandInfo normalizeRi;
   // Equation 7 in Azofeifa 2017: calculate r_i^k
   for (int k = 0; k < K; k++) { // computing the responsibility terms per model
      perModelRi = setModels[k]->calculateRi(z, coverage);
      normalizeRi.forward += perModelRi.forward;
      normalizeRi.reverse += perModelRi.reverse;

   }
   // plus Noise!
   perStrandInfo noiseRi = noise.calculateRi(z, coverage); 
   normalizeRi.forward += noiseRi.forward;
   normalizeRi.reverse += noiseRi.reverse;

   return normalizeRi;
}

void ModelContainer::updateExpectations(double z, perStrandInfo coverage, perStrandInfo normalizeRi) {
   // now we need to add the sufficient statistics 
   for (int k = 0; k < K ; k++) {
      setModels[k]->updateExpectations(z, coverage, normalizeRi);
   }
   // plus noise!
   noise.updateExpectations(coverage, normalizeRi);
}

