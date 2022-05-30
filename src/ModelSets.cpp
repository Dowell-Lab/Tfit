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
#include "helper.h"

/**********************  Model Wrapper (for EM) *************/

ModelWrapper::ModelWrapper() {
   type = EMPTY;         
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
   if (gene != NULL) {
      gene->resetSufficiencyStats();
   } else if (bidir != NULL) {
      bidir->resetSufficiency();
   }
}


/***************** sets of models *************/

ModelContainer::ModelContainer()
	: noise() {
  K = 0;
  ll = nINF;
  pi = 0;	
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
  pi = 0;	
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
   output += "ll: " + tfit::prettyDecimal(ll,4) + " pi: " + tfit::prettyDecimal(pi,3);
   return output; 
}


