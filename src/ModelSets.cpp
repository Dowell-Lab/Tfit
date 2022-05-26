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

   resetSufficiencyStats();
   weight = 0;
}

ModelWrapper::ModelWrapper(FullModel *v_gene) {
   type = FMOD;
   gene  = v_gene;
   bidir = NULL;
}

ModelWrapper::ModelWrapper(Bidirectional *v_bidir) {
   type = BIDIR;
   gene  = NULL;
   bidir = v_bidir;
}

void ModelWrapper::resetSufficiencyStats() {
   if (type == FMOD) {  // FullModel
      gene->rTerms.reset();
   } else if (type == BIDIR) {
      bidir->rTerms.reset();
   }
}

/***************** sets of models *************/

ModelContainer::ModelContainer()
	: noise() {
  K = 0;
  setModels = NULL;
  ll = nINF;
  pi = 0;	
}

ModelContainer::ModelContainer(int v_k, ModTypes type) {
  K = v_k;
  for (int i = 0; i < K; i++) {
     // Allocate model of correct type
  }
  ll = nINF;
  pi = 0;	
}


