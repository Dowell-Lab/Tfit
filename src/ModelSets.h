/**
 * @file ModelSets.h
 * @author Robin Dowell 
 * @version 0.1
 * @date 2022-05-26
 * 
 */
#ifndef ModelSets_h 
#define ModelSets_h 

#include <string>
#include <vector>

#include "Models.h"

enum ModTypes {EMPTY, FMOD, BIDIR};
extern double nINF;     // currently defined in template_matching?

/**
 * @brief This class provides a unified front for the two different
 * models we can use within the EM (Bidirectional vs FullModel).
 * For the EM, we'll want to keep some intermediate information
 * during the M/E steps.  This class is designed to encompass those 
 * extra quantities.
 * 
 */
class ModelWrapper {
  public:
  ModTypes type;  

  // Only one of the below will be non-NULL at any one time.
  FullModel *gene;
  Bidirectional *bidir;

  // Constructor
  ModelWrapper();
  ModelWrapper(FullModel *);
  ModelWrapper(Bidirectional *);
  // Destructor
  ~ModelWrapper();

  // Functions
  std::string write_out();

  void setPriors();
  void initializeBounds(double mu_seed);
  void resetSufficiencyStats();
  double getResponsibility();

  double getMu();

  void updateParameters(double N, int K);
};

/**
 * @brief The container holds K total models.   This is what
 * the EM will be using to run on and returning from each run.
 * 
 * Note that somewhere else we've got to decide what K is.
 * 
 */
class ModelContainer {		
  public:
  UniformModel noise;

  int K;
  std::vector<ModelWrapper *> setModels;
  double ll,pi; // log likelihood, 

  // Constructor
  ModelContainer();
  ModelContainer(int, ModTypes);
  //Destructor
  ~ModelContainer();

  // Functions;
  std::string write_out();
  void initializeComponents(dInterval *data);
  void Seed_SetBounds();

  std::vector<double> RandomSeeds();
  void SortByMu();

  void resetSufficiencyStats();
  double getAllResponsibilities();
};

#endif