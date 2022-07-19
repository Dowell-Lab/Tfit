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
#include "Data.h"   // gInterval

enum ModTypes {EMPTY, FMOD, BIDIR};
extern double nINF;     // currently defined in template_matching?

/**
 * @brief This class provides a unified front for the two different
 * models we can use within the EM (Bidirectional vs FullModel).
 */
class ModelWrapper {
  public:
  ModTypes type;  

  // Only one of the below will be non-NULL at any one time.
  // Where should this be checked?
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
  double calculateRi(double z, char strand);
  void updateExpectations(double i, perStrandInfo coverage, perStrandInfo normalizeRi);
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

  int K;  //<! number of models in the set
  std::vector<ModelWrapper *> setModels;
  double ll; // log likelihood 

  bidirConstraints limits;    // constraints on the bidir models

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

  void resetAllSufficiencyStats();
  double getAllResponsibilities();
  perStrandInfo calculateAllRi(double i, dInterval *data);
  void updateExpectations(double i, dInterval *data, perStrandInfo normalizeRi);

};

#endif
