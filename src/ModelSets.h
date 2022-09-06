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
#include "Data.h"   // bed12 

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
  void initializeBounds(double v_mu, double v_sigma, 
     double v_lambda, double v_minX, double v_maxX);
  void resetSufficiencyStats();
  double getResponsibility();

  double getMu();

  void updateParameters(double N, int K);
  perStrandInfo calculateRi(double z, perStrandInfo coverage);
  void updateExpectations(double i, perStrandInfo coverage, perStrandInfo normalizeRi);

  void initializeBounds(double v_mu, double v_sigma, 
            double v_lambda, double v_weight, double minX, double maxX);
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
  void initializeWithPriors(dInterval *data);
  // void useSeeds2SetBounds();
  void SortByMu();

  void resetAllSufficiencyStats();
  double getAllResponsibilities();
  perStrandInfo calculateAllRi(double i, perStrandInfo coverage);
  void updateExpectations(double z, perStrandInfo coverage, perStrandInfo normalizeRi);
};

#endif
