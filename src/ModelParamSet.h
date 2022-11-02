/**
 * @file ModelParamSet.h
 * @author Robin Dowell
 * @brief Parameters associated with the model/program.
 * @version 0.1
 */
#ifndef ModelParamSet_H 
#define ModelParamSet_H 

#include <map>
#include <string>
#include <vector>

#include "Models.h"

/**
 * @brief This replaces the parameters string in segment_fits.
 * 
 * Current output is: 
 * mus+"\t"+sigmas+"\t"+lambdas+"\t" + pis+"\t" + fps+ "\t" + ws + "\t" + fbs+"\t" +ras;
 * 
 * @remark Currently keeping all parameters as doubles, as Joey did.  Reality is that
 * they are not really all doubles.  Should be checking those that are probabilities for 
 * proper bounds.
 *
 * This could also be reused as PRIORs.  In that scenario we'd need to be able to 
 * read/write them (something more readable than Joey's current K-models format,
 * perhaps JSON?!?!).
 * We'd need a way to specify a parameter as "not specified" in the current priors 
 * (i.e. to permit subsets of priors only).
 * 
 */
class ModelParams {
public:
  double mu;    // Technically an int
  double sigma;
  double lambda;
  double footprint;
  double pi;  // This one is a probability
  double omega[3]; // weights (this mimics Joey's but is ugly)

  double negL;  // reverse strand extent of elongation
  double posL;  // forward strand extent of elongation

  // Constructors & Destructors
  ModelParams();
  ModelParams(double,double,double,double,double,double);

    // Getters and Setters
  void SetfromBidirectional(Bidirectional model);
  void SetfromFullModel(FullModel model);
  double getBedStart();
  double getBedEnd();

  // Functions
  std::string write();        
  void read(std::string);
  std::vector<std::string> fetch_as_strings();  // Needed for writeAsKModels?
  std::string writeAsJSON();  // Output in JSON format
  void readFromJSON(std::string);


};

/**
 * @brief A collection of ModelParams, i.e. K models.
 *  Should this provide iterator behavior? 
 */
class ModelParamSet {
  public:
  std::vector<ModelParams *> collection;  // Contains K models
  double log_likelihood;

  //Constructors
  ModelParamSet();
  ModelParamSet(int);
  ~ModelParamSet();

  // Functions
  std::string write();
  void read_from_K_models(std::string);  // This will parse/read in a K_models string
  std::string writeAsKmodels();
  void readFromKmodels();
  std::string writeJSON();
  void readJSON();
  
  private:
  int K;    // Number of models
};


#endif