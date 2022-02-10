/**
 * @file single_model.h
 * @author Robin Dowell
 * @brief Parameters associated with a single instance of the model.
 * @version 0.1
 */
#ifndef single_model_H
#define single_model_H 

#include <map>
#include <string>
#include <vector>

/**
 * @brief This is going to replace the parameters string in segment_fits.
 * This is redundant with other classes in the model, so it likely won't be here long.
 * 
 * Current output is: 
 * mus+"\t"+sigmas+"\t"+lambdas+"\t" + pis+"\t" + fps+ "\t" + ws + "\t" + fbs+"\t" +ras;
 * 
 * Unclear what fbs and ras are .. so leaving those out for now.
 */
class EMGparameters {
public:
  double mu;
  double sigma;
  double lambda;
  double pi;
  double footprint;
  double omega;

  // Constructors
  EMGparameters();
  EMGparameters(double,double,double,double,double,double);

  // Functions
  std::string write();        
  void read(std::string);
  double getStart();
  double getEnd();
  // What about a function that outputs this as BED.  We don't keep
  // currently its chromosome.  We also aren't sure if these will be relative
  // or genomic coordinates yet. 
};

class Set_EMGparameters {
  public:
  std::vector<EMGparameters *> collection;  // Contains K models

   //Constructors
   Set_EMGparameters();
   Set_EMGparameters(int);

   // Functions
  std::string write();
  void read(std::string);
  
  private:
  int K;    // Number of models
};


#endif