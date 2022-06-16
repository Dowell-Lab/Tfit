/**
 * @file Models.h
 * @author Robin Dowell 
 * @version 0.1
 * @date 2022-05-12
 * 
 */
#ifndef Models_h
#define Models_h

#include <string>
#include <vector>

#include "Distro.h"
#include "ModelSupport.h"   
#include "Data.h"   // gInterval


class BasicModel {
  public:
  double weight;    // Necessary for EM
  double pi;  //strand bias
  Responsibilities sufficiencyStats;   // The current read

  // Priors
  Priors betaPi; // prior on pi, beta
  Priors DirichletW; // prior on weight, Dirichlet

  //Constructor
  BasicModel();

  std::string write_out();
  double pdf(double z, char s);   // pdf on strand s for position z

  void setPriorPi(double);
  void setPriorWeight(double);
  void updateParameters(double,double);
  void updateExpectations(perStrandInfo coverage, perStrandInfo normalizedRi);

  double getResponsibility();   
  void resetSufficiency();
  double calculateRi(double z, char strand);

};

class Bidirectional: public BasicModel {
  public:
  Normal loading;
  Exponential initiation;
  double footprint;		// the ad hoc footprint parameter 

  // Priors
  Priors GammaSigma; // prior: variance in loading, gamma
  Priors GammaLambda; // prior: rate of initiation, gamma

  // Convenience Variables for calculating updates and responsibilities
  sumExpected sumOverN; // Sums of assorted variables from i = 1 to N

  // Constructor
  Bidirectional();
  Bidirectional(double, double, double, double, double);

  // Functions
  std::string write_out();

  // Get and Set the EMG parameters
  double getMu();
  double getSigma();
  double getLambda();
  void setMu(double);
  void setSigma(double);
  void setLambda(double);
  void setPriorSigma(double, double);
  void setPriorLambda(double, double);

  // Calculate the pdf and expected values for a single read z with strand s
  double pdf(double z, char s);   // standard (Grushka 1972) formula for pdf
  double pdf_alt(double z, char s);   // Kalambet 2011 pdf for numerical stability
  double ExpX(double z, char strand);
  double ExpY(double z, char s);
  double ExpX2(double z, char strand);
  double ExpY2(double z, char s);

  // Functions of all Models
  std::vector<double> generate_data(int n);
  void updateExpectations(double z, perStrandInfo coverage, perStrandInfo normalizeRi);
  void calcExpectedVals (double position, char strand, double rtimescov);
  void updateParameters(double,double);

  // Supporting Functions (could be private?)
  double millsRatio(double);    
  int indicatorStrand(char s);
  double applyFootprint (double z, char s);


};

/**
 * @brief Adds sufficiency stats, weight, and strand bias
 * to a Uniform distribution.
 *
 * Does the pi value get fixed {0,1} if this is used for
 * elongation?
 */
class UniformModel: public BasicModel {
  public:
  Uniform uni;

  // Constructor
  UniformModel();
  UniformModel(double v_a, double v_b);

  // Functions
  std::string write_out();
  double pdf(double x, char s);

  // Functions of all Models
  void updateParameters(double,double);

  // These are functions specifically for when its NOISE?
  // Noise Expectation: pG(b-a)/S where S is length of genome, 
  //  G is # reads mapped, p is prob noise
  void setPi(double v_pi);
  double calculateLikelihood(dInterval *data);
  void setBounds(double v_lower, double v_upper);

};

class FullModel {
  public:
  Bidirectional bidir;
  UniformModel forwardElongation; // s = '+' 
  UniformModel reverseElongation;  // s = '-' 

  // Constructor
  FullModel();

  // Functions
  std::string write_out();
  double pdf(double z, char s);

  // Functions of all Models
  void resetSufficiencyStats();
  double getResponsibility();   
  void updateParameters(double,double);
  double calculateRi(double z, char strand);

  void updateExpectations(double i, perStrandInfo coverage, perStrandInfo normalizeRi);

};


#endif
