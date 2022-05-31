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


class BasicModel {
  public:
  double weight;
  Responsibilities sufficiencyStats;   // The current read

  //Constructor
  BasicModel();

  std::string write_out();
  double pdf(double z, char s);   // pdf on strand s for position z

  void resetSufficiency();
};

class Bidirectional: public BasicModel {
  public:
  Normal loading;
  Exponential initiation;
  double pi;		// strand bias
  double footprint;		// the ad hoc footprint parameter 

  // Constructor
  Bidirectional();
  Bidirectional(double, double, double, double, double);

  // Functions
  std::string write_out();

  double pdf(double z, char s);   // standard (Grushka 1972) formula for pdf
  double pdf_alt(double z, char s);   // Kalambet 2011 pdf for numerical stability

  double ExpX(double z, char strand);
  double ExpY(double z, char s);
  double ExpX2(double z, char strand);
  double ExpY2(double z, char s);

  double getMu();

  std::vector<double> generate_data(int n);

// private:
  double millsRatio(double);    
  int indicatorStrand(char s);
  double applyFootprint (double z, char s);
};

class UniformModel: public BasicModel {
  public:
  Uniform uni;
  double pi;    // Why does this need strand bias?

  // Constructor
  UniformModel();
  UniformModel(double v_a, double v_b);

  // Functions
  std::string write_out();
  double pdf(double x, char s);
  double getResponsibility();   

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

  void resetSufficiencyStats();
  double getResponsibility();   
};


#endif
