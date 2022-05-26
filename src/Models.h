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


/**
 * @brief Every model component must be able to report on its 
 * responsiblities, e.g. Prob(k|data point, params) in the EM 
 * algorithm. 
 * 
 */
class Responsibilities {
  public:
  double r_forward, r_reverse; //responsibilities per strand

  //Constructor
  Responsibilities();

  //Functions
  void reset();
};

class HyperParameters {
  //FOR SIGMA ; variance in loading, gamma
	double ALPHA_0, BETA_0;
	//FOR LAMBA ; rate of initiation, gamma
	double ALPHA_1, BETA_1;
	//FOR W ; weight , Dirichlet
	double ALPHA_2;
	//FOR PI ; strand prob. , beta
	double ALPHA_3;

  //Constructor
  HyperParameters();
};

class Bidirectional {
  public:
  Normal loading;
  Exponential initiation;
  double pi;		// strand bias
  double footprint;		// the ad hoc footprint parameter 

  Responsibilities rTerms;
  Responsibilities currentRterms;

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

  std::vector<double> generate_data(int n);

// private:
  double millsRatio(double);    
  int indicatorStrand(char s);
  double applyFootprint (double z, char s);


};

class FullModel {
  public:
  Bidirectional bidir;
  Uniform forwardElongation;
  Uniform reverseElongation;
  double w_forward, w_reverse;

  Responsibilities rTerms;
  Responsibilities currentRterms;

  // Constructor
  FullModel();

  // Functions
  std::string write_out();
};

class NoiseModel {
  public:
  Uniform noise;
  Responsibilities rTerms;
  Responsibilities currentRterms;

  // Constructor
  NoiseModel();
  NoiseModel(double v_a, double v_b);

  // Expectation: pG(b-a)/S where S is length of genome, G is # reads mapped, p is prob noise

  // Functions
  std::string write_out();
  double pdf(double x, char s);

};

#endif
