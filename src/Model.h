/**
 * @file Model.h
 * @author Robin Dowell 
 * @version 0.1
 * @date 2022-05-12
 * 
 */
#ifndef Model_H
#define Model_H

#include <string>
#include <vector>

#include "Distro.h"

class Bidirectional {
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

  double pdf(double z, char s);
  double ExpX(double z, char strand);
  double ExpY(double z, char s);
  double ExpX2(double z, char strand);
  double ExpY2(double z, char s);

  std::vector<double> generate_data(int n);

// private:
  double normalPDF(double);
  double millsRatio(double);
  double normalCDF(double x);
  int indicatorStrand(char s);
  double applyFootprint (double z, char s);

};

class NoiseModel {
  Uniform noise;

  // Constructor
  NoiseModel();
  NoiseModel(double v_a, double v_b);

  // Expectation: pG(b-a)/S where S is length of genome, G is # reads mapped, p is prob noise

  // Functions
  std::string write_out();
  double pdf(double x, char s);

};

class FullModel {
  public:
  Bidirectional bidir;
  Uniform forwardElongation;
  Uniform reverseElongation;
  double w_forward, w_reverse;

  // Constructor
  FullModel();

  // Functions
  std::string write_out();
};


#endif
