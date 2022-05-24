/**
 * @file Distro.h
 * @author Robin Dowell 
 * @version 0.1
 * @date 2022-05-24
 * 
 */
#ifndef Distro_H
#define Distro_H

#include <string>
#include <vector>

class Normal {
  public:
	double mu, sigma;

  // Constructor
  Normal();     // Note that default is the standard Normal(0,1)
  Normal(double, double);

  // Functions
  std::string write_out();

  std::vector<double> generate_data(int n);
  double pdf(double x);
  double cdf(double x);

};

class Exponential {
  public:
	double lambda;

  // Constructor
  Exponential();
  Exponential(double);

  // Functions
  std::string write_out();

  std::vector<double> generate_data(int n);
};

class Uniform {
public:
	double lower, upper;

  // Constructor
  Uniform();
  Uniform(double v_a, double v_b);

  // Functions
  std::string write_out();
  std::vector<double> generate_data(int n);

};

#endif
