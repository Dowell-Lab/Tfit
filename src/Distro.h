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

/**
 * @brief A normal distribution
 * 
 */
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

/**
 * @brief An exponential distribution
 * 
 */
class Exponential {
  public:
	double lambda;

  // Constructor
  Exponential();
  Exponential(double);

  // Functions
  std::string write_out();

  std::vector<double> generate_data(int n);

  double pdf (double x);
  double cdf (double x);

  double ExpX();   // Expected value = 1/lambda
};

/**
 * @brief A Uniform distribution
 * 
 */
class Uniform {
public:
	double lower, upper;

  // Constructor
  Uniform();
  Uniform(double v_a, double v_b);

  // Functions
  std::string write_out();
  std::vector<double> generate_data(int n);

  double pdf(double x);
  double cdf(double x);

};

#endif
