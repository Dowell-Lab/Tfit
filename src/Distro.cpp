/**
 * @file Distro.cpp
 * @author Robin Dowell 
 * @brief Contains the model objects
 * @version 0.1
 * @date 2022-05-24
 * 
 */
#include "Distro.h"

#include <iostream>
#include <cmath>

#include "helper.h"

/******************** Normal Distribution *********************/
Normal::Normal() {
   mu = 0;
   sigma =1;  
}

Normal::Normal(double v_mu, double v_sigma) {
  mu = v_mu;
  sigma = v_sigma; 
}

/**
 * @return std::string distribution parameters as a string
 */
std::string Normal::write_out() {
   return ("N(" + tfit::prettyDecimal(mu,2) + "," + tfit::prettyDecimal(sigma,2) + ")");
}

/**
 * @brief A generating function for creating data that fits this
 * distribution.
 * 
 * @param n     Number of samples desired
 * @return std::vector<double>   n total samples from the distribution
 */
std::vector<double> Normal::generate_data(int n) {
  Random num_generator;
  std::vector<double> output;
  for(int i = 0; i < n; i++) {
    output.push_back(num_generator.fetchNormal(mu,sigma));
  }
  return output;
}

/** 
 * @brief probability density function (pdf) for normal distro
 * Note Normal() allocator defaults to the standard normal.
 * @param x  position on which to evaluate pdf
 * @return double  f(x|mu, sigma^2)
 */
double Normal::pdf(double x) { 
    return 1.0 / (sigma * std::sqrt(2.0 * M_PI)) * exp(-(pow((x-mu)/sigma, 2)/2.0));
}

/**
 * @brief cummulative distribution function (cdf) for normal distro
 * Note Normal() allocator defaults to the standard normal.
 * @param x  position on which to evaluate pdf
 * @return double phi(x)
 */
double Normal::cdf(double x){ 
    return 0.5*(1+std::erf((x-mu)/(sigma*std::sqrt(2))));
}

/******************** Exponential Distribution *********************/

Exponential::Exponential() {
  lambda = 1.0; 
}

Exponential::Exponential(double v_lambda) {
   lambda = v_lambda; 
}

/**
 * @return std::string distribution parameters as a string
 */
std::string Exponential::write_out() {
   return ("E(" + tfit::prettyDecimal(lambda,2) + ")");
}
/**
 * @brief A generating function for creating data that fits this
 * distribution.
 * 
 * @param n     Number of samples desired
 * @return std::vector<double>   n total samples from the distribution
 */

std::vector<double> Exponential::generate_data(int n) {
  Random num_generator;
  std::vector<double> output;
  for(int i = 0; i < n; i++) {
    output.push_back(num_generator.fetchExponential(lambda));
  }
  return output;
}

/** 
 * @brief probability density function (pdf) for exponential distro
 * @param x  position on which to evaluate pdf
 * @return double  f(x|lambda)
 */
double Exponential::pdf (double x) {
   return (lambda * exp(-1.0 *lambda * x));
}

/**
 * @brief cummulative distribution function (cdf) for exponential distro
 * @param x  position on which to evaluate pdf
 * @return double phi(x)
 */
double Exponential::cdf (double x) {
   return (1.0 - exp(-1.0 * lambda * x));
}

/**
 * @brief Expected value Exp[X|lambda]
 * 
 * @return double expected value of X
 */
double Exponential::ExpX() {
   return (1.0 / lambda);  
}

/******************** Uniform Distribution *********************/

Uniform::Uniform() {
  lower = 0;
  upper = 1;
}

Uniform::Uniform(double v_lower, double v_upper) {
  lower = v_lower; 
  upper = v_upper;
}

/**
 * @return std::string distribution parameters as a string
 */
std::string Uniform::write_out() {
   return ("U(" + tfit::prettyDecimal(lower,2) + "," + tfit::prettyDecimal(upper,2) + ")");
}

/**
 * @brief A generating function for creating data that fits this
 * distribution.
 * 
 * @param n     Number of samples desired
 * @return std::vector<double>   n total samples from the distribution
 */
std::vector<double> Uniform::generate_data(int n) {
  Random num_generator;
  std::vector<double> output;
  for(int i = 0; i < n; i++) {
    output.push_back(num_generator.fetchUniform(lower,upper));
  }
  return output;
}

/** 
 * @brief probability density function (pdf) for uniform distro
 * @param x  position on which to evaluate pdf
 * @return double  f(x|lower,upper)
 */
double Uniform::pdf(double x) {
  if ((x > lower) && (x < upper)) {
	   return (1.0 / abs(upper - lower));
  }
  return 0.0;
}

/**
 * @brief cummulative distribution function (cdf) for uniform distro
 * @param x  position on which to evaluate pdf
 * @return double phi(x)
 */
double Uniform::cdf(double x) {
  if (x < lower) {
    return 0.0; 
  } else if (x > upper) {
    return 1.0;
  } else {
    return ((x-lower)/abs(upper-lower));
  }
}
