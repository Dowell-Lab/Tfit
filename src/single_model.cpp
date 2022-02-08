/**
 * @file single_model.cpp
 * @author Robin Dowell 
 * @brief Routines necessary for loading data and segments.
 * @version 0.1
 * @date 2022-02-08
 * 
 */
#include "single_model.h"

#include <math.h>   
#include <stdio.h>
#include <unistd.h>

#include <string>
#include <vector>

#include "split.h"

parameters::parameters() {
  // Empty constructor.
  mu = sigma = lambda = pi = footprint = omega = 0.0;    
}

/**
 * @brief Construct a new parameters::parameters object
 * 
 * @param mu        The center of the bidirectional, position of loading
 * @param sigma     The variance on the loading position (mu)
 * @param lambda    The length of the exponential.  EMG = N(mu,sigma) + Exp(lambda)
 * @param pi        The strand bias
 * @param footprint     The offset parameter (aka Beta)
 * @param omega     The pausing probability. 
 */
parameters::parameters(double mu, double sigma, double lambda, double pi, double footprint, double omega) {
    mu = mu;
    sigma = sigma;
    lambda = lambda;
    pi = pi;
    footprint = footprint;
    omega = omega;
}
/**
 * @brief Convert the parameters into an ouput string.
 * 
 * @return std::string 
 */
std::string parameters::write() {
  std::string output =    std::to_string(mu) +"\t" + std::to_string(sigma) +"\t" + 
    std::to_string(lambda) + "\t" + std::to_string(pi) +"\t" + std::to_string(footprint) 
    + "\t" + std::to_string(omega) + "\n";
   return output;
}
/**
 * @brief This is the opposite of the write() function.  Expects
 * the write() function's output as its input.
 * 
 * @param single  Tab delimited list of parameters, 
 * in order: mu, sigma, lambda, pi, footprint, and omega.
 */
void parameters::read(std::string single) {
   vector<std::string> split_tab = split_by_tab(single, ""); // note delim variable ignored!

   mu = std::stod(split_tab[0]);
   sigma = std::stod(split_tab[1]);
   lambda = std::stod(split_tab[2]);
   pi = std::stod(split_tab[3]);
   footprint = std::stod(split_tab[4]);
   omega = std::stod(split_tab[5]);

   // Should we check that we don't have more parameters?
}

/**
 * @brief The start coordinate for a bed file with this call.
 * 
 * @return std::string 
 */
double parameters::getStart() {
    return max(mu-(sigma+lambda), 0.0);
}

/**
 * @brief The end coordinate for a bed file with this call.
 * 
 * @return std::string 
 */
double parameters::getEnd() {
   return (mu + (sigma+lambda));
}

/******************************************************/
/* Class: Set_parameters */
/*************************/

Set_parameters::Set_parameters() {
   K = 0;
   // Should we do something with the collection (set to NULL)?
}

Set_parameters::Set_parameters(int K) {
    K = K;
   // Should we do something with the collection (set to NULL)?
}
