/**
 * @file Distros.cpp
 * @author Robin Dowell 
 * @brief Contains the basic mathematical distribution functions for the
 * Normal, Exponential, EMG, and Uniform
 * @version 0.1
 * @date 2022-05-22
 * 
 */
#include "Distros.h"

Uniform::Uniform() {
   a = b = 0;	
}

std::string Uniform::write_out() {
	std::string output = "U(" + std::to_string(a) + "," + to_string(b) + ")";
    return output;
}

EMG2::EMG2() {
  mu = sigma = lambda = 0.0;
}

std::string EMG2::write_out() {
	std::string output = "EMG(" + std::to_string(mu) + "," + to_string(sigma) 
            + "," + std::to_string(lambda) + ")";
    return output;
}


