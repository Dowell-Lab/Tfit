/**
 * @file helper.h
 * @author Robin Dowell 
 * @brief A series of helper functions (random numbers, etc)
 * @version 0.1
 * @date 2022-02-28
 * 
 */
#ifndef helper_H 
#define helper_H 

#include <string>
#include <vector>
#include <random>

class Random {
	std::mt19937 mt;

    bool testing;
    int  tseed;

    // Constructors
    Random();
    Random(int);

    // Wrapper functions
    double FetchUniform(double,double);
    double FetchNormal(double,double);  // return a random number in this interval, normal dist.
    double FetchProbability();
};

#endif
