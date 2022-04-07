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
#include <map>
#include <random>

/**
 * @brief Allows for all spots of random calls to become
 * predictable when necessary for testing.
 * 
 */
class Random {
    public:
	std::mt19937 mt;

    // Constructors
    Random();

    // Wrapper functions
    double fetchUniform(double,double);
    double fetchNormal(double,double);  // return a random number in this interval, normal dist.
    double fetchProbability();
};

/**
 * @brief Keeps a mapping function that is bidirectional,
 * i.e. map<name,index> and map<index,name>
 * Will be used in converting strings to indexes.
 * 
 */
class Bimap {
  public:
   int num_elements;     // Total number of elements

   std::map<int,std::string> index2str;
   std::map<std::string, int> str2index;   

   // Constructors
   Bimap(); 

   // Functions:
   void addIdentifier(std::string);    // Add new name to indexes
   int lookupIndex(std::string);    // Given a name, give correct index
   std::string lookupName(int);  // Given an index, what is the corresponding name?

   std::string print_index_names();  

};

#endif
