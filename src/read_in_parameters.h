/**
 * @file read_in_parameters.h
 * @author Joey Azofeifa
 * @brief  Parameters are kept in a large map/hash. 
 * @version 0.1
 * @date 2016-05-20
 */
#ifndef read_in_parameters_H
#define read_in_parameters_H

#include <unistd.h>

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

/**
 * @brief The parameters class.
 * Keeps a singular map ("p") which is a hash of all parameters.
 */
class params{
public:
	map<string, string> p;
	int N;
	string module;
	bool bidir;
	bool model;
	bool CONFIG;
	bool select;

	map<string, string> p2;
	map<string, string> p3;
	map<string, string> p4;
	map<string, string> p5;
	map<string, string> p6;

	bool EXIT;
	
	// Should these be private?
	const char * isIntGroup[8] = {"-pad", "-minK", "-maxK", "-rounds", 
	"-mi", "-MLE", "-elon", "-merge"};

	const char * isDecGroup[17]  = {  "-br","-ns", "-ct", "-max_noise", 
	"-r_mu", "-ALPHA_0", "-ALPHA_1", "-ALPHA_2", "-BETA_0", "-BETA_1", 
	"-bct", "-ms_pen" , "-lambda", "-sigma", "-foot_print", "-pi", "-w" }; 

	const char * isPathGroup[8] = {"-config", "-i", "-j", "-k", "-tss", 
"-log_out", "-o", "-q"};
	
	//Constructors
	params();

	//Functions
	void display(int,int);
	void help();
	string get_header(int);
	vector<string> validate_parameters();
};

/* Deprecated: These don't appear to have code anywhere 
void fillInOptions(char*,params);
params * readInParameters(char**);
void fill_in_bidir_boostrap(params *);
*/
const std::string currentDateTime();
int read_in_parameters( char**, params *, int );

#endif
