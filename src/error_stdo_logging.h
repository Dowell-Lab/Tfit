/**
 * @file error_stdo_logging.h
 * @author Joey Azofeifa
 * @brief 
 * @version 0.1
 * @date 2016-05-20
 * 
 */
#ifndef error_stdo_logging_H
#define error_stdo_logging_H

#include <fstream>
#include <iostream>
#include <string>

using namespace std;

/**
 * @brief Presumably holds file info for logging?
 * 
 */
class Log_File{
public:
	int job_ID;
	string job_name;
	int rank; 
	string log_out_dir;
	ofstream FHW;
	Log_File();
	Log_File(int, int, string, string);
	void write(string, int);
};


#endif
