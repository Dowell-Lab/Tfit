/**
 * @file split.h
 * @author Joey Azofeifa
 * @brief 
 * @version 0.1
 * @date 2016-05-20
 * 
 */
#ifndef split_H
#define split_H

#include <string>
#include <vector>

using namespace std;

vector<string> string_split(string , const char);
vector<string> splitter(string ELE, string D);
string strip(string ELE, string D);
string join(vector<string>, string);
vector<string> splitter2(string ELE, string D);
vector<string> split_by_bar(string ELE, string D);
vector<string> split_by_colon(string ELE, string D);
vector<string> split_by_tab(string ELE);
vector<string> split_by_comma(string ELE, string D);
vector<string> split_by_dash(string ELE, string D);

#endif
