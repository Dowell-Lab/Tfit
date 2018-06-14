#include "split.h"
#include <iostream>
#include <sstream>
using namespace std;

/** Splits a given string by the provided delimiter.
 * @param s String to split
 * @param delimiter Character by which to split the string.
 * @return Set of tokens within the string.
 */
vector<string> string_split(string s, const char delimiter)
{
  size_t start=0;
  size_t end=s.find_first_of(delimiter);
  
  vector<string> output;
  
  while (end <= string::npos)
    {
      output.emplace_back(s.substr(start, end-start));
      
      if (end == string::npos)
	break;
      
      start=end+1;
      end = s.find_first_of(delimiter, start);
    }
  
  return output;
}

/** Splits a string by another string.
 * @param ELE The string to split.
 * @param D The delimiter string.
 * @return Set of tokens within the string.
 */
vector<string> splitter(string ELE, string D){
	int j = 0;
	vector<string> results;
	const char *d =D.c_str();
	while (not ELE.empty() and j != ELE.size()){
		if (ELE[j] == *d){
			results.push_back(ELE.substr(0,j));
			ELE=ELE.substr(j+1,ELE.size());
			j=0;
		}
		j++;
	}
	results.push_back(ELE.substr(0,j));

	return results;
}

/** Splits a string by another string. Evidently, a string must have a tab within it.
 * @param line The line to be tokenized.
 * @param delim Set of delimiters to use in parsing the line.
 * @return Set of tokens.
 */
vector<string> splitter2(string line, string delim){
	vector<string> tokens;
	istringstream iss(line);
	string token;
	while(std::getline(iss, token, '\t' )){   // but we can specify a different one
		tokens.push_back(token);
	}
	return tokens;
}

/** Splits a string by the vertical bar character |
 * @param line Input string
 * @param delim Additinoal delimiters to use (unused).
 * @return Set of tokens.
 */
vector<string> split_by_bar(string line, string delim){
	vector<string> tokens;
	istringstream iss(line);
	string token;
	while(std::getline(iss, token, '|' )){   // but we can specify a different one
		tokens.push_back(token);
	}
	return tokens;
}

/** Splits a string by a colon character.
 * @param line string
 * @param delim Additional delimiters to use (unused)
 * @return Set of tokens. 
 */
vector<string> split_by_colon(string line, string delim){
	vector<string> tokens;
	istringstream iss(line);
	string token;
	while(std::getline(iss, token, ':' )){   // but we can specify a different one
		tokens.push_back(token);
	}
	return tokens;
}

/** Splits a string by a tab character.
 * @param line string
 * @param delim Additional delimiters to use (unused)
 * @return Set of tokens. 
 */
vector<string> split_by_tab(string line, string delim){
	vector<string> tokens;
	istringstream iss(line);
	string token;
	while(std::getline(iss, token, '\t' )){   // but we can specify a different one
		tokens.push_back(token);
	}
	return tokens;
}

/** Splits a string by a comma character.
 * @param line string
 * @param delim Additional delimiters to use (unused)
 * @return Set of tokens. 
 */
vector<string> split_by_comma(string line, string delim){
	vector<string> tokens;
	istringstream iss(line);
	string token;
	while(std::getline(iss, token, ',' )){   // but we can specify a different one
		tokens.push_back(token);
	}
	return tokens;
}

/** Splits a string by a hyphen character.
 * @param line string
 * @param delim Additional delimiters to use (unused)
 * @return Set of tokens. 
 */
vector<string> split_by_dash(string line, string delim){
	vector<string> tokens;
	istringstream iss(line);
	string token;
	while(std::getline(iss, token, '-' )){   // but we can specify a different one
		tokens.push_back(token);
	}
	return tokens;
}

/** Removes all instances of given delimiters from the string.
 * @param ELE element to strip.
 * @param D set of delimiters to remove.
 * @return Stripped string.
 */
string strip(string ELE, string D){
	const char *d 	= D.c_str();
	string result 	= "";
	for (int i = 0; i < ELE.size(); i++){
		if (ELE[i]==*d){
			break;
		}else{
			result+=ELE[i];
		}
	}
	return result;
}

/** Concatenates a set of strings with the given delimiter.
 * @param toBeJoined Set of strings to concatenate.
 * @param delim Delimiter string to be placed between elements of the previously specified set.
 * @return Concatenated string.
 */
string join(vector<string> toBeJoined, string delim){
	typedef vector<string>::iterator vs_it;
	string result="";
	for (vs_it i =toBeJoined.begin(); i != toBeJoined.end(); i++){
		result=result + delim + *i;
	}
	return result;
}	
