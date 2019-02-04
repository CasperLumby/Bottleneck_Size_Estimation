//Code taken from iain as posted on http://stackoverflow.com/questions/865668/how-to-parse-command-line-arguments-in-c

//Included dependencies
#include "inputParser.hpp"
#include <iostream>
#include <algorithm>

using namespace std;

//Constructors
InputParser::InputParser(int &argc, char **argv){
	for (int i=1; i < argc; ++i) {
		this->tokens.push_back(string(argv[i]));
	}
}


/// @author iain
const string InputParser::getCmdOption(const string &option) const {
	vector<string>::const_iterator itr;
	itr =  find(this->tokens.begin(), this->tokens.end(), option);
	if (itr != this->tokens.end() && ++itr != this->tokens.end()){
		return *itr;
	}
		return "";
}
/// @author iain
bool InputParser::cmdOptionExists(const string &option) const {
	return find(this->tokens.begin(), this->tokens.end(), option)
		!= this->tokens.end();
}


