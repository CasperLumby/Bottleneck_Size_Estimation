//Code taken from iain as posted on http://stackoverflow.com/questions/865668/how-to-parse-command-line-arguments-in-c

//Include guard
#ifndef INPUTPARSER_HPP
#define INPUTPARSER_HPP


//Included dependencies
#include <string>
#include <vector>

class InputParser{

public:
	InputParser (int &argc, char **argv);
        const std::string getCmdOption(const std::string &option) const;
        bool cmdOptionExists(const std::string &option) const;
    
private:
        std::vector <std::string> tokens;
};

#endif 
