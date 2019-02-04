//Include guard
#ifndef SEQUENCE_H
#define SEQUENCE_H


//Forward declared dependencies

//Included dependencies
#include <vector>
#include <string>
#include <gsl/gsl_rng.h>

class Sequence{
private:
    //Internal variables
    std::vector<char> seq;
    
public:
    Sequence();
    Sequence(int length);
    Sequence(std::vector<char>& s);
    Sequence(std::string s);
    ~Sequence();
    
    void setSeq(std::vector<char>& s);
	std::string getSeqAsString();
    void setBase(int index, char b);
    void addBase(char b);
	void removeBase(int index);
    char getBase(int index);
    int getLength();
    void print();
    std::string printToString();
    std::string printToStringNoHyphens();
    Sequence subset(int start, int finish);
    void addBaseFront(char b);
};

#endif 
