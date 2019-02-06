
//Forward declared dependencies

//Included dependencies
#include "sequence.h"
#include "iostream"

using namespace std;

Sequence::Sequence() {}
Sequence::Sequence(int length) : seq(vector<char>(length,'-')){ }
Sequence::Sequence(vector<char>& s) : seq(s) { }
Sequence::Sequence(string s) {
    
    vector<char> SEQ;
    for(unsigned int i=0; i<s.size();i++) {
        SEQ.push_back(s[i]);
    }
    seq=SEQ;
}
Sequence::~Sequence() {}


//Getters and setters
void Sequence::setSeq(vector<char> & s) { seq = s; }
string Sequence::getSeqAsString() {

//	string s = &seq[0]; //Not entirely sure why this works, but it does
	string s(seq.begin(), seq.end());

	return s;
}
void Sequence::setBase(int index, char b) { seq[index]=b; }
void Sequence::addBase(char b) { seq.push_back(b); }
void Sequence::removeBase(int index) { seq.erase(seq.begin()+index); }
char Sequence::getBase(int index) { return seq[index]; }
int Sequence::getLength() { return (int) seq.size(); }
void Sequence::print() {
    
    for(unsigned int i=0;i<seq.size();i++) { cout << seq[i]; }
    cout << "\n";
}

string Sequence::printToString() {
    
    string output;
    for(unsigned int i=0;i<seq.size();i++) { output += seq[i]; }
    
    return output;
}

//If sequence is -ACGT-- print only ACGT
string Sequence::printToStringNoHyphens() {
    
	string output;
	for(unsigned int i=0;i<seq.size();i++) { 
	
		if(seq[i] != '-') {
			output += seq[i];
		}
	}
    
	return output;
}
//Non-end-inclusive, i.e. gives [start,finish)
Sequence Sequence::subset(int start, int finish) {
    
    Sequence result;
    for(int i=start; i<finish;i++) {
        result.addBase(seq[i]);
    }
    return result;
}

void Sequence::addBaseFront(char b) {
    seq.insert(seq.begin(),b);
}

