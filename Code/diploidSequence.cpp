//
//  diploidSequence.cpp
//  TransmissionProject
//
//  Created by Casper Lumby on 28/10/2015.
//  Copyright © 2015 Casper Lumby. All rights reserved.
//

#include <stdio.h>
#include "iostream"
#include "diploidSequence.h"

using namespace std;

DiploidSequence::DiploidSequence() {}
DiploidSequence::DiploidSequence(int length) : majorSeq(Sequence(length)), minorSeq(Sequence(length)) { }
DiploidSequence::DiploidSequence(Sequence & ma, Sequence & mi) : majorSeq(ma), minorSeq(mi) { }
DiploidSequence::DiploidSequence(string ma, string mi) : majorSeq(Sequence(ma)), minorSeq(Sequence(mi)) { }
DiploidSequence::~DiploidSequence() {}


//Getters and setters
void DiploidSequence::setMajor(Sequence & ma) { majorSeq = ma; }
void DiploidSequence::setMinor(Sequence & mi){ minorSeq = mi; }
Sequence DiploidSequence::getMajor() { return majorSeq; }
Sequence DiploidSequence::getMinor() { return minorSeq; }
void DiploidSequence::setMajor(int index, char ma) { majorSeq.setBase(index, ma); }
void DiploidSequence::setMinor(int index, char mi) { minorSeq.setBase(index, mi); }
void DiploidSequence::print() {
    
    for(int i=0;i<majorSeq.getLength();i++) {
        cout << majorSeq.getBase(i) << "/" << minorSeq.getBase(i) << ",";
    }
    cout << "\n";
    
}

string DiploidSequence::printToString() {
    
	string output;
	for(int i=0;i<majorSeq.getLength();i++) { 
		output += majorSeq.getBase(i);
		output += "/";
		output += minorSeq.getBase(i);
		if(i<majorSeq.getLength()-1) {
			output += " , "; //Don't add after the last pair of bases
		}
	}
    
    return output;
}

int DiploidSequence::getLength() {
	
	return majorSeq.getLength();
}

void DiploidSequence::removeBase(int index) {

	majorSeq.removeBase(index);	
	minorSeq.removeBase(index);
}
