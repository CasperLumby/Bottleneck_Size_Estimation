//
//  diploidSequence.h
//  TransmissionProject
//
//  Created by Casper Lumby on 28/10/2015.
//  Copyright © 2015 Casper Lumby. All rights reserved.
//

#ifndef diploidSequence_h
#define diploidSequence_h


#include "sequence.h"
#include <string>

class DiploidSequence{
    
private:
    Sequence majorSeq;
    Sequence minorSeq;
    
public:
    DiploidSequence();
    DiploidSequence(int length);
    DiploidSequence(Sequence & ma, Sequence & mi);
    DiploidSequence(std::string ma, std::string mi);
    ~DiploidSequence();
    
    void setMajor(Sequence & mas);
    void setMinor(Sequence & mis);
    void setMajor(int index, char ma);
    void setMinor(int index, char mi);
    Sequence getMajor();
    Sequence getMinor();
    void print();
	std::string printToString();
	int getLength();
	void removeBase(int index);
    
};

#endif /* diploidSequence_h */
