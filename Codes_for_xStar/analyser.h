//Include guard
#ifndef ANALYSER_H
#define ANALYSER_H


//Forward declared dependencies

//Included dependencies
//#include "anaParam.h"
#include "data.h"
#include "anaParam.h"


class Analyser{
    
public:
    virtual void loadData(Data *d) = 0;
    virtual void runAnalysis(AnaParam *ap) = 0;
};

#endif 
