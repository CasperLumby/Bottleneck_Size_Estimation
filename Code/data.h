//Include guard
#ifndef DATA_H
#define DATA_H


//Forward declared dependencies

//Included dependencies
#include "path.h"
#include "simParam.h"
#include <vector>


class Data{
    
public:
    virtual void readData(Path * p) = 0;
    virtual void simulateData(SimParam *sp) = 0;
    
};

#endif 
