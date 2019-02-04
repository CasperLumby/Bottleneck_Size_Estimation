//Include guard
#ifndef ANALYSERPHMGWH_HPP
#define ANALYSERPHMGWH_HPP


//Forward declared dependencies


//Included dependencies
#include "analyser.h"
#include "analyserPH.hpp"
#include "dataPHMGgen.hpp"
#include <vector>
#include <gsl/gsl_rng.h>
#include "model.hpp"
#include "anaParam.h"
#include <string>

class AnalyserPHMGWH : public Analyser {
private:
    DataPHMGgen *data;
	bool warningFound;


public:
    AnalyserPHMGWH(); //Constructor
    ~AnalyserPHMGWH(); //Deconstructor
    void loadData(Data *d);
    void runAnalysis(AnaParam *ap);
    
    int getNumGenes();
    
    
    
};


#endif
