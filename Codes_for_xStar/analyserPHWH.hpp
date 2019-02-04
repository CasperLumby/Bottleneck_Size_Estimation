//
//  analyserPHWH.hpp
//
//

#ifndef analyserPHWH_hpp
#define analyserPHWH_hpp

#include <stdio.h>
#include "data.h"
#include "analyser.h"
#include "dataPHgen.hpp"
#include "gsl/gsl_rng.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

class AnalyserPHWH : public Analyser {
private:
    DataPHgen *data;
    
public:
    AnalyserPHWH(); //Constructor
    ~AnalyserPHWH(); //Deconstructor
    void loadData(Data *d);
    void runAnalysis(AnaParam *ap);
    
        
	// * = slightly deprecated, but kept for completeness

	void optimiseqMean(int numFHs, gsl_rng *r, int C, std::vector<std::vector<double> >* qMean, std::vector<double>* pL); 

	//Likelihood calculators
    double computeLmean(std::vector<double>* qMean, gsl_rng* r, double C, int timeIndex, bool print);

	double computeLmeanDirMult(std::vector<double>* qMean, double C, int timeIndex, bool print);
    
};


#endif /* analyserPH_hpp */
