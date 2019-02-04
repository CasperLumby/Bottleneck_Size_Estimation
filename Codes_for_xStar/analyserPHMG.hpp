//Include guard
#ifndef ANALYSERPHMG_HPP
#define ANALYSERPHMG_HPP


//Forward declared dependencies


//Included dependencies
#include "analyser.h"
#include "analyserPH.hpp"
#include "dataPHMG.hpp"
#include <vector>
#include <gsl/gsl_rng.h>
#include "model.hpp"
#include "anaParam.h"
#include <string>

class AnalyserPHMG : public Analyser {
private:
    DataPHMG *data;
	bool warningFound;


	//HapFitG related variables
	std::vector<std::vector<double> > hapFitGMG;
    
public:
    AnalyserPHMG(); //Constructor
    ~AnalyserPHMG(); //Deconstructor
    void loadData(Data *d);
    void runAnalysis(AnaParam *ap);
    
    int getNumGenes();
    
    
	std::vector<std::vector<double> > computeHapFitGrowth(std::vector<std::vector<Sequence> > &fullHaps, int deltaDays);

	std::vector<Model> analyseSelection(AnaParam *ap, std::vector<std::vector<Sequence> > &fullHaps, gsl_rng *r, std::vector<std::vector<double> >& qBmean, std::vector<gsl_matrix* >& qBvar);
	
	int analyseLeonard(std::string method);
    
	void analyseSelectionCombination(AnaParam *ap, std::vector<std::vector<Sequence> > &fullHaps, gsl_rng *r, Model &M);
	
	//FilePathFolder argument is optional. If not supplied, required frequency haps limit will not be printed.
	//This is useful if used to generate haplotypes only.
	//Repindex is optional
	std::vector<std::vector<Sequence> > filterHaps3(gsl_rng* r, double C, std::string filePathFolder = "", int repIndex = -1); 
    
};


#endif
