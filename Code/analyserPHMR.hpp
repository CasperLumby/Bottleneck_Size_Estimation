//Include guard
#ifndef ANALYSERPHMR_HPP
#define ANALYSERPHMR_HPP


//Forward declared dependencies


//Included dependencies
#include "analyser.h"
#include "analyserPH.hpp"
#include "analyserPHMG.hpp"
#include "dataPHMG.hpp"
#include "dataPHMR.hpp"
#include <vector>
#include <gsl/gsl_rng.h>
#include "modelMR.hpp"
#include "anaParam.h"

class AnalyserPHMR : public Analyser {
private:
    DataPHMR *data;
	std::vector<std::vector<std::vector<double> > > hapFitGMR;
    
public:
    AnalyserPHMR(); //Constructor
    ~AnalyserPHMR(); //Deconstructor
    void loadData(Data *d);
    void runAnalysis(AnaParam *ap);
   

//    int getNumReplicates();
    
    
    
	std::vector<ModelMR> analyseSelection(AnaParam *ap, std::vector<std::vector<std::vector<Sequence> > > &fullHaps, gsl_rng *r, std::vector<std::vector<std::vector<double> > > &qBmean, std::vector<std::vector<gsl_matrix*> >& qBvar);
    
	void analyseSelectionCombination(AnaParam *ap, std::vector<std::vector<std::vector<Sequence> > > &fullHaps, gsl_rng *r, ModelMR &M);
    
	std::vector<std::vector<int> > findAllPos();
	std::vector<std::vector<std::vector<int> > > createMapFromAllPosToReps(std::vector<std::vector<int> > &allPos);
	std::vector<DiploidSequence> createCollapsedFullHapsReplicates(std::vector<std::vector<std::vector<Sequence> > > &fullHaps, std::vector<std::vector<std::vector<int> > > &map);
	std::vector<std::vector<int> > findSharedPos();
	std::vector<std::vector<std::vector<int> > > createMap(std::vector<std::vector<int> > &sharedPos);

	//FilePathFolder argument is optional. If not supplied, required frequency haps limit will not be printed.
	//This is useful if used to generate haplotypes only.
	std::vector<std::vector<std::vector<Sequence> > > filterHaps3(gsl_rng* r, double C, std::string filePathFolder = ""); 
 
};



#endif
