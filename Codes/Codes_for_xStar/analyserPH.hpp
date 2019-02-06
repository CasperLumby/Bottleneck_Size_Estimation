//
//  analyserPH.hpp
//  FluTransmissionProject
//
//  Created by Casper Lumby on 30/11/2015.
//  Copyright © 2015 Casper Lumby. All rights reserved.
//

#ifndef analyserPH_hpp
#define analyserPH_hpp

#include <stdio.h>
#include "data.h"
#include "analyser.h"
#include "dataPH.hpp"
//#include "analyserFH.h"
#include "gsl/gsl_rng.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

class AnalyserPH : public Analyser {
private:
    DataPH *data;
    
public:
    AnalyserPH(); //Constructor
    ~AnalyserPH(); //Deconstructor
    void loadData(Data *d);
    void runAnalysis(AnaParam *ap);
    
    //Haplotype reconstruction methods
    struct par {
        int i1;
        int i2;
    };
    void constructFullHaps(std::vector<Sequence> *phs, std::vector<Sequence> *fhs);
    void OverlapStepShared1 (std::vector<std::vector<char> >& haps);
    void ResetHaps (std::vector<int> incl, std::vector< std::vector<char> >& haps);
    void DuplicateStep (std::vector< std::vector<char> >& haps);
    void OverlapStepShared2 (std::vector< std::vector<char> >& haps);
    void OverlapStepNoShare (std::vector< std::vector<char> >& haps);
    void GetStartFinish(std::vector<par>& sf, std::vector< std::vector<char> >& haps);
    void BuildPaths(std::vector<par>& sf, std::vector< std::vector<char> > haps, std::vector< std::vector<char> >& new_haps);
    void AddPath (int prev, int start, std::vector<char> v, std::vector<par>& sf, std::vector< std::vector<char> >& haps, std::vector< std::vector<char> >& new_haps);
    
    
	// * = slightly deprecated, but kept for completeness

	//Analysis methods
	//Mean but no var
	double computeLSelection(AnaParam *ap, int Nt, std::vector<double> &qBfhFreqs, std::vector<double> &fhHapFit); //*
	double computeLNeutral(AnaParam *ap, int Nt, std::vector<double> &qBfhFreqs); //*

	//Mean and var
	double computeLSelTVar(AnaParam *ap, int Nt, std::vector<double> &qBfhFreqs, gsl_matrix* qBvar, std::vector<double> &fhHapFitT);
	double computeLSelTSelGVar(AnaParam *ap, int Nt, std::vector<double> &qBfhFreqs, gsl_matrix* qBvar, std::vector<double> &fhHapFitT, std::vector<double> &fhHapFitG);
	double computeLNeutralVar(AnaParam *ap, int Nt, std::vector<double> &qBfhFreqs, gsl_matrix* qBvar);
	double computeLNeutralVarAdvanced(AnaParam *ap, int Nt, std::vector<double> &qBfhFreqs, gsl_matrix* qBvar);
	double computeLNeutralSelGVar(AnaParam *ap, int Nt, std::vector<double> &qBfhFreqs, gsl_matrix* qBvar, std::vector<double> &hapFitG);

	//Optimisation algorithms
	std::vector<double> findOptimalFreqFHBefore(int numFHs, double *L, gsl_rng *r, int C, std::vector<double>* logFactStore); //*
	std::vector<double> findOptimalFreqFHBeforeFromMMS(std::vector<std::vector<Sequence> > &MSSfhs, std::vector<Sequence> *fhs, gsl_rng *r, int C, std::vector<double>* logFactStore); //*
	std::vector<double> findOptimalFreqFHBefore(int numFHs, gsl_rng *r, int C, std::vector<double>* logFactStore); //*
	std::vector<double> findOptimalFreqFHAfter(int numFHs, gsl_rng *r, int C, std::vector<double>* logFactStore); //*
	void optimiseqBmeanFromPre(int numFHs, gsl_rng *r, int C, std::vector<double>* qBmean, double* pL); 
	void optimiseqBmeanAndVarFromPre(int numFHs, gsl_rng *r, int C, std::vector<double>* qBmean, gsl_matrix* qBvar, double* pL); 
	void optimiseqBmeanAndVarSA(int numFHs, gsl_rng *r, int C, std::vector<double>* qBmean, gsl_matrix* qBvar, double* pL); 

	//Likelihood calculators
	void getBeforeFrequenciesAndProbabilities(int numFHS, gsl_rng*r, int C, std::vector<double>* logFactStore,std::vector<std::vector<double> > *freqs, std::vector<double> *probs, int numPoints); //*
    void getBeforeLandscape(int numFHS, gsl_rng*r, int C, std::vector<double>* logFactStore,std::vector<std::vector<double> > *freqs, std::vector<double> *likelihoods, int numPoints, double diffFromMax); //*
    double computeLb(std::vector<double> &freqs, gsl_rng *r, int C, std::vector<double>* logFactStore); //*
    double computeLa(std::vector<double> &freqs, gsl_rng *r, int C, std::vector<double>* logFactStore); //*
    double computeLmeanPreOnly(std::vector<double>* qBmean, gsl_rng* r, double C, bool print);
    double computeLmeanAndVarPreOnly(std::vector<double>* qBmean, gsl_matrix* qBvar, gsl_rng* r, double C);


	//FilePathFolder argument is optional. If not supplied, required frequency haps limit will not be printed.
	//This is useful if used to generate haplotypes only.
	//Rep and gene indices also optional. Note, if only one int is supplied, then geneIndex = int, repIndex = -1, hence the ordering.
	std::vector<Sequence> filterHaps3(gsl_rng* r, double C, std::string filePathFolder = "", int geneIndex = -1, int repIndex = -1); 
    

	void analyseLeonard(std::vector<double>& logLNts);
	void analyseLeonardExact(std::vector<double>& logLNts);

    
};


#endif /* analyserPH_hpp */
