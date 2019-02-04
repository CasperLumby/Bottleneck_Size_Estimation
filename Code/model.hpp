//
//  model.hpp
//  FluTransmissionProject
//
//  Created by Casper Lumby on 27/01/2016.
//  Copyright © 2016 Casper Lumby. All rights reserved.
//

#ifndef model_hpp
#define model_hpp

#include <stdio.h>
#include "diploidSequence.h"
#include <gsl/gsl_matrix.h>


class Model {
    
    
public:
    
    struct epistasis {
    private:
        int dim; //dim={2,3,4...}, e.g. dim=2 for two way epistasis
        std::vector<int> positions; //Positions for the dim-way interactions
        
    public:
        void addPositions(int pos);
        void setPositions(std::vector<int> &pos);
        void setDim(int d);
	int getDim();
        std::vector<int> getPositions();
	int getPosition(int index);
	void print();
    };
    
    struct modelSingleGene {
    private:
        int numParamsToBeFitted;
        DiploidSequence model; //E.g. {A,-,-,C,-,T,-}
        std::vector<double> selCoefsBest; //E.g. {2,0,0,1.5,0,-1.2,0}
        std::vector<double> selCoefsNew;
        std::vector<epistasis> epiModel; //E.g. vector(two way epistasis, three way epistasis)
        std::vector<double> epiCoefsBest; //E.g. {1.5, -3}
        std::vector<double> epiCoefsNew;
	std::vector<double> qBmeanBest; //Pre transmission frequency
	std::vector<double> qBmeanNew;        
	gsl_matrix* qBvarBest;
	gsl_matrix* qBvarNew;
	bool qBvarBestAllocated, qBvarNewAllocated;
	bool selectionPresent;

    public:
	modelSingleGene();
	~modelSingleGene();
	
        void countNumParamsToBeFitted();
        int getNumParamsToBeFitted();
        void setModel(DiploidSequence & ds);
	void setEpiModel(std::vector<epistasis> &em);
        DiploidSequence getModel(); //Might be worth changing name to something more accurate
        void addEpiModel(epistasis & epi);
        void setSelCoefsNew(std::vector<double> &SCN);
        void setSelCoefsBest(std::vector<double> &SCB);
        void setSelCoefsNew(int index, double value);
        void setSelCoefsBest(int index, double value);
        std::vector<double> getSelCoefsNew();
        double getSelCoefsNew(int index);
        void setSelCoefsNewToBest();
        void setSelCoefsBestToNew();
        void setEpiCoefsNewToBest();
        void setEpiCoefsBestToNew();
        std::vector<double> getSelCoefsBest();
        void setEpiCoefsBest(std::vector<double> &ECB);
        void setEpiCoefsNew(std::vector<double> &ECN);
	std::vector<double> getEpiCoefsNew();
	std::vector<double> getEpiCoefsBest();
        
        std::vector<epistasis> getEpiModel();
        int getEpiModelLength();
        void setEpiCoefsNew(int index, double value);
        void setEpiCoefsBest(int index, double value);
        double getEpiCoefsBest(int index);
        double getEpiCoefsNew(int index);
	bool getSelectionPresent();
	void setSelectionPresent(bool b);
        bool epistasisPresent();
	std::vector<double> computeHapFitTransmission(std::vector<Sequence>& fullHaps);
        
		void setqBmeanBest(std::vector<double> &q);
		void setqBmeanNew(std::vector<double> &q);
		std::vector<double>& getqBmeanBest();
		std::vector<double>& getqBmeanNew();
		void print();
		void printNew();
		void setqBvarBest(gsl_matrix* var);
		void setqBvarNew(gsl_matrix* var);
		gsl_matrix* getqBvarBest();
		gsl_matrix* getqBvarNew();
		void deallocateqB();
    };
    
protected:
    int numParamsToBeFitted;
    double BICbest;
    int NtBest, NtNew;
     std::vector<modelSingleGene> modelAllGenes;
	std::vector<std::vector<int> > selCoefsToBeFitted; //Outer: coefficient number, e.g. {1,2,3,4..}. Inner: Gene and pos, e.g. {1,4} or {3,2}
	std::vector<std::vector<int> > epiCoefsToBeFitted; //Outer: coefficient number. Inner: gene and epi number, e.g. {1,1} and {1,2} or {3,1}
	int numqBtoBeFitted; //Total number of haplotype frequencies to be fitted (across all genes)
	int numNtToBeFitted;
    
	//Number of accepts during update
	std::vector<int> accepts;
	std::vector<int> tries;
	std::vector<double> acceptanceRates;
	int currentlyBeingUpdated; //0 = Nt, 1 = qB, 2 = sel, 3 = epi
	int acceptsOverall;
	int triesOverall;
	double acceptanceRateOverall;
	double selectionCap; //Assumed positive
	

public:
    
    //Constructors
    Model();
    ~Model();
    
   
 
    void countNumParamsToBeFitted();
    int getNumParamsToBeFitted();
    double getBICbest();
    
    
    std::vector<double> computeHapFitTransmission(int gene, std::vector<Sequence> &fullHaps);
    void updateCoefs(gsl_rng *r, double delta);
    void update(gsl_rng *r, double delta, int maxNt);
    void updateSelCoefs(gsl_rng *r, double delta);
    void updateEpiCoefs(gsl_rng *r, double delta);
    void updateNt(gsl_rng *r, double delta, int maxNt);
	void updateqB(gsl_rng *r, double delta);
    void initialiseAndRandomiseCoefs(gsl_rng *r, std::vector<std::vector<Sequence> >& haps);
    bool epistasisPresent();
	bool getSelectionPresent(int gene);
    int getNumGenes();
    int getNtNew();
    int getNtBest();
    void setNtNew(int n);
    void setNtBest(int n);
    void updateBICbest(double BICnew);
	void setqBmeanNew(std::vector<std::vector<double> > &qB);
	void setqBmeanBest(std::vector<std::vector<double> > &qB);
	void setqBmeanNew(int gene, std::vector<double> &qB);
	void setqBmeanBest(int gene, std::vector<double> &qB);
	std::vector<std::vector<double> > getqBmeanNew();    
	std::vector<std::vector<double> > getqBmeanBest();
	std::vector<double>& getqBmeanNew(int gene);
	std::vector<double>& getqBmeanBest(int gene);
	void setqBvarBest(int gene, gsl_matrix* var);
	void setqBvarNew(int gene, gsl_matrix* var);
	gsl_matrix* getqBvarBest(int gene);
	gsl_matrix* getqBvarNew(int gene);
	void setBICbest(double BICnew);
	void print();
	void printNew();
	void resetAccepts();
	double updateAcceptanceRates(); //Returns the overall acceptance rate
	void setSelCoefsToBeFitted(std::vector<std::vector<int> >& s);
	void setSelectionCap(double value);

    
    //Means of setting and getting models
    void addModelSingleGene(modelSingleGene &MSG);
	void setModelSingleGene(int index, modelSingleGene &MSG);
    modelSingleGene getModelSingleGene(int index);


	bool isASubsetOf(Model &M);

	void deallocateqB();
    
};

#endif /* model_hpp */
