//Include guard
#ifndef SIMPARAM_H
#define SIMPARAM_H


//Forward declared dependencies

//Included dependencies

#include <vector>
#include "sequence.h"
#include "simParam.h"
#include <iostream>
#include <fstream>
#include <sstream>

/*
 * Should really clean up this class.
 * In particular, should be split into two classes - one for multi and one for single genes.
 */
class SimParam {
protected:
    int Nb,Na,numLoci,Nt,numGenes, numGenesWithData, numReps, numLociPerGene, geneLength;
    unsigned long int seed;
	std::ofstream* outputFile;
	bool withinHostSelectionPresent, useExternalReadStatistic, removeMonomorphicSim;
	int filterHaplotypesMethod; //No filtering = 0, 1 = by cutoff, 2 = by repeated reduction, 3 = by ensuring that all phs are represented, i.e. a variation of 1

	//Single gene - Transmission
	Sequence sel;
	std::vector<double> selMagVec;
	std::vector<Sequence> seqVec; //Full haplotypes
	std::vector<std::vector<int> > epiPosVec;
	std::vector<double> epiMagVec;
 
	//Single gene - Growth
	Sequence selG;
	std::vector<double> selMagVecG;
	std::vector<std::vector<int> > epiPosVecG;
	std::vector<double> epiMagVecG;   

	//Multi gene - Transmission
	std::vector<Sequence> selMulti;
	std::vector<std::vector<double> > selMagVecMulti;
	std::vector<std::vector<Sequence> > seqVecMulti;
	std::vector<std::vector<std::vector<int> > > epiPosVecMulti;
	std::vector<std::vector<double> > epiMagVecMulti;
    
	//Multi gene - Growth
	std::vector<Sequence> selMultiG;
	std::vector<std::vector<double> > selMagVecMultiG;
	std::vector<std::vector<std::vector<int> > > epiPosVecMultiG;
	std::vector<std::vector<double> > epiMagVecMultiG;
    

public:
    //Constructors
    SimParam();
    
    virtual ~SimParam();
    

	/*
	 * Setters
	*/	
	void setNa(int NA);
	void setNb(int NB); 
	void setNumLoci(int num); 
	void setNt(int NT);
	void setNumGenes(int num);
	void setNumGenesWithData(int num);
	void setSeed(unsigned long int s);
	void setOutputFile(std::ofstream* o);
	void setNumReps(int num);
	void setNumLociPerGene(int num);
	void setWithinHostSelectionPresent(bool b);
	void setUseExternalReadStatistic(bool b);
	void setGeneLength(int gl);
	void setFilterHaplotypesMethod(int fhm);
	void setRemoveMonomorphicSim(bool b);

	//Transmission
	void setSel(Sequence s);
	void setSelMagVec(std::vector<double> & SMV);
	void setSeqVec(std::vector<Sequence> & SV);
	void setSelMulti(std::vector<Sequence> & SM);
	void setSelMagVecMulti(std::vector<std::vector<double> > & SMVM);
	void setSeqVecMulti(std::vector<std::vector<Sequence> > & SVM);
	void setEpiPosVec(std::vector<std::vector<int> > &EPV);
	void setEpiMagVec(std::vector<double> &EMV);
	void setEpiPosVecMulti(std::vector<std::vector<std::vector<int> > > & EPVM);
	void setEpiMagVecMulti(std::vector<std::vector<double> > &EMVM);
	
	//Growth
	void setSelG(Sequence s);
	void setSelMagVecG(std::vector<double> & SMV);
	void setSelMultiG(std::vector<Sequence> & SM);
	void setSelMagVecMultiG(std::vector<std::vector<double> > & SMVM);
	void setEpiPosVecG(std::vector<std::vector<int> > &EPV);
	void setEpiMagVecG(std::vector<double> &EMV);
	void setEpiPosVecMultiG(std::vector<std::vector<std::vector<int> > > & EPVM);
	void setEpiMagVecMultiG(std::vector<std::vector<double> > &EMVM);
	


  	 /*
	 * Getters
	 */
	int getNa();
	int getNb();
	int getNumLoci();
	int getNt();
	int getNumGenes();
	int getNumGenesWithData();
	unsigned long int getSeed();
	std::ofstream* getOutputFile();
	char getSel(int i);
	int getNumReps();
	int getNumLociPerGene();
	bool getWithinHostSelectionPresent();
	bool getUseExternalReadStatistic();
	int getGeneLength();
        int getFilterHaplotypesMethod();
	bool getRemoveMonomorphicSim();

	//Transmission
	Sequence* getSel();
	std::vector<double>* getSelMagVec();
	std::vector<Sequence> * getSeqVec();
	std::vector<std::vector<int> > *getEpiPosVec();
	std::vector<double> *getEpiMagVec();
	std::vector<Sequence>* getSelMulti();
	std::vector<std::vector<Sequence> >* getSeqVecMulti();
	std::vector<std::vector<double> >* getSelMagVecMulti();
	std::vector<std::vector<std::vector<int> > > *getEpiPosVecMulti();
	std::vector<std::vector<double> > *getEpiMagVecMulti();

	//Growth
	Sequence* getSelG();
	std::vector<double>* getSelMagVecG();
	std::vector<std::vector<int> > *getEpiPosVecG();
	std::vector<double> *getEpiMagVecG();
	std::vector<Sequence>* getSelMultiG();
	std::vector<std::vector<double> >* getSelMagVecMultiG();
	std::vector<std::vector<std::vector<int> > > *getEpiPosVecMultiG();
	std::vector<std::vector<double> > *getEpiMagVecMultiG();
     
};

#endif 
