//Include guard
#ifndef ANAPARAM_H
#define ANAPARAM_H


//Forward declared dependencies

//Included dependencies
#include "sequence.h"
#include "model.hpp"
#include <string>

class AnaParam{
protected:
	double c;
	int maxNt;
	std::vector<double> * logFactStore;
	unsigned long int seed;
	std::vector<Sequence> selVec;
	std::vector<std::vector<double> > selMagVec;
	std::string pathToFolder;
	int numPrePoints;
	bool useGreedyMMS;
	bool useSharedBottleneck; //Shared bottleneck across replicates
	std::vector<bool> withinHostSelectionPresent;
	std::string withinHostSelectionFolder;
	double BICpenalty;
	double selectionCap;
	bool selectionCapPresent;
	bool useIntegralApproach;
	int maxNumFittedParams;
	int terminateEarly; //Default = 0, i.e don't terminate early
	int filterHaplotypesMethod; //No filtering = 0, 1 = by cutoff, 2 = by repeated reduction, 3 = by ensuring that all phs are represented, i.e. a variation of 1
	bool noVar; //True == don't use variance in qB when computing likelihoods (but compute qBmean and qBvar simultaneously)
	bool meanOnly; //True == don't use variance in qB and compute qBmean alone
	int growthFactor;
	bool analyseSL; //Whether or not to analyse using single locus methods
	
	//Within-host selection - single rep
	std::vector<Sequence> selGvec;
	std::vector<std::vector<double> > selGmagVec;
	std::vector<std::vector<Model::epistasis> > epiGvec;
	std::vector<std::vector<double> > epiGmagVec;
	std::vector<std::vector<double> > hapFitG;
	std::vector<std::vector<std::vector<double> > > hapFitGmultiRep;
    
public:
	//Constructors
	AnaParam();
	~AnaParam();
    
	//Getters and setters
	int getMaxNt();
	void setSeed(int s);
	unsigned long int getSeed();
	std::vector<Sequence>* getSelVec();
	std::vector<std::vector<double> >* getSelMagVec();
	double getC();
	void setC(double C);
	std::vector<double> * getLogFactStore();
	int getNumPrePoints();
    
	void loadOutputFolder(std::string & fullPathToFolder);
	void loadWithinHostSelectionFolder(std::string& fullPathToFolder);
	std::string getWithinHostSelectionFolder();
	void loadWithinHostSelection(std::string & fullPathToFolder, std::vector<std::vector<int> >& physicalPos);
private:
	void loadWithinHostSelectionGene(std::string & fullPathToFolder, std::string& gene, std::vector<int>& physicalPos, bool* selFound);
	void loadWithinHostSelectionGeneOld(std::string & fullPathToFolder, std::string& gene, std::vector<int>& physicalPos);
public:
	void loadUseGreedyMMS(bool value);
	std::string getOutputFolder();
	bool getUseGreedyMMS();
	void setUseSharedBottleneck(bool value);
	bool getUseSharedBottleneck();
	bool getWithinHostSelectionPresent(int gene);
	bool getWithinHostSelectionPresent();
	void setWithinHostSelectionPresent(int gene, bool b);
	void addWithinHostSelectionPresent(bool b);
	
	std::vector<Sequence> getSelGvec();
	void setSelGvec(std::vector<Sequence>& SGV);
	std::vector<std::vector<double> > getSelGmagVec();
	void setSelGmagVec(std::vector<std::vector<double> >& SGMV);
	std::vector<std::vector<Model::epistasis> > getEpiGvec();
	void setEpiGvec(std::vector<std::vector<Model::epistasis> >& EGV);
	std::vector<std::vector<double> > getEpiGmagVec();
	void setEpiGmagVec(std::vector<std::vector<double> >& EGMV);
	
	void setHapFitG(std::vector<std::vector<double> > & HFG);
	void setHapFitG(std::vector<std::vector<std::vector<double> > > & HFG); //Multi rep version
	std::vector<double> getHapFitG(int gene);
	std::vector<double> getHapFitG(int rep, int gene); //Multi rep version
	void setBICpenalty(double p);
	double getBICpenalty();
	void setSelectionCap(double value);
	double getSelectionCap();
	bool getSelectionCapPresent();
	bool getUseIntegralApproach();
	void setUseIntegralApproach(bool b);
	void setMaxNumFittedParams(int num);
	int getMaxNumFittedParams();
	void setMaxNt(int mNt);
	void setTerminateEarly(int t);
	int getTerminateEarly();
	void setFilterHaplotypesMethod(int fhm);
	int getFilterHaplotypesMethod();
	void setNoVar(bool b);
	bool getNoVar();
	void setMeanOnly(bool b);
	bool getMeanOnly();
	void setGrowthFactor(int gf);
	int getGrowthFactor();
	void setAnalyseSL(bool b);
	bool getAnalyseSL();
    
};

#endif 
