//Include guard
#ifndef SIMPARAMPH_H
#define SIMPARAMPH_H


//Forward declared dependencies

//Included dependencies
#include "simParam.h"


class SimParamPH : public SimParam {
protected:
    int dim;
    double C;
	std::string pathToFolder, repID, format, filterData;
	bool repScenario;
	int geneIndex;
	int numGenerations; //Number of generations per day of WH replications
	int deltaDays; //Number of days between donor and recipient sampling
	int growthFactor;

public:
    
    //Constructors
    SimParamPH();
   
    //Deconstructors
    ~SimParamPH();
    

	//Setters
	void setDim(int d);
	void setC(double c);
	void setPathToFolder(std::string path);
	void setRepID(std::string ID);
	void setRepScenario(bool b);
	void setGeneIndex(int g);
	void setFormat(std::string f);
	void setFilterData(std::string fd);
	void setGrowthFactor(int gf);
	void setNumGenerations(int ng);
	void setDeltaDays(int dd);

	//Getters
    int getDim();
    double getC();
	std::string getPathToFolder();
	std::string getRepID();
	bool getRepScenario();
	int getGeneIndex();
	std::string getFormat();
	std::string getFilterData();
	int getGrowthFactor();
	int getNumGenerations();
	int getDeltaDays();
    
};

#endif 
