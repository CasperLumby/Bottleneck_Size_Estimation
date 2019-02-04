//Include guard
#ifndef DATAPHMGGEN_H
#define DATAPHMGGEN_H


//Forward declared dependencies


//Included dependencies
#include "data.h"
#include "pathPHMG.hpp"
#include "dataPHgen.hpp"
#include "diploidSequence.h"
#include "model.hpp"
#include <vector>
#include <string>


//Full haplotype multi gene
class DataPHMGgen : public Data {
private:
	//Internal variables
	PathPHMG path; //Path to data
	std::vector<bool> genePresent; //Describes if data is present for each gene (e.g. dim=8 for flu)
	std::vector<DataPHgen> genes;
    

	    
public:
	DataPHMGgen(); //Constructor
	~DataPHMGgen(); //Deconstructor
    
	void readData(Path * p);
	void simulateData(SimParam *sp);
	int getNumGenes();
	DataPHgen* getGene(int index);
	bool getGenePresent(int g);
	void setGenePresent(int g, bool b);
    
};

#endif
