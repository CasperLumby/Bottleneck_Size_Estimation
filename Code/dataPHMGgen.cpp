#include "dataPHMGgen.hpp"
#include <stdio.h>
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "misc.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace std;

//Class related methods + constructors
DataPHMGgen::DataPHMGgen() { } //Constructor does nothing
DataPHMGgen::~DataPHMGgen() {} //Deconstructor

void DataPHMGgen::readData(Path * p) {
    
    path = * static_cast<PathPHMG *>(p);
    
    //Add data gene by gene
    for(int i=0; i<path.getNumOfPHPaths(); i++) {

	DataPHgen dphgen;
        PathPH pph = path.getPHPath(i);
        dphgen.readData(&pph);

	if(dphgen.getNumOfSOPHs() > 0) {
		cout << "Data added for gene " << i << "\n";
	}
        

        
	//Check if gene contains any data
	if(dphgen.getNumOfSOPHs() > 0) {
		genePresent.push_back(true);
	} else {

		cout << "Gene data not present for gene: " << i << "\n";
		genePresent.push_back(false);
	}
        	
	//Add the gene    
	genes.push_back(dphgen);
	//cout << i << " added: ";
 
    }
    
}

void DataPHMGgen::simulateData(SimParam *sp) { } //Implement later

int DataPHMGgen::getNumGenes() {	return (int) genes.size(); }
DataPHgen* DataPHMGgen::getGene(int index) { return &(genes[index]); }

bool DataPHMGgen::getGenePresent(int g) { return genePresent[g]; }
void DataPHMGgen::setGenePresent(int g, bool b) { genePresent[g] = b; }

