//
//  pathPHMG.cpp
//  FluTransmissionProject
//
//  Created by Casper Lumby on 11/12/2015.
//  Copyright © 2015 Casper Lumby. All rights reserved.
//

#include "pathPHMG.hpp"
#include <sstream>
#include <iostream>

using namespace std;

//Define default constructor
PathPHMG::PathPHMG() : physicalPosMatter(false), importHaps(false), importHapsMahan(false), filterData(true), minReadDepth(100) { /* do nothing */ }
PathPHMG::~PathPHMG() {} //Deconstructor

PathPH PathPHMG::getPHPath(int index) { return PHpaths[index]; }

void PathPHMG::addPHPath(PathPH& p) { PHpaths.push_back(p); }

int PathPHMG::getNumOfPHPaths() { return (int) PHpaths.size(); }

bool PathPHMG::getPhysicalPosMatter() { return physicalPosMatter; }
void PathPHMG::setPhysicalPosMatterToTrue() { physicalPosMatter = true; }

//Sets the frequency at which imported haps are called.
void PathPHMG::setImportHapsFreq(string& freq) {

	importHaps = true;
	importHapsFreq = freq;

}

void PathPHMG::setImportHapsMahan(string& path) {

	importHapsMahan = true;
	importHapsMahanPath = path;
}

void PathPHMG::setImportHapsMahanFilename(string& filename) {

	importHapsMahanFilename = filename;
}

void PathPHMG::setMinReadDepth(int minRD) { minReadDepth = minRD; }
int PathPHMG::getMinReadDepth() { return minReadDepth; }

void PathPHMG::setFilterData(bool fd) { filterData = fd; }
bool PathPHMG::getFilterData() { return filterData; }

void PathPHMG::loadH5N1FilePaths(string& folder) {
    
    	PathPH pHA;
	if(physicalPosMatter == true) { //Only if relevant
		pHA.setPhysicalPosMatterToTrue();
	}
	pHA.loadH5N1HAGeneFilePaths(folder);
   	PHpaths.push_back(pHA);
    
    PathPH pNA;
	if(physicalPosMatter == true) { //Only if relevant
		pNA.setPhysicalPosMatterToTrue();
	}
    pNA.loadH5N1NAGeneFilePaths(folder);
    PHpaths.push_back(pNA);
    
    PathPH pM1;
	if(physicalPosMatter == true) { //Only if relevant
		pM1.setPhysicalPosMatterToTrue();
	}
    pM1.loadH5N1M1GeneFilePaths(folder);
    PHpaths.push_back(pM1);
}

void PathPHMG::loadFluFilePaths(string &folder) {

	vector<string> genes {"HA", "MP", "NA", "NP", "NS", "PA", "PB1", "PB2"};

	//Add filepaths for each gene in turn
	for(unsigned int i=0; i<genes.size(); i++) {

		PathPH p;
		if(physicalPosMatter == true) { //Only if relevant
			p.setPhysicalPosMatterToTrue();
		}	
		if(importHaps == true) {
			p.setImportHapsFreq(importHapsFreq);
		}
		if(importHapsMahan == true) {
	
			stringstream ss;
			ss << importHapsMahanPath << "/test_" << i << "/" << importHapsMahanFilename;
			string filePathGene = ss.str();
			p.setImportHapsMahanFile(filePathGene);
			
		}
		p.setMinReadDepth(minReadDepth);
		p.setFilterData(filterData);
		p.loadFluGeneFilePaths(folder, genes[i]);
		PHpaths.push_back(p);
	} 
}

