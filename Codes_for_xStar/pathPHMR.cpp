//
//  pathPHMR.cpp
//  FluTransmissionProject
//
//  Created by Casper Lumby on 11/12/2015.
//  Copyright © 2015 Casper Lumby. All rights reserved.
//

#include "pathPHMR.hpp"

using namespace std;

//Define default constructor
PathPHMR::PathPHMR() : importHaps(false) { /* do nothing */ }
PathPHMR::~PathPHMR() {} //Deconstructor

PathPHMG PathPHMR::getPHMGPath(int index) { return PHMGpaths[index]; }

void PathPHMR::addPHMGPath(PathPHMG& p) { PHMGpaths.push_back(p); }

int PathPHMR::getNumOfPHMGPaths() { return (int) PHMGpaths.size(); }

void PathPHMR::setImportHapsFreq(string& freq) {

	importHaps = true;
	importHapsFreq = freq;

}

void PathPHMR::loadFluReplicateFilePaths(vector<string> & repFolders) {
    
	//Loop over replicates
	for(unsigned int i=0; i<repFolders.size(); i++) {

		PathPHMG currentRepPath;
		currentRepPath.setPhysicalPosMatterToTrue(); //Relevant for replicate analysis
		
		//Import full haps
		if(importHaps == true) {

			currentRepPath.setImportHapsFreq(importHapsFreq);	
		}

		currentRepPath.loadFluFilePaths(repFolders[i]);
		PHMGpaths.push_back(currentRepPath);
	}

}

