
//Forward declared dependencies

//Included dependencies
#include "pathPH.h"
#include <stdio.h>
#include <iostream>
#include <sstream>
#include "misc.h"
#include <sys/stat.h>

using namespace std;

//Define default constructor
PathPH::PathPH() : physicalPosMatter(false), importHaps(false), importHapsMahan(false), minReadDepth(100), filterData(true), timesFile("") { /* do nothing */ }
PathPH::~PathPH() {} //Deconstructor

vector<string> PathPH::getHapPaths() { return hapPaths; }
vector<string> PathPH::getContribPaths() { return contribPaths; }
vector<string> PathPH::getLociPaths() { return lociPaths; }
string PathPH::getMultiLocusTrajFilePath() { return multiLocusTrajFilePath; }
string PathPH::getSingleLocusTrajFilePath() { return singleLocusTrajFilePath; }
bool PathPH::getPhysicalPosMatter() { return physicalPosMatter; } //False as standard
void PathPH::setPhysicalPosMatterToTrue() { physicalPosMatter = true; }
string PathPH::getPhysicalPosPath() { return physicalPosPath; }
void PathPH::setMinReadDepth(int minRD) { minReadDepth = minRD; }
int PathPH::getMinReadDepth() { return minReadDepth; }
void PathPH::setFilterData(bool fd) { filterData = fd; }
bool PathPH::getFilterData() { return filterData; }

//Sets the frequency at which imported haps are called.
void PathPH::setImportHapsFreq(string& freq) {

	importHaps = true;
	importHapsFreq = freq;

}

void PathPH::setImportHapsMahanFile(string& path) {

	importHapsMahan = true;
	importHapsMahanFile = path;
	
}

string PathPH::getImportHapsFile() {

	if(importHaps == true) {

		return importHapsFile;

	} else {

		cout << "Import haps file is not loaded, so can't be obtained from PathPH. Exiting.\n";
		exit(1);

	}

}

string PathPH::getImportHapsMahanFile() {

	if(importHapsMahan == true) {

		return importHapsMahanFile;

	} else {

		cout << "Import haps Mahan file is not loaded, so can't be obtained from PathPH. Exiting.\n";
		exit(1);

	}
}

string PathPH::getTimesFile() {

	return timesFile;
}

bool PathPH::getImportHaps() { return importHaps; }
bool PathPH::getImportHapsMahan() { return importHapsMahan; }

/*
 * H5N1 data. One load method for each gene. Not the most compact code. Consider simplifying
 */
void PathPH::loadH5N1HAGeneFilePaths(string& folder) {
    
    hapPaths.clear();
    contribPaths.clear();
    lociPaths.clear();
    
    stringstream ssHap;
    stringstream ssCon;
    stringstream ssLoc;
    string currentHapFile;
    string currentConFile;
    string currentLocFile;
    
    
    for(int i=0; i<9999999;i++) {
        
        //Current haplotype file
        ssHap.str("");
        ssHap  << folder << "HA/Hap_data" << i <<".dat";
        currentHapFile = ssHap.str();
        
        //Current contributions file
        ssCon.str("");
        ssCon  << folder << "HA/Contribs" << i <<".dat";
        currentConFile = ssCon.str();
        
        //Current loci file
        ssLoc.str("");
        ssLoc  << folder << "HA/Loci" << i <<".dat";
        currentLocFile = ssLoc.str();
        
        if(fileExists(currentHapFile)) {
            
            hapPaths.push_back(currentHapFile);
            contribPaths.push_back(currentConFile);
            lociPaths.push_back(currentLocFile);
            
        } else { // no more files to add, so break out of loop
            break;
        }
        
    }


	//Add path to physical positions file, if requested
	if(physicalPosMatter==true) {
		stringstream ssPos;
		ssPos << folder << "HA/Positions.dat";
		physicalPosPath = ssPos.str();
	}
 
    cout << "H5N1 HA filepaths loaded.\n";
}

void PathPH::loadH5N1NAGeneFilePaths(string& folder) {
    
    hapPaths.clear();
    contribPaths.clear();
    lociPaths.clear();
    
    stringstream ssHap;
    stringstream ssCon;
    stringstream ssLoc;
    string currentHapFile;
    string currentConFile;
    string currentLocFile;
    
    
    for(int i=0; i<9999999;i++) {
        
        //Current haplotype file
        ssHap.str("");
        ssHap  << folder << "NA/Hap_data" << i <<".dat";
        currentHapFile = ssHap.str();
        
        //Current contributions file
        ssCon.str("");
        ssCon  << folder << "NA/Contribs" << i <<".dat";
        currentConFile = ssCon.str();
        
        //Current loci file
        ssLoc.str("");
        ssLoc  << folder << "NA/Loci" << i <<".dat";
        currentLocFile = ssLoc.str();
        
        if(fileExists(currentHapFile)) {
            
            hapPaths.push_back(currentHapFile);
            contribPaths.push_back(currentConFile);
            lociPaths.push_back(currentLocFile);
            
        } else { // no more files to add, so break out of loop
            break;
        }
        
    }
    
	//Add path to physical positions file, if requested
	if(physicalPosMatter==true) {
		stringstream ssPos;
		ssPos << folder << "NA/Positions.dat";
		physicalPosPath = ssPos.str();
	}
 

    cout << "H5N1 NA filepaths loaded.\n";
}

void PathPH::loadH5N1M1GeneFilePaths(string& folder) {
    
    hapPaths.clear();
    contribPaths.clear();
    lociPaths.clear();
    
    stringstream ssHap;
    stringstream ssCon;
    stringstream ssLoc;
    string currentHapFile;
    string currentConFile;
    string currentLocFile;
    
    
    for(int i=0; i<9999999;i++) {
        
        //Current haplotype file
        ssHap.str("");
        ssHap  << folder << "MP/Hap_data" << i <<".dat";
        currentHapFile = ssHap.str();
        
        //Current contributions file
        ssCon.str("");
        ssCon  << folder << "MP/Contribs" << i <<".dat";
        currentConFile = ssCon.str();
        
        //Current loci file
        ssLoc.str("");
        ssLoc  << folder << "MP/Loci" << i <<".dat";
        currentLocFile = ssLoc.str();
        
        if(fileExists(currentHapFile)) {
            
            hapPaths.push_back(currentHapFile);
            contribPaths.push_back(currentConFile);
            lociPaths.push_back(currentLocFile);
            
        } else { // no more files to add, so break out of loop
            break;
        }
        
    }
 
	//Add path to physical positions file, if requested
	if(physicalPosMatter==true) {
		stringstream ssPos;
		ssPos << folder << "MP/Positions.dat";
		physicalPosPath = ssPos.str();
	}
    
    cout << "H5N1 MP filepaths loaded.\n";
}


void PathPH::loadFluGeneFilePaths(string &folder, string &gene) {
    
	hapPaths.clear();
	contribPaths.clear();
	lociPaths.clear();
    
	stringstream ssHap;
	stringstream ssCon;
	stringstream ssLoc;
	string currentHapFile;
	string currentConFile;
	string currentLocFile;
    
	bool filesFound = false;
	for(int i=0; i<9999999;i++) {
        
		//Current haplotype file
		ssHap.str("");
		ssHap  << folder << gene <<"/Hap_data" << i <<".dat";
		currentHapFile = ssHap.str();
        
		//Current contributions file
		ssCon.str("");
		ssCon  << folder << gene << "/Contribs" << i <<".dat";
		currentConFile = ssCon.str();
        
		//Current loci file
		ssLoc.str("");
		ssLoc  << folder << gene << "/Loci" << i <<".dat";
		currentLocFile = ssLoc.str();
        
		if(fileExists(currentHapFile)) {
            
			filesFound = true;
			hapPaths.push_back(currentHapFile);
			contribPaths.push_back(currentConFile);
			lociPaths.push_back(currentLocFile);
            
		} else { // no more files to add, so break out of loop
			break;
		}
        
	}

	//Add path to single_traj file if existing
	stringstream ssSingle;
	ssSingle << folder << gene << "/Single_locus_trajectories.out";
	singleLocusTrajFilePath = ssSingle.str();

	//Add path to multi-traj file if existing
	stringstream ssMulti;
	ssMulti << folder << gene << "/Multi_locus_trajectories.out";
	multiLocusTrajFilePath = ssMulti.str();

	//Add path to physical positions file, if requested
	if(physicalPosMatter==true) {
		stringstream ssPos;
		ssPos << folder << gene << "/Positions.dat";
		physicalPosPath = ssPos.str();
	}

	//Import full haps
	if(importHaps == true) {

		//Define input file
		stringstream ssImportHaps;
		ssImportHaps << folder << gene << "/InferredFilter" << importHapsFreq << ".dat";
		importHapsFile = ssImportHaps.str();
		
	}


	if(filesFound == true) {
		cout << "Flu " << gene << " gene filepaths loaded.\n";
	}

	//If existing, create path to Times.in file. If non-existing, string is empty ("")
	stringstream ssTimes;
	ssTimes << folder << gene << "/Times.in";
	timesFile = ssTimes.str();
	if(fileExists(timesFile) != true) {

		//File doesn't exists, so change string to ""
		timesFile = "";
	}
	
}


