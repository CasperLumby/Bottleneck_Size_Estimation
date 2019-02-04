#include "dataPHMGgen.hpp"
#include "dataPHMR.hpp"
#include "pathPHMG.hpp"
#include "pathPHMR.hpp"
#include "simParam.h"
#include "anaParam.h"
#include "analyserPHMGWH.hpp"
#include "analyserPHMR.hpp"
#include "sequence.h"
#include "misc.h"
#include "simParamPH.h"
#include "inputParser.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <sys/stat.h>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif


using namespace std;

int main (int argc, char* argv[]) {

	//Get start time
	int timeStart = omp_get_wtime();

	//Initiate the input parser from the command line inputs
	InputParser input(argc, argv);

	//Get analysis seed
	int seed = -1;
	const string seedString = input.getCmdOption("-as");
	if (!seedString.empty()){
		cout << "Analysis seed input was: " << seedString << "\n";
		seed = atoi(seedString.c_str());
	} else {
		seed =  1; //If not defined, set to 1
	}

	//Get C value
	double c = -1;
	const string cString = input.getCmdOption("-C");
	if (!cString.empty()){
		cout << "C input was: " << cString << "\n";
		c = atof(cString.c_str());
	} else {
		c = 200; //Default
	}




	//Get dataset folder in full form, e.g.
	// /rds-d3/user/ckl35/hpc-work/Datasets/Moncla/FolderA/FolderB/FolderC
	//Disregard final "/"!
	string inputFolder;
	stringstream ssInputF;
	const string datasetF = input.getCmdOption("-if");
	if(!datasetF.empty()) {

		//Add dataset to input
		ssInputF << datasetF << "/";
		inputFolder = ssInputF.str();

	} else {

		cout << "Error in getting the correct input filepath!\n";
		return -1;
	}

	//Create output folder based on full path
	//e.g. /rds-d3/user/ckl35/hpc-work/Output/Moncla/Analysis/FolderA/FolderB/FolderC
	string outputFolder;
	stringstream ssOutputF;
	const string outputFolderStringF = input.getCmdOption("-of");
	if(!outputFolderStringF.empty()) {

		
		//Create relevant output folders, i.e. split tempFolders string by "/"
		stringstream ssTempFolder;
		ssTempFolder << outputFolderStringF << "/Seed_" << seed << "/";
		string folder;
		while(getline(ssTempFolder, folder, '/')) {

			ssOutputF << folder << "/";
			string outputFolderCurrent = ssOutputF.str();
			//If output folder doesn't exist, create it
			if(fileExists(outputFolderCurrent)!= true) {
				cout << "Ouput folder: " << outputFolderCurrent << " created.\n";
				mkdir(outputFolderCurrent.c_str(),0775);		
			}
		}
		outputFolder = ssOutputF.str();

	}
	cout << "Outputfolder is : " << outputFolder << "\n";



	const string reps = input.getCmdOption("-reps");
	bool repsPresent=false; //As standard, data is single replicate data
	if(!reps.empty()) {

		if(reps.compare("True")==0 || reps.compare("true")==0) {

			repsPresent = true;
		}	
	}

	//Flag to use imported haps rather than letting the program infer full haplotypes itself.
	//Can be useful in order to limit the number of full haplotypes used.
	//The import method assumes that the filtered haplotypes are placed in the same folder as the
	//data itself. The import method expects a number, e.g. either 01 for 0.1% or 1 for 1%.
	//The import method expects that the haplotype files are named e.g. InferredFilter01Haps_Gene_0.dat
	//for single replicate data or InferredFilter01Haps_Rep0_Gene0.dat for replicate data.
	const string importHapsFreqString = input.getCmdOption("-importHaps");
	string importHapsFreq = "";
	bool importHaps = false;
	if(!importHapsFreqString.empty()) {

		importHapsFreq = importHapsFreqString.c_str();
		importHaps = true;
	} else {

		importHaps = false; //Default, don't use imported haplotypes
	}

	cout << "Analysing dataset " << datasetF << " using seed " << seed << "\n";

	AnaParam ap = AnaParam();
	ap.setC(c);
	ap.setSeed(seed);
	ap.loadOutputFolder(outputFolder);
	
	if(repsPresent == false) {

		PathPHMG pPHMG;
		if(importHaps == true) {
			pPHMG.setImportHapsFreq(importHapsFreq);
		}
		pPHMG.loadFluFilePaths(inputFolder);
		DataPHMGgen dPHMGgen;
		dPHMGgen.readData(&pPHMG);

		
		AnalyserPHMGWH aPHMGWH;
		aPHMGWH.loadData(&dPHMGgen);
		aPHMGWH.runAnalysis(&ap);

	} else {
		
		//Implement later
	}

	//Print elapsed time to file 
	int timeFinish = omp_get_wtime();
	cout << "Elapsed time: " << timeFinish - timeStart << "\n";

	stringstream ssTime;
	ssTime << outputFolder << "Elapsed_Time.dat";
	string filePath = ssTime.str();
	ofstream outputFile;
	outputFile.open(filePath.c_str());
	outputFile << timeFinish - timeStart << "\n";
	outputFile.close();

	return 0;
}

