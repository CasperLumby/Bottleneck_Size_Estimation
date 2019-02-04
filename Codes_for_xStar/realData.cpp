#include "dataPHMG.hpp"
#include "dataPHMR.hpp"
#include "pathPHMG.hpp"
#include "pathPHMR.hpp"
#include "simParam.h"
#include "anaParam.h"
#include "analyserPHMG.hpp"
#include "analyserPHMR.hpp"
#include "sequence.h"
#include "misc.h"
#include "dataPH.hpp"
#include "simParamPH.h"
#include "analyserPH.hpp"
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

	//Get BIC penalty value
	double BICpenalty = -1;
	const string BICString = input.getCmdOption("-BIC");
	if (!BICString.empty()){
		cout << "BIC penalty input was: " << BICString << "\n";
		BICpenalty = atof(BICString.c_str());
	} else {
		BICpenalty = 10; //Default
		cout << "BICpenalty set to default: " << BICpenalty << "\n";
	}

	//Maximum number of fitted parameters
	int maxNumFittedParams = -1;
	const string maxNumFittedParamsString = input.getCmdOption("-maxParams");
	if (!maxNumFittedParamsString.empty()){
		cout << "Max number of fitted parameter input was: " << maxNumFittedParamsString << "\n";
		maxNumFittedParams = atoi(maxNumFittedParamsString.c_str());
	} else {
		maxNumFittedParams = 5; //Default
	}

	//Get maxNt used in analysis
	int maxNt=-1;
	const string maxNtString = input.getCmdOption("-maxNt");
	if (!maxNtString.empty()) {
		cout << "maxNt input was: " << maxNtString << "\n";
		maxNt = atoi(maxNtString.c_str());
	} else {
		maxNt=10000; //Default
	}


	


	//Get output name (output folder location specified by input)
	const string outputNameString = input.getCmdOption("-o"); //E.g. date or other relevant info


	//Get dataset folder
	//(defined as folders after /rds-d3/user/ckl35/hpc-work/Datasets/)
	//Disregard final "/"!
	stringstream ssInput, ssOutput;
	ssInput << "/rds-d3/user/ckl35/hpc-work/Datasets/";
	ssOutput << "/rds-d3/user/ckl35/hpc-work/Output/";
	const string dataset = input.getCmdOption("-i");
	if(!dataset.empty()) {

		//Add dataset to input
		ssInput << dataset << "/";

		//Create relevant output folders, i.e. split dataset string by "/"
		stringstream ssDataset;
		ssDataset << dataset << "/C_" << c << "_BICpenalty_" << BICpenalty << "_" << outputNameString << "/" << "/Seed_" << seed << "/";
		string folder;
		while(getline(ssDataset, folder, '/')) {

			ssOutput << folder << "/";
			string outputFolder = ssOutput.str();
			//If output folder doesn't exist, create it
			if(fileExists(outputFolder)!= true) {
				cout << "Ouput folder: " << outputFolder << " created.\n";
				mkdir(outputFolder.c_str(),0775);		
			}
		}

	} else {

		cout << "Error in getting the correct filepath!\n";
		return -1;
	}
	string inputFolder = ssInput.str();
	string outputFolder = ssOutput.str();


	//Get folder with information on within host selection, if present
	//(defined as folders after /scratch/ckl35/Datasets/)
	//Folder should contain e.g. "HA", "NA", etc. subfolders
	stringstream ssWithinHost;
	ssWithinHost << "/rds-d3/user/ckl35/hpc-work/Datasets/";
	const string whFolder = input.getCmdOption("-wh");
	bool usePhysicalPositions;
	if(!whFolder.empty()) {

		//Add within host folder to string stream
		ssWithinHost << whFolder << "/";

		cout << "Within host selection folder input: " << ssWithinHost.str() << "\n";

		usePhysicalPositions = true;
		cout << "Need to use physical positions as within host selection information used. Ensure Positions.dat exists.\n";

	} else { //No information inputted

		ssWithinHost.str(""); //If no information, make string empty
		usePhysicalPositions = false;
	}
	string withinHostFolder = ssWithinHost.str();

	//Get minimum required read depth. Setting to 0 results in unfiltered outcome.
	int minReadDepth = -1;
	const string minReadDepthString = input.getCmdOption("-minRD");
	if (!minReadDepthString.empty()) {
		cout << "min read depth input was: " << minReadDepthString << "\n";
		minReadDepth = atoi(minReadDepthString.c_str());
	} else {
		minReadDepth = 100; //Default
	}

	//Follow-up: If min read depth is non zero, we need to use physical positions
	if(minReadDepth > 0) {
		usePhysicalPositions = true;
	}
	
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


	//Filter haps within the C++ code
	const string filterString = input.getCmdOption("-filter");
	int filterHapsMethod=0; //Default is 0, i.e. no filtering
	if(!filterString.empty()) {

		filterHapsMethod = atoi(filterString.c_str());
	}

	//Get termination point
	int terminationPoint = 0; //Default is not to terminate
	const string terminationString = input.getCmdOption("-term");
	if (!terminationString.empty()){
		cout << "Termination point input was: " << terminationString << "\n";
		terminationPoint = atoi(terminationString.c_str());
	}

	const string noVarString = input.getCmdOption("-noVar");
	bool noVar=false; //As standard, we use variance
	if(!noVarString.empty()) {

		if(noVarString.compare("True")==0 || noVarString.compare("true")==0) {

			noVar = true;
		} else if (noVarString.compare("false")==0 || noVarString.compare("False")==0) {

			noVar = false;
		} else {

			cout << "Error with noVar input: Input has to be true or false.\n"; exit(1);
		}	
	}

	const string meanOnlyString = input.getCmdOption("-meanOnly");
	bool meanOnly=false; //As standard, we use variance
	if(!meanOnlyString.empty()) {

		if(meanOnlyString.compare("True")==0 || meanOnlyString.compare("true")==0) {

			meanOnly = true;
		} else if (meanOnlyString.compare("false")==0 || meanOnlyString.compare("False")==0) {

			meanOnly = false;
		} else {

			cout << "Error with meanOnly input: Input has to be true or false.\n"; exit(1);
		}	
	}

	cout << "Analysing dataset " << dataset << " using seed " << seed << "\n";

	AnaParam ap = AnaParam();
	ap.setC(c);
	ap.setBICpenalty(BICpenalty);
	ap.setMaxNumFittedParams(maxNumFittedParams);
	ap.setSeed(seed);
	ap.loadOutputFolder(outputFolder);
	ap.setFilterHaplotypesMethod(filterHapsMethod);
	ap.setTerminateEarly(terminationPoint);
	ap.setMaxNt(maxNt);
	ap.setNoVar(noVar);
	ap.setMeanOnly(meanOnly);
	
	if(repsPresent == false) {

		PathPHMG pPHMG;
		if(usePhysicalPositions == true) {
			pPHMG.setPhysicalPosMatterToTrue();
		}
		if(importHaps == true) {
			pPHMG.setImportHapsFreq(importHapsFreq);
		}
		pPHMG.setMinReadDepth(minReadDepth);
		pPHMG.loadFluFilePaths(inputFolder);
		DataPHMG dPHMG;
		dPHMG.readData(&pPHMG);

		if(withinHostFolder.compare("") != 0) {

//Temporarily commented out? Need to decide.			dPHMG.loadWithinHostSelection(withinHostFolder);
			ap.loadWithinHostSelectionFolder(withinHostFolder);
		}
		ap.loadUseGreedyMMS(false);
		AnalyserPHMG aPHMG;
		aPHMG.loadData(&dPHMG);
		aPHMG.runAnalysis(&ap);

	} else {

		//Get infoFile from inputfolder
		stringstream ss;
		ss << inputFolder << "RepInfo.dat";
		string infoFilePath = ss.str();

		//Get replicate sub folders
		vector<string> repFolders;
		ifstream infoFile;
		infoFile.open(infoFilePath.c_str());
		string line;
		while(getline(infoFile, line)) {

			ss.str("");
			ss << inputFolder << line << "/";
			string repFolder = ss.str();
			repFolders.push_back(repFolder);
		}


		PathPHMR pPHMR;
		if(importHaps == true) {
			pPHMR.setImportHapsFreq(importHapsFreq);
		}
		pPHMR.loadFluReplicateFilePaths(repFolders);
		DataPHMR dPHMR;
		dPHMR.readData(&pPHMR);

		if(withinHostFolder.compare("") !=0) {

//Temporarily commented out? Need to decide.			dPHMR.loadWithinHostSelection(withinHostFolder);
			ap.loadWithinHostSelectionFolder(withinHostFolder);

		}

		ap.loadUseGreedyMMS(false);
		AnalyserPHMR aPHMR;
		aPHMR.loadData(&dPHMR);
		aPHMR.runAnalysis(&ap);
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

