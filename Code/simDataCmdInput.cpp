#include "dataPHMG.hpp"
#include "simParam.h"
#include "anaParam.h"
#include "analyserPHMG.hpp"
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
#include <string.h>
#include <gsl/gsl_rng.h>

#ifdef _OPENMP
        #include <omp.h>
#else
        #define omp_get_thread_num() 0
#endif

using namespace std;

int main (int argc, char* argv[]) {

	//Get start time
	int timeStart = omp_get_wtime();
	
	//Create param objects
	SimParamPH spPHMG;
	AnaParam ap;

	//Initiate the input parser from the command line inputs
	InputParser input(argc, argv);

	//Get simulation mode ("simple" or "real", real is default))
	int simMode = -1;
	const string modeString = input.getCmdOption("-mode");
	if (!modeString.empty()){
		cout << "Simulation mode input was: " << modeString << "\n";
		if(modeString.compare("real")==0 || modeString.compare("Real")==0) {
			simMode = 0;
		
		} else if(modeString.compare("simple")==0 || modeString.compare("Simple")==0) {
			simMode = 1;

		} else {
			cout << "Input to simulation mode was wrong. Should be \"real\" or \"simple\". Input was: " << modeString << "\n";
			exit(1);
		}
	} else {
		simMode = 0; //Default is real
	}

	//Get bottleneck size
	int Nt = -1;
	const string NtString = input.getCmdOption("-Nt");
	if (!NtString.empty()){
		cout << "Bottleneck input was: " << NtString << "\n";
		Nt = atoi(NtString.c_str());
	} else {
		Nt = 100; //Default
	}

	//Get simulation seed
	int simSeed = -1;
	const string simSeedString = input.getCmdOption("-ss");
	if (!simSeedString.empty()){
		cout << "Simulation seed input was: " << simSeedString << "\n";
		simSeed = atoi(simSeedString.c_str());
	} else {
		simSeed = 1; //Default
	}

	//Get analysis seed
	int anaSeed = -1;
	const string anaSeedString = input.getCmdOption("-as");
	bool anaSeedSpecified = false; //Keeps track of whether anaSeed was externally specified or randomly allocated based on simSeed
	if (!anaSeedString.empty()){
		cout << "Analysis seed input was: " << anaSeedString << "\n";
		anaSeed = atoi(anaSeedString.c_str());
		anaSeedSpecified = true;
	} else {

		//If not defined, get random seed based on simSeed
		gsl_rng * rng;
		const gsl_rng_type * T;
		gsl_rng_env_setup();
		T = gsl_rng_default;
		rng = gsl_rng_alloc (T);
		gsl_rng_set(rng, simSeed);

		anaSeed = gsl_rng_uniform_int(rng,1000000); //In theory, should use full integer range, but this should suffice. See https://www.gnu.org/software/gsl/manual/html_node/Sampling-from-a-random-number-generator.html for more info.
		cout << "Analysis seed was: " << anaSeed << "\n";
	}

	//Get number of genes
	int numGenes = -1;
	const string numGenesString = input.getCmdOption("-ng");
	if (!numGenesString.empty()){
		cout << "Number of genes input was: " << numGenesString << "\n";
		numGenes = atoi(numGenesString.c_str());
	} else {
		numGenes = 8; //Default
	}

	
	//Get number of genes with actual data
	int numGenesWithData = -1;
	const string numGenesWithDataString = input.getCmdOption("-ngwd");
	if (!numGenesWithDataString.empty()){
		cout << "Number of genes with data input was: " << numGenesWithDataString << "\n";
		numGenesWithData = atoi(numGenesWithDataString.c_str());
	} else {
		numGenesWithData = 8; //Default
	}

	//Get geneLength (shared by all genes)
	int geneLength = -1; //Undefined
	const string geneLengthString = input.getCmdOption("-gl");
	if (!geneLengthString.empty()){
		cout << "Gene length: " << geneLengthString << "\n";
		geneLength = atoi(geneLengthString.c_str());
	}
	
	//Get number of full haps
	int numHaps = -1;
	const string numHapsString = input.getCmdOption("-nh");
	if (!numHapsString.empty()){
		cout << "Number of full haplotypes input was: " << numHapsString << "\n";
		numHaps = atoi(numHapsString.c_str()); //Should be < 2^numLoci, e.g. 2^4=16 possible full haplotypes to choose from for four loci
	} else {
		numHaps = 8; //Default
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

	//Get simC value
	double simc = -1;
	const string simcString = input.getCmdOption("-simC");
	if (!simcString.empty()){
		cout << "simC input was: " << simcString << "\n";
		simc = atof(simcString.c_str());
	}
	//Get anaC value
	double anac = -1;
	const string anacString = input.getCmdOption("-anaC");
	if (!anacString.empty()){
		cout << "anaC input was: " << anacString << "\n";
		anac = atof(anacString.c_str());
	}  

	//Get BIC penalty value
	double BICpenalty = -1;
	const string BICString = input.getCmdOption("-BIC");
	if (!BICString.empty()){
		cout << "BIC penalty input was: " << BICString << "\n";
		BICpenalty = atof(BICString.c_str());
	} else {
		BICpenalty = 10; //Default
	}


	//Get selection cap
	double selCap = -1;
	const string selCapString = input.getCmdOption("-selCap");
	if(!selCapString.empty()) {
		cout << "Selection cap input was: " << selCapString << "\n";
		selCap = atof(selCapString.c_str());
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

	//Not relevant for realistic simulations	
	//Get sampling depth before
	int Nb = -1;
	const string NbString = input.getCmdOption("-Nb");
	if (!NbString.empty()){
		cout << "Nb input was: " << NbString << "\n";
		Nb = atoi(NbString.c_str());
	} else {
		Nb = 10000; //Default
	}

	//Get sampling depth after transmission
	int Na = -1;
	const string NaString = input.getCmdOption("-Na");
	if (!NaString.empty()){
		cout << "Na input was: " << NaString << "\n";
		Na = atoi(NaString.c_str());
	} else {
		Na = 10000; //Default
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

	//Get number of generations used for replications per 24 hour
	int numGenerations = 1; //Standard
	const string numGenerationsString = input.getCmdOption("-numGenerations");
	if (!numGenerationsString.empty()) {
		cout << "numGenerations input was: " << numGenerationsString << "\n";
		numGenerations = atoi(numGenerationsString.c_str());
	}

	//Get number of days between donor and recipients being sampled
	int deltaDays = 1; //Standard
	const string deltaDaysString = input.getCmdOption("-deltaDays");
	if (!deltaDaysString.empty()) {
		cout << "deltaDays input was: " << deltaDaysString << "\n";
		deltaDays = atoi(deltaDaysString.c_str());
	}


	//Get the growth factor used in simulation
	int growthFactorSim = 22; //Standard
	const string growthFactorSimString = input.getCmdOption("-gfSim");
	if (!growthFactorSimString.empty()) {
		cout << "growthFactorSim input was: " << growthFactorSimString << "\n";
		growthFactorSim = atoi(growthFactorSimString.c_str());
	}

	//Get the growth factor used in analysis
	int growthFactorAna = growthFactorSim; //Standard has sim and ana the same
	const string growthFactorAnaString = input.getCmdOption("-gfAna");
	if (!growthFactorAnaString.empty()) {
		cout << "growthFactorAna input was: " << growthFactorAnaString << "\n";
		growthFactorAna = atoi(growthFactorAnaString.c_str());
	}

	//Filtering of simulated data
	const string filterDataSimString = input.getCmdOption("-filterDataSim");
	string filterDataSim = "Transmission"; //Default
	if(!filterDataSimString.empty()) {
		cout << "filterDataSimString input was: " << filterDataSimString << "\n";
		if(filterDataSimString.compare("Transmission")==0 || filterDataSimString.compare("transmission")==0) {
			filterDataSim = "Transmission";
		} else if(filterDataSimString.compare("Samfire")==0 || filterDataSimString.compare("samfire")==0 || filterDataSimString.compare("SAMFIRE")==0) {
			filterDataSim = "Samfire";
		} else if(filterDataSimString.compare("None")==0 || filterDataSimString.compare("none")==0) {
			filterDataSim = "None";
		} else {

			cout << "Error with filterDataSim input: Input has to be 'Transmission' or 'Samfire'.\n"; exit(1);
		}	
	}


	//Get simOnly
	bool simOnly = false;
	const string simOnlyString = input.getCmdOption("-simOnly");
	if(!simOnlyString.empty()) {
		cout << "simOnly input was: " << simOnlyString << "\n";
		if(simOnlyString.compare("true")==0 || simOnlyString.compare("True")==0) {
			simOnly = true;
		} else if (simOnlyString.compare("false")==0 || simOnlyString.compare("False")==0) {
			simOnly = false;
		} else {

			cout << "Error with simOnly input: Input has to be true or false.\n"; exit(1);
		}	
	} else {
		simOnly = false; //Default
	}

	//Get termination point
	int terminationPoint = 0; //Default is not to terminate
	const string terminationString = input.getCmdOption("-term");
	if (!terminationString.empty()){
		cout << "Termination point input was: " << terminationString << "\n";
		terminationPoint = atoi(terminationString.c_str());
	}

	//Filter haps within the C++ code
	const string filterString = input.getCmdOption("-filter");
	int filterHapsMethod=0; //Default is 0, i.e. no filtering
	if(!filterString.empty()) {

		filterHapsMethod = atoi(filterString.c_str());
	}

	//Get analyseSL. If true, the program uses single locus methods (alongside multi locus methods) to infer transmission bottleneck
	bool analyseSL = false; //As default, don't use SL methods for analysing data
	const string analyseSLString = input.getCmdOption("-analyseSL");
	if(!analyseSLString.empty()) {
		cout << "analyseSL input was: " << analyseSLString << "\n";
		if(analyseSLString.compare("true")==0 || analyseSLString.compare("True")==0) {
			analyseSL = true;
		} else if (analyseSLString.compare("false")==0 || analyseSLString.compare("False")==0) {
			analyseSL = false;
		} else {

			cout << "Error with analyseSL input: Input has to be true or false.\n"; exit(1);
		}	
	}
	

	//Get export format, e.g. "Mahan" or "Haps"
	string format = "";
	const string formatString = input.getCmdOption("-format");
	if(!formatString.empty()) {
		cout << "Format input was: " << formatString << "\n";
		format = formatString; 
	} else {
		format = ""; //Default, undefined.
	}

	//Get selection information in the form "--A--_0_0_2_0_0_-----_0_0_0_0_0" etc
	const string selInfoString = input.getCmdOption("-sel");
	if (!selInfoString.empty()){
		cout << "Selection information input was: " << selInfoString << "\n";
	}
	vector<string> 	selInfoStringVec = split(selInfoString, '_'); //Split string by underscore
	
	vector<Sequence> selVec;
	vector<vector<double> > selMagVec;
	int currentGene = 0;
	int index = 0;
	while(currentGene < numGenesWithData) {
		
		Sequence currentSeq = Sequence(selInfoStringVec[index].c_str()); //Get the composition
		cout << "Current seq: "; currentSeq.print();
		selVec.push_back(currentSeq);
		index++;

		vector<double> selMag; //Get selection coefficients
		for(int i=0; i<currentSeq.getLength(); i++) { //Loop over number of loci
			selMag.push_back(atof(selInfoStringVec[index].c_str()));
			index++;
		}

		cout << "selMag:\n"; printDoubleVector(selMag);
		selMagVec.push_back(selMag);

		currentGene++;
	}

  

	//For simulation:
	//Get growth selection information in the form "--A--_0_0_2_0_0_-----_0_0_0_0_0" etc
	//Currently, these coefficients needs to have magnitudes of 24hrVal/numGenerationsPer24Hours. As standard, numGenerationsPer24hours=1.
	const string selGInfoString = input.getCmdOption("-selG");
	if (!selGInfoString.empty()){
		cout << "Selection for growth input was: " << selGInfoString << "\n";
		vector<string> 	selGInfoStringVec = split(selGInfoString, '_'); //Split string by underscore
	
		vector<Sequence> selGVec;
		vector<vector<double> > selGMagVec;
		currentGene = 0;
		index = 0;
		while(currentGene < numGenesWithData) {
		
			Sequence currentSeqG = Sequence(selGInfoStringVec[index].c_str()); //Get the composition
			//cout << "Current seqG: "; currentSeqG.print();
			selGVec.push_back(currentSeqG);
			index++;

			vector<double> selGMag; //Get selection coefficients
			for(int i=0; i<currentSeqG.getLength(); i++) { //Loop over number of loci
				selGMag.push_back(atof(selGInfoStringVec[index].c_str()));
				index++;
			}

			//cout << "selGMag:\n"; printDoubleVector(selGMag);
			selGMagVec.push_back(selGMag);

			currentGene++;
		}

		spPHMG.setWithinHostSelectionPresent(true);
		spPHMG.setSelMultiG(selGVec);
		spPHMG.setSelMagVecMultiG(selGMagVec);


		//If epistasis simulations are needed, make sure to implement.
		//Could maybe have an epistasis flag, e.g. -epiG 3_1_3_4_2.5 indicating 3 positions (1,3,4) with coefficient 2.5
		//Alternatively, for multiple epistasis coefficients, e.g. -epiG 2_3_1_3_4_2.5_2_1_3_1.2 indicating 2 coefficients
		//with 3 (1,3,4) and 2 (1,3) positions respectively and coefficients 2.5 and 1.2 respectively.
	} else {
		spPHMG.setWithinHostSelectionPresent(false);
	}
	

	//For analysis:
	//Get growth selection information in the form "--A--_0_0_2_0_0_-----_0_0_0_0_0" etc
	const string selGAnaInfoString = input.getCmdOption("-selGAna");
	if (!selGAnaInfoString.empty()){
		cout << "Selection for growth input (analysis) was: " << selGAnaInfoString << "\n";
		vector<string> 	selGAnaInfoStringVec = split(selGAnaInfoString, '_'); //Split string by underscore
	
		vector<Sequence> selGAnaVec;
		vector<vector<double> > selGAnaMagVec;
		currentGene = 0;
		index = 0;
		while(currentGene < numGenesWithData) {
		
			Sequence currentSeqG = Sequence(selGAnaInfoStringVec[index].c_str()); //Get the composition
			//cout << "Current seqG: "; currentSeqG.print();
			bool selPresent = false;
			for(int i=0; i<currentSeqG.getLength(); i++) {
				if(currentSeqG.getBase(i) != '-') { selPresent = true; break; }
			}
			ap.addWithinHostSelectionPresent(selPresent);
			selGAnaVec.push_back(currentSeqG);
			index++;

			vector<double> selGAnaMag; //Get selection coefficients
			for(int i=0; i<currentSeqG.getLength(); i++) { //Loop over number of loci
				selGAnaMag.push_back(atof(selGAnaInfoStringVec[index].c_str()));
				index++;
			}

			//cout << "selGAnaMag:\n"; printDoubleVector(selGAnaMag);
			selGAnaMagVec.push_back(selGAnaMag);

			currentGene++;
		}

		ap.setSelGvec(selGAnaVec);
		ap.setSelGmagVec(selGAnaMagVec);

		//If epistasis simulations are needed, make sure to implement.
		//Could maybe have an epistasis flag, e.g. -epiGAna 3_1_3_4_2.5 indicating 3 positions (1,3,4) with coefficient 2.5
		//Alternatively, for multiple epistasis coefficients, e.g. -epiGAna 2_3_1_3_4_2.5_2_1_3_1.2 indicating 2 coefficients
		//with 3 (1,3,4) and 2 (1,3) positions respectively and coefficients 2.5 and 1.2 respectively.

	}

	//Create output folder
	stringstream ssOutput;
	ssOutput << "/rds-d3/user/ckl35/hpc-work/Output/Simulated_Data/Analysis/";
	const string outputFolderString = input.getCmdOption("-o");
	if(!outputFolderString.empty()) {

		
		//Create relevant output folders, i.e. split tempFolders string by "/"
		stringstream ssTempFolder;
		if(anaSeedSpecified == true) { //AnaSeed specified so use this to denote output

			ssTempFolder << outputFolderString << "/Nt_" << Nt << "/Seed_" << anaSeed << "/";

		} else { //AnaSeed not specified, so use simSeed to denite output

			ssTempFolder << outputFolderString << "/Nt_" << Nt << "/Seed_" << simSeed << "/";
		}
		string folder;
		while(getline(ssTempFolder, folder, '/')) {

			ssOutput << folder << "/";
			string outputFolder = ssOutput.str();
			//If output folder doesn't exist, create it
			if(fileExists(outputFolder)!= true) {
				cout << "Ouput folder: " << outputFolder << " created.\n";
				mkdir(outputFolder.c_str(),0775);		
			}
		}

	} else { //Create output folder based on input - caution: Filename may get too long which can lead to issues

		//Create relevant output folders, i.e. split tempFolders string by "/"
		stringstream ssTempFolder;
		
		for(unsigned int i=0; i<selVec.size(); i++) {

			string selString = selVec[i].getSeqAsString();
			cout << "selString: " << selString << "\n";
			ssTempFolder << selString << "_";
			for(unsigned int j=0; j<selMagVec[i].size(); j++) {

				ssTempFolder << selMagVec[i][j] << "_";
			}
		}
		if(anaSeedSpecified == true) { //AnaSeed specified so use this to denote output
		
			ssTempFolder << "/Nt_" << Nt << "/Seed_" << anaSeed << "/";
		
		} else { //AnaSeed not specified, so use simSeed to denite output

			ssTempFolder << "/Nt_" << Nt << "/Seed_" << simSeed << "/";
		}

		
		string folder;
		while(getline(ssTempFolder, folder, '/')) {

			ssOutput << folder << "/";
			string outputFolder = ssOutput.str();
			//If output folder doesn't exist, create it
			if(fileExists(outputFolder)!= true) {
				cout << "Ouput folder: " << outputFolder << " created.\n";
				mkdir(outputFolder.c_str(),0775);		
			}
		}

	}
	string outputFolder = ssOutput.str();	

	//Related to analysis
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

	//Related to analysis
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




	//Create file to write input parameters to
	stringstream ssInputParams;
	ssInputParams << outputFolder << "InputParams.out";
	string filePathInputParams = ssInputParams.str();
	ofstream outputFileInputParams;
	outputFileInputParams.open(filePathInputParams.c_str());
	outputFileInputParams << "Parameters used for simulation/analysis:\n";
	outputFileInputParams << "SimMode: " << simMode << "\n";
	outputFileInputParams << "Nt: " << Nt << "\n";
	outputFileInputParams << "SimSeed: " << simSeed << "\n";
	outputFileInputParams << "AnaSeed: " << anaSeed << "\n";
	outputFileInputParams << "NumGenes: " << numGenes << "\n";
	outputFileInputParams << "NumGenesWithData: " << numGenesWithData << "\n";
	if(geneLength != -1) { outputFileInputParams << "GeneLength: " << geneLength << "\n"; }
	outputFileInputParams << "NumHaps: " << numHaps << "\n";
	outputFileInputParams << "C: " << c << "\n";
	if(selCap != -1) { outputFileInputParams << "SelCap: " << selCap << "\n"; }
	outputFileInputParams << "BIC penalty: " << BICpenalty << "\n";
	outputFileInputParams << "Max number of fitted params: " << maxNumFittedParams << "\n";
	outputFileInputParams << "Nb: " << Nb << "\n";
	outputFileInputParams << "Na: " << Na << "\n";
	outputFileInputParams << "maxNt (analysis): "  << maxNt << "\n";
	outputFileInputParams << "numGenerations: "  << numGenerations << "\n";
	outputFileInputParams << "deltaDays: "  << deltaDays << "\n";
	outputFileInputParams << "growthFactorSim: "  << growthFactorSim << "\n";
	outputFileInputParams << "growthFactorAna: "  << growthFactorAna << "\n";
	outputFileInputParams << "filterDataSim: " << filterDataSim << "\n";
	outputFileInputParams << "simOnly: " << simOnly << "\n";
	outputFileInputParams << "term: " << terminationPoint << "\n";
	outputFileInputParams << "filterHapsMethod: " << filterHapsMethod << "\n";
	outputFileInputParams << "analyseSL: " << analyseSL << "\n";
	outputFileInputParams << "Format: " << format << "\n";
	outputFileInputParams << "SelInfoString: " << selInfoString << "\n";
	outputFileInputParams << "SelGInfoString: " << selGInfoString << "\n";
	outputFileInputParams << "SelGAnaInfoString: " << selGAnaInfoString << "\n";
	outputFileInputParams << "OutputFolderString: " << outputFolderString << "\n";
	outputFileInputParams << "noVarString: " << noVarString << "\n";
	outputFileInputParams << "meanOnlyString: " << meanOnlyString << "\n";
	outputFileInputParams.close();
	

	
	//Set up simParam and simulate data
	spPHMG.setNb(Nb);
	spPHMG.setNa(Na);
	if(simc!=-1) { //Use simulation specific C value
		spPHMG.setC(simc);
	} else { //simC not specifically defined, so use joint C for sim and ana
		spPHMG.setC(c);
	}
	spPHMG.setNt(Nt);
	spPHMG.setNumGenes(numGenes);
	spPHMG.setNumGenesWithData(numGenesWithData);
	if(geneLength != -1) { spPHMG.setGeneLength(geneLength); }
	spPHMG.setSeed(simSeed);
	spPHMG.setSelMulti(selVec);
	spPHMG.setSelMagVecMulti(selMagVec);
	spPHMG.setDim(numHaps);
	spPHMG.setNumGenerations(numGenerations);
	spPHMG.setDeltaDays(deltaDays);
	spPHMG.setGrowthFactor(growthFactorSim);
	spPHMG.setFormat(format); //Define output format, e.g. "Mahan"
	spPHMG.setPathToFolder(outputFolder);
	spPHMG.setFilterData(filterDataSim);
	DataPHMG dphmg;
	if(simMode == 0) {
        	dphmg.simulateDataRealistic(&spPHMG);
	} else if(simMode == 1) {
		dphmg.simulateData(&spPHMG);
	}
        

	if(simOnly == false) { //Also analyse data
	
		cout << "Preparing for analysis.\n";

		//Set up anaparam and analyse data
		if(anac!=-1) { //Use analysis specific C value
			ap.setC(anac);
		} else { //anaC not specifically defined, so use joint C for sim and ana
			ap.setC(c);
		}
		ap.setBICpenalty(BICpenalty);
		ap.setMaxNumFittedParams(maxNumFittedParams);
		ap.setSeed(anaSeed);
		ap.setMaxNt(maxNt);
		ap.setGrowthFactor(growthFactorAna);
		if(selCap != -1) {
			ap.setSelectionCap(selCap);
		}
		ap.setAnalyseSL(analyseSL);
		ap.setFilterHaplotypesMethod(filterHapsMethod);
		ap.setTerminateEarly(terminationPoint);
        	ap.loadOutputFolder(outputFolder);
		ap.setNoVar(noVar);
		ap.setMeanOnly(meanOnly);	
		AnalyserPHMG aPHMG;
        	aPHMG.loadData(&dphmg);
	        aPHMG.runAnalysis(&ap);

	} else {

		cout << "Finishing without carrying out analysis.\n";
	}
	
	//Print elapsed time to file 
	int timeFinish = omp_get_wtime();
	cout << "Elapsed time: " << timeFinish - timeStart << "\n";
	
	stringstream ssTime;
	ssTime << ssOutput.str() << "Elapsed_Time.dat";
	string filePathTime = ssTime.str();
	ofstream outputFileTime;
	outputFileTime.open(filePathTime.c_str());
	outputFileTime << timeFinish - timeStart << "\n";
	outputFileTime.close();
  
	return 0;
}

