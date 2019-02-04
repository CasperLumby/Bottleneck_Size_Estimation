#include "dataPHMR.hpp"
#include "simParam.h"
#include "anaParam.h"
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
#include <string.h>
#include "sys/types.h"
#include "stdlib.h"
#include "stdio.h"

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
	if (!anaSeedString.empty()){
		cout << "Analysis seed input was: " << anaSeedString << "\n";
		anaSeed = atoi(anaSeedString.c_str());
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




	//Get number of genes (most of the time, doesn't need changing)
	int numGenes = -1;
	const string numGenesString = input.getCmdOption("-ng");
	if (!numGenesString.empty()){
		cout << "Number of genes input was: " << numGenesString << "\n";
		numGenes = atoi(numGenesString.c_str());
	} else {
		numGenes = 8; //Default
	}
	
	//Get number of genes with data
	int numGenesWithData = -1;
	const string numGenesWithDataString = input.getCmdOption("-ngwd");
	if (!numGenesWithDataString.empty()){
		cout << "Number of genes with data input was: " << numGenesWithDataString << "\n";
		numGenesWithData = atoi(numGenesWithDataString.c_str());
	} else {
		numGenesWithData = 8; //Default
	}

	//Perhaps genLength input here?


	//Get number of loci per gene
	int numLociPerGene = -1;
	const string numLociPerGeneString = input.getCmdOption("-nlpg");
	if (!numLociPerGeneString.empty()){
		cout << "Number of loci per gene input was: " << numLociPerGeneString << "\n";
		numLociPerGene = atoi(numLociPerGeneString.c_str());
	} else {
		numLociPerGene = 5; //Default
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


	//Get number of replicates
	int numReps = -1;
	const string numRepsString = input.getCmdOption("-nr");
	if (!numRepsString.empty()){
		cout << "Number of replicate input was: " << numRepsString << "\n";
		numReps = atoi(numRepsString.c_str());
	} else {
		numReps = 3; //Default
	}
	
	//Get selection information in the form "--A--_0_0_2_0_0_-----_0_0_0_0_0" etc.
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


	//Perhaps SelG input here?
	

	//Perhaps SelGAna input here?

	//Get C value
	double c = -1;
	const string cString = input.getCmdOption("-C");
	if (!cString.empty()){
		cout << "C input was: " << cString << "\n";
		c = atof(cString.c_str());
	} else {
		c = 200; //Default
	}
	
	//Get information on using shared bottlenecks
	bool useSharedBottleneck = false;
	const string useSharedBottleneckString = input.getCmdOption("-usb");
	if (!useSharedBottleneckString.empty()){
		cout << "Use shared bottleneck input was: " << useSharedBottleneckString << "\n";
		string useSharedBottleneckCstring = useSharedBottleneckString.c_str();
		useSharedBottleneck = stringToBool(useSharedBottleneckCstring);
	} else {
		useSharedBottleneck = false; //Default
	}

	
	//Create output folder
	stringstream ssOutput;
	ssOutput << "/scratch/ckl35/Output/Simulated_Data/Analysis/";
	const string outputFolderString = input.getCmdOption("-o");
	if(!outputFolderString.empty()) {

		
		//Create relevant output folders, i.e. split tempFolders string by "/"
		stringstream ssTempFolder;
		ssTempFolder << outputFolderString << "/Nt_" << Nt << "/Seed_" << simSeed << "/";
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

	} else { //No alternative

		cout << "No output folder name given. Exiting.\n";
		exit(1);
	}
        string outputFolder = ssOutput.str();

	//Create file to write input parameters to
	stringstream ssInputParams;
	ssInputParams << outputFolder << "InputParams.out";
	string filePathInputParams = ssInputParams.str();
	ofstream outputFileInputParams;
	outputFileInputParams.open(filePathInputParams.c_str());
	outputFileInputParams << "Parameters used for simulation/analysis:\n";
	outputFileInputParams << "Nt: " << Nt << "\n";
	outputFileInputParams << "SimSeed: " << simSeed << "\n";
	outputFileInputParams << "AnaSeed: " << anaSeed << "\n";
	outputFileInputParams << "NumGenes: " << numGenes << "\n";
	outputFileInputParams << "NumGenesWithData: " << numGenesWithData << "\n";
//	if(geneLength != -1) { outputFileInputParams << "GeneLength: " << geneLength << "\n"; }
	outputFileInputParams << "NumHaps: " << numHaps << "\n";
	outputFileInputParams << "C: " << c << "\n";
	if(selCap != -1) { outputFileInputParams << "SelCap: " << selCap << "\n"; }
	outputFileInputParams << "BIC penalty: " << BICpenalty << "\n";
	outputFileInputParams << "Max number of fitted params: " << maxNumFittedParams << "\n";
	outputFileInputParams << "SelInfoString: " << selInfoString << "\n";
//	outputFileInputParams << "SelGInfoString: " << selGInfoString << "\n";
//	outputFileInputParams << "SelGAnaInfoString: " << selGAnaInfoString << "\n";
	outputFileInputParams << "OutputFolderString: " << outputFolderString << "\n";
	outputFileInputParams.close();
	





	SimParamPH spPHMR;
	spPHMR.setC(c);
	spPHMR.setNt(Nt);
	spPHMR.setNumGenes(numGenes);
	spPHMR.setNumGenesWithData(numGenesWithData); //Implicit that total number of genes is 8 - only need to specify how many we have data for
	spPHMR.setNumLociPerGene(numLociPerGene);
	spPHMR.setSeed(simSeed);
	spPHMR.setSelMulti(selVec);
	spPHMR.setSelMagVecMulti(selMagVec);
	spPHMR.setDim(numHaps);
	spPHMR.setPathToFolder(outputFolder);
	spPHMR.setNumReps(numReps);
        DataPHMR dphmr;
        dphmr.simulateData(&spPHMR);

	AnaParam ap;
	ap.setC(c);
	ap.setBICpenalty(BICpenalty);
	ap.setMaxNumFittedParams(maxNumFittedParams);
	ap.setSeed(anaSeed);
	if(selCap != -1) {
		ap.setSelectionCap(selCap);
	}
        ap.loadOutputFolder(outputFolder);
	ap.setUseSharedBottleneck(useSharedBottleneck);
        AnalyserPHMR aPHMR;
        aPHMR.loadData(&dphmr);
        aPHMR.runAnalysis(&ap);

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

