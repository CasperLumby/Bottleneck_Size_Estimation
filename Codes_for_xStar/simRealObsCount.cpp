#include "dataPHMG.hpp"
#include "dataPHMR.hpp"
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



/*
 * Do we need this?
 */
	//Get number of genes
	int numGenes = -1;
	const string numGenesString = input.getCmdOption("-ng");
	if (!numGenesString.empty()){
		cout << "Number of genes input was: " << numGenesString << "\n";
		numGenes = atoi(numGenesString.c_str());
	} else {
		numGenes = 8; //Default
		cout << "numGenes set to default: " << numGenes << "\n";
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
		cout << "numHaps set to default: " << numHaps << "\n";
	}

	//Get C value
	double c = -1;
	const string cString = input.getCmdOption("-C");
	if (!cString.empty()){
		cout << "C input was: " << cString << "\n";
		c = atof(cString.c_str());
	} else {
		c = 200; //Default
		cout << "C set to default: " << c << "\n";
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

	//Get dataset folder
	//(defined as folders after /rds-d3/user/ckl35/hpc-work/Datasets/)
	//Disregard final "/"!
	bool inputGiven = false;
	string inputFolder;
	stringstream ssInput;
	ssInput << "/rds-d3/user/ckl35/hpc-work/Datasets/";
	const string dataset = input.getCmdOption("-i");
	if(!dataset.empty()) {

		//Add dataset to input
		ssInput << dataset << "/";
		inputFolder = ssInput.str();
		inputGiven = true;

	}	


	//Get dataset folder in full form, e.g.
	// /rds-d3/user/ckl35/hpc-work/Datasets/Moncla/FolderA/FolderB/FolderC
	//Disregard final "/"!
	stringstream ssInputF;
	const string datasetF = input.getCmdOption("-if");
	if(!datasetF.empty()) {

		//Add dataset to input
		ssInputF << datasetF << "/";
		inputFolder = ssInputF.str();
		inputGiven = true;

	}

	if(inputGiven == false) {

		cout << "Error in getting the correct filepath!\n";
		return -1;
	}

	//Flag to use imported haplotypes and frequencies from Mahan format.
	//Reads in one set of full haplotypes and two sets of frequencies
	//(before and after transmission). String to folder containing
	//haplotypes and frequencies must be supplied. The folder must contain
	//8 subfolders, "/test_0/", ... , "/test_7/". Each subfolder must contain
	//a file named "outcome_1.txt" which is of the form:
	//	1	h1	qB1	qA1
	//	2	h2	qB2	qA2
	//	..
	//	k	hk	qBk	qAk
	//
	//where each column is tab delimited.
	const string importHapsMahanString = input.getCmdOption("-importHapsMahan");
	bool importHapsMahan = false;
	string importHapsMahanPath = "";
	if(!importHapsMahanString.empty()) {

		importHapsMahan = true;
		stringstream ss;
		ss << importHapsMahanString << "/"; //Add '/' in case folder doesn't end in '/' itself
		importHapsMahanPath = ss.str();
	} else {

		importHapsMahan = false; //Default, don't use imported haplotypes
	}

	//Flag to use a different filename than "outcome_1.txt" for Mahan's haplotype information
	const string importHapsMahanFileString = input.getCmdOption("-importHapsMahanFile");
	string importHapsMahanFilename = "";
	if(!importHapsMahanFileString.empty()) {

		importHapsMahanFilename = importHapsMahanFileString;

	} else {
		importHapsMahanFilename = "outcome_1.txt"; //Standard filename
	}



	//Create output folder - defined as folders after "/rds-d3/user/ckl35/hpc-work/Output/"
	string outputFolder;
	stringstream ssOutput;
	ssOutput << "/rds-d3/user/ckl35/hpc-work/Output/";
	const string outputFolderString = input.getCmdOption("-o");
	if(!outputFolderString.empty()) {

		
		stringstream ssTempFolder;
		if(importHapsMahan == true) { //No transmission in this case, so no Nt
			ssTempFolder << outputFolderString << "/Seed_" << simSeed << "/";
		} else { //Default
			ssTempFolder << outputFolderString << "/Nt_" << Nt << "/Seed_" << simSeed << "/";
		}
		
		//Create relevant output folders, i.e. split tempFolders string by "/"
		string folder;
		while(getline(ssTempFolder, folder, '/')) {

			ssOutput << folder << "/";
			string outputFolderCurrent = ssOutput.str();
			//If output folder doesn't exist, create it
			if(fileExists(outputFolderCurrent)!= true) {
				cout << "Ouput folder: " << outputFolderCurrent << " created.\n";
				mkdir(outputFolderCurrent.c_str(),0775);		
			}
		}

		outputFolder = ssOutput.str();

	}


	//Create output folder based on full path
	//e.g. /rds-d3/user/ckl35/hpc-work/Output/Moncla/Analysis/FolderA/FolderB/FolderC
	stringstream ssOutputF;
	const string outputFolderStringF = input.getCmdOption("-of");
	if(!outputFolderStringF.empty()) {

		
		//Create relevant output folders, i.e. split tempFolders string by "/"
		stringstream ssTempFolder;
		if(importHapsMahan == true) { //No transmission in this case, so no Nt
			ssTempFolder << outputFolderStringF << "/Seed_" << simSeed << "/";
		} else { //Default
			ssTempFolder << outputFolderStringF << "/Nt_" << Nt << "/Seed_" << simSeed << "/";
		}
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

	//Filter data (read in data, not simulated data) (not haplotypes) - this is a data flag, not an analysis flag
	const string filterDataString = input.getCmdOption("-filterData");
	bool filterData = true; //Default
	if(!filterDataString.empty()) {
		cout << "filterDataString input was: " << filterDataString << "\n";
		if(filterDataString.compare("true")==0 || filterDataString.compare("True")==0) {
			filterData = true;
		} else if (filterDataString.compare("false")==0 || filterDataString.compare("False")==0) {
			filterData = false;
		} else {

			cout << "Error with filterData input: Input has to be true or false.\n"; exit(1);
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

	//Get export format, e.g. "Mahan" or "Haps"
	string format = "";
	const string formatString = input.getCmdOption("-format");
	if(!formatString.empty()) {
		cout << "Format input was: " << formatString << "\n";
		format = formatString; 
	} else {
		format = ""; //Default, undefined.
	}

	//Get remove monomorphic sites in simulated data
	bool removeMonomorphicSim = true;
	const string removeMonomorphicSimString = input.getCmdOption("-rmMonoSim");
	if(!removeMonomorphicSimString.empty()) {
		cout << "removeMonomorphicSim input was: " << removeMonomorphicSimString << "\n";
		if(removeMonomorphicSimString.compare("true")==0 || removeMonomorphicSimString.compare("True")==0) {
			removeMonomorphicSim = true;
		} else if (removeMonomorphicSimString.compare("false")==0 || removeMonomorphicSimString.compare("False")==0) {
			removeMonomorphicSim = false;
		} else {

			cout << "Error with removeMonomorphicSim input: Input has to be true or false.\n"; exit(1);
		}	
	} else {
		removeMonomorphicSim = true; //Default
	}

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
	if(geneLength != -1) { outputFileInputParams << "GeneLength: " << geneLength << "\n"; }
	outputFileInputParams << "NumHaps: " << numHaps << "\n";
	outputFileInputParams << "C: " << c << "\n";
	outputFileInputParams << "BIC penalty: " << BICpenalty << "\n";
	outputFileInputParams << "Max number of fitted params: " << maxNumFittedParams << "\n";
	outputFileInputParams << "filterDataSim: " << filterDataSim << "\n";
	outputFileInputParams << "OutputFolderString: " << outputFolderString << "\n";
	outputFileInputParams << "Reps: " << reps << "\n";
	outputFileInputParams << "ImportHaps: " << importHaps << "\n";
	if(importHaps == true) { outputFileInputParams << "ImportHapsFreq: " << importHapsFreq << "\n"; }
	if(importHapsMahan == true) { outputFileInputParams << "ImportHapsMahanString: " << importHapsMahanString << "\n"; }
	if(filterHapsMethod > 0) { outputFileInputParams << "Filter haplotypes method: " << filterHapsMethod << "\n"; }
	outputFileInputParams << "simOnly: " << simOnly << "\n";
	outputFileInputParams << "Format: " << format << "\n";
	outputFileInputParams << "removeMonomorphicSim: " << removeMonomorphicSim << "\n";
	outputFileInputParams.close();
	
	//Generate a sel model for the data
	vector<Sequence> selVec;
	vector<vector<double> > selMagVec;


	//Get selection information
	int selGene;
	const string selGeneString = input.getCmdOption("-selGene");
	if (!selGeneString.empty()){
		cout << "SelGene input was: " << selGeneString << "\n";
		selGene = atoi(selGeneString.c_str());
	} else {
		selGene = -1; //Neutral
	}

	int selLocus;
	const string selLocusString = input.getCmdOption("-selLocus");
	if (!selLocusString.empty()){
		cout << "SelLocus input was: " << selLocusString << "\n";
		selLocus = atoi(selLocusString.c_str());
	} else {
		selLocus = -1; //Neutral
	}

	double selStrength;
	const string selStrengthString = input.getCmdOption("-selStrength");
	if (!selStrengthString.empty()){
		cout << "selStrength input was: " << selStrengthString << "\n";
		selStrength = atof(selStrengthString.c_str());
	} else {
		selStrength = 0; //Neutral
	}

	//Read in external data to get observation counts
	DataPHMG dPHMG;
	DataPHMR dPHMR;
	if(repsPresent == false) {
		PathPHMG pPHMG;
		if(importHaps == true) {
			pPHMG.setImportHapsFreq(importHapsFreq);
		}
		if(importHapsMahan == true) { //Perhaps do else if
			pPHMG.setImportHapsMahan(importHapsMahanPath);
			pPHMG.setImportHapsMahanFilename(importHapsMahanFilename);
		}
		pPHMG.setFilterData(filterData);
		pPHMG.loadFluFilePaths(inputFolder);
		dPHMG.readData(&pPHMG);

		//If imported haps are used, we need to correct the data for these
		//Discussion point: We definitely need to do this when analysing real data based on imported haplotypes,
		//but what about the case where we generate simulated data based on real data and imported haplotypes?
		//Surely here, we just need the total read depths for each partial haplotype set + the loci the cover, i.e.
		//we do not care about what the actual partial haplotypes look like, nor if they are represented by an
		//imported full haplotype. Why? Because for generating data we simply generate data based on the imported
		//haplotypes and the loci covered by the partial haplotype set in question. As such, presumably, here we
		//shouldn't need the below?
		if(importHaps == true) {

			cout << "Imported haps triggered, so updating data to account for this:\n";
			for(int g=0; g<dPHMG.getNumGenes(); g++) {
			
				if(dPHMG.getGenePresent(g) == true) {

					DataPH* dPH = dPHMG.getGene(g);
					cout << "Printing data for gene " << g << ":\n";
					dPH->print();

					cout << "Printing imported full haps:\n";
					vector<Sequence> importedFullHaps = dPH->getImportedFullHaps();
					for(unsigned int i=0; i<importedFullHaps.size(); i++) {
						importedFullHaps[i].print();
					}
			
					//Need to do further processing here
					//Essentially, two things can happen: 
					//1) Some loci in the full haplotypes are monomorphic. This means that this locus needs 
					//deleting from all partial haps. Need to update physical pos as well, if present.
					//2) Some partial haps may no longer be represented by a full haplotype. In this case we
					//need to shorten the haplotype until it can be merged with a shorter haplotype.

					if(importedFullHaps.size() > 0) { //Only update data if there actually are imported full haplotypes 
				
						dPH->updateDataWithRespectToImportedFullHaps(); //Also computes contribs.

						//The above might move partial haplotypes from one set to a shorter subset.
						//In this process, some sophs become illegal. This is the case if NaTot = 0,
						//which is not defined for our inference process.
						//To this end, we rerun the "removeLowCountObservations" method. This takes care of that.
						//It also updates lociCovered and physicalPos.
						//Note that removing possible monomorphic SNPs should not be necessary, as the imported
						//full haplotypes should already be polymorphic.
						dPH->removeLowCountObservations();

					} else { //Remove all data, as no full haplotypes to match

						cout << "No imported full haps, so removing all partial haplotype for gene " << g << " if present.\n";
						dPH->clear();
						dPH->clearPhysicalPos(); //If no data, can't have physical pos either
						dPHMG.setGenePresent(g,false);	
					}

					//Physical positions is implemented somewhat badly. In particular, each data entity has its own copy, but the upper
					//levels do not make use of the lower levels, e.g. if the lower levels are changed, the upper levels don't.
					//At some future restructuring of the code, this should be attended to. For the time being, here is a simple function
					//that copies data from lower structures up to higher ones.
					if(dPHMG.physicalPosPresent()==true) {
						dPHMG.updatePhysicalPos(g, dPH);
					}


				} 
			}


//			cout << "Printing data after modification by imported haps:\n";
//			for(int g=0; g<dPHMG.getNumGenes(); g++) {
//			
//				DataPH* dPH = dPHMG.getGene(g);
//				cout << "Printing haplotypes for gene " << g << ":\n";
//				vector<Sequence> modifiedFullHaps = dPH->getImportedFullHaps();
//				for(unsigned int i=0; i<modifiedFullHaps.size(); i++) {
//					modifiedFullHaps[i].print();
//				}
//				cout << "Printing data for gene " << g << ":\n";
//				dPH->print();
//				cout << "\n\n";
//			}
		}


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
		dPHMR.readData(&pPHMR);

		//If imported haps are used, we need to correct the data for these
		if(importHaps == true) {

			cout << "Imported haps triggered, so updating data to account for this:\n";

			for(int r=0; r<dPHMR.getNumReplicates(); r++) {

				DataPHMG* dPHMG = dPHMR.getReplicate(r);
				for(int g=0; g<dPHMG->getNumGenes(); g++) {
			
					if(dPHMG->getGenePresent(g) == true) {

						DataPH* dPH = dPHMG->getGene(g);
						cout << "Printing data for gene " << g << ":\n";
						dPH->print();

						cout << "Printing imported full haps:\n";
						vector<Sequence> importedFullHaps = dPH->getImportedFullHaps();
						for(unsigned int i=0; i<importedFullHaps.size(); i++) {
							importedFullHaps[i].print();
						}
			
						//Need to do further processing here
						//Essentially, two things can happen: 
						//1) Some loci in the full haplotypes are monomorphic. This means that this locus needs 
						//deleting from all partial haps. Need to update physical pos as well, if present.
						//2) Some partial haps may no longer be represented by a full haplotype. In this case we
						//need to shorten the haplotype until it can be merged with a shorter haplotype.

						if(importedFullHaps.size() > 0) { //Only update data if there actually are imported full haplotypes 
					
							dPH->updateDataWithRespectToImportedFullHaps(); //Also computes contribs.
	
							//The above might move partial haplotypes from one set to a shorter subset.
							//In this process, some sophs become illegal. This is the case if NaTot = 0,
							//which is not defined for our inference process.
							//To this end, we rerun the "removeLowCountObservations" method. This takes care of that.
							//It also updates lociCovered and physicalPos.
							//Note that removing possible monomorphic SNPs should not be necessary, as the imported
							//full haplotypes should already be polymorphic.
							dPH->removeLowCountObservations();

						} else { //Remove all data, as no full haplotypes to match

							cout << "No imported full haps, so removing all partial haplotype for gene " << g << " if present.\n";
							dPH->clear();
							dPH->clearPhysicalPos(); //If no data, can't have physical pos either
							dPHMG->setGenePresent(g,false);	
						}

						//Physical positions is implemented somewhat badly. In particular, each data entity has its own copy, but the upper
						//levels do not make use of the lower levels, e.g. if the lower levels are changed, the upper levels don't.
						//At some future restructuring of the code, this should be attended to. For the time being, here is a simple function
						//that copies data from lower structures up to higher ones.
						if(dPHMG->physicalPosPresent()==true) {
							dPHMG->updatePhysicalPos(g, dPH);
						}


					} 
				}

				//Physical positions are implemented somewhat badly. See description above.
				if(dPHMR.physicalPosPresent()==true) {
					dPHMR.updatePhysicalPos(r,dPHMG);
				}
			}
		}


	}

	
	if(repsPresent == false) {
	
		//Loop over genes
		for(int i=0; i<dPHMG.getNumGenes(); i++) {

			cout << "Gene " << i << " status: " << dPHMG.getGenePresent(i) << "\n";
			if(dPHMG.getGenePresent(i) == true) {


				//Find the haplotype length for gene i
				DataPH* dPH = dPHMG.getGene(i);
				vector<Sequence> allSeqs = dPH->getSequences();
				int hapLength = allSeqs[0].getLength();
				cout << "Hap length: " << hapLength << "\n";

				Sequence selSeq = Sequence(hapLength); //-,-,-,-,-,-
				vector<double> selMag (hapLength, 0.0); //0,0,0,0,0,0
			
				if(i==selGene) { //Gene under selection
			
				
					//Find allele for locus under selection
					char allele = '%';
					for(unsigned int j=0; j<allSeqs.size(); j++) {
			
						if(allSeqs[j].getBase(selLocus) != '-') {

							allele = allSeqs[j].getBase(selLocus); //WLOG choose the first allele found
							cout << "Selection chosen to act on " << allele << " at locus " << selLocus << " with strength " << selStrength << ".\n";
							break;
						}
					}
	
					selSeq.setBase(selLocus, allele); //Set selection act on A's at the mid point of gene 0
					selMag[selLocus] = selStrength;
				}

				selVec.push_back(selSeq);
				selMagVec.push_back(selMag);

			} else {
		
				Sequence emptySeq;
				vector<double> emptyVec;
				selVec.push_back(emptySeq);
				selMagVec.push_back(emptyVec);
			}
		}

	} else {	
	
		vector<vector<int> > allPos = dPHMR.findAllPos();
		vector<vector<vector<int> > > mapFromAllPosToReps = dPHMR.findMapFromAllPosToReps(allPos);
		vector<vector<bool> > mapFromAllPosToShared = dPHMR.findMapFromAllPosToShared(allPos, mapFromAllPosToReps);

		
		//Loop over genes
		for(unsigned int i=0; i<allPos.size(); i++) {

			if(allPos[i].size() > 0) { //Gene has positions

				Sequence selSeq = Sequence(allPos[i].size()); //-,-,-,-,-,-
				vector<double> selMag (allPos[i].size(), 0.0); //0,0,0,0,0,0
			
				if((int) i==selGene) { //Gene under selection
			
					//Find shared selLocus
					int sharedSelLocus = -1;
					int sharedPosCounter = 0;
					for(unsigned int j=0; j<allPos[i].size(); j++) {

						if(mapFromAllPosToShared[i][j] == true) { //Pos is shared between all replicates
							sharedPosCounter++;
						}
						if(sharedPosCounter == selLocus+1) {

							sharedSelLocus = j;
							break;
						}
					}

					//Consider replicate 0 WLOG
					DataPHMG* dPHMG = dPHMR.getReplicate(0);
					DataPH* dPH = dPHMG->getGene(i);
					vector<Sequence> allSeqs = dPH->getSequences();
					
			
					//Find allele for locus under selection
					char allele = '%';
					for(unsigned int j=0; j<allSeqs.size(); j++) {
			
						if(allSeqs[j].getBase(mapFromAllPosToReps[i][0][sharedSelLocus]) != '-') {

							allele = allSeqs[j].getBase(mapFromAllPosToReps[i][0][sharedSelLocus]); //WLOG choose the first allele found
							cout << "Selection chosen to act on " << allele << " at locus " << sharedSelLocus << " with strength " << selStrength << ".\n";
							break;
						}
					}

					selSeq.setBase(sharedSelLocus, allele); 
					selMag[sharedSelLocus] = selStrength;
				}

				selVec.push_back(selSeq);
				selMagVec.push_back(selMag);

			} else {
		
				Sequence emptySeq;
				vector<double> emptyVec;
				selVec.push_back(emptySeq);
				selMagVec.push_back(emptyVec);
			}
		}	
	}



	
	//Set up simParam and simulate data using read statistic from external data set
	spPHMG.setC(c);
	spPHMG.setNt(Nt);
	spPHMG.setNumGenes(numGenes);
	if(geneLength != -1) { spPHMG.setGeneLength(geneLength); }
	spPHMG.setSeed(simSeed);
	spPHMG.setSelMulti(selVec);
	spPHMG.setSelMagVecMulti(selMagVec);
	spPHMG.setDim(numHaps);
	spPHMG.setPathToFolder(outputFolder);
	spPHMG.setUseExternalReadStatistic(true);
	spPHMG.setFilterHaplotypesMethod(filterHapsMethod);
	spPHMG.setFormat(format); //Define output format, e.g. "Mahan"
	spPHMG.setRemoveMonomorphicSim(removeMonomorphicSim);
	spPHMG.setFilterData(filterDataSim);
	

	if(repsPresent == false) {
		dPHMG.simulateData(&spPHMG);
	} else {
		dPHMR.simulateRealData(&spPHMG);
	} 
	

	if(simOnly == false) { //Also analyse data
	
		//Set up anaparam and analyse data
		ap.setC(c);
		ap.setBICpenalty(BICpenalty);
		ap.setMaxNumFittedParams(maxNumFittedParams);
		ap.setSeed(anaSeed);
		ap.loadOutputFolder(outputFolder);
		ap.setFilterHaplotypesMethod(filterHapsMethod);
		if(repsPresent == false) {
			AnalyserPHMG aPHMG;
			aPHMG.loadData(&dPHMG);
			aPHMG.runAnalysis(&ap);
		} else {

			AnalyserPHMR aPHMR;
			aPHMR.loadData(&dPHMR);
			aPHMR.runAnalysis(&ap);
	
		}
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

