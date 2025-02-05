//
//  dataPH.cpp
//  FluTransmissionProject
//
//  Created by Casper Lumby on 16/11/2015.
//  Copyright © 2015 Casper Lumby. All rights reserved.
//

#include "dataPH.hpp"
#include <stdio.h>
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "misc.h"
#include "analyserPH.hpp"
#include <math.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <algorithm> //Apparently necessary for sort command
#include <sys/stat.h> //For folder management

using namespace std;

//Class related methods + constructors
DataPH::DataPH() : importFullHaps(false), importHapsMahan(false) { } //Constructor does nothing
DataPH::~DataPH() {} //Deconstructor

void DataPH::readData(Path * p) {

	//Load paths from path object
	path = static_cast<PathPH *>(p);
	vector<string> phPaths = path->getHapPaths();
	vector<string> lociPaths = path->getLociPaths();
	string slPath = path->getSingleLocusTrajFilePath();

	if(phPaths.size() > 0) { //Only try and add data if present. If no data, DataPHMG takes care of this.

		//Find the maximum loci number
		int lociMax = -1;
		vector<int> lociIncluded;
		vector<int> includedFiles;
		for(unsigned int i=0; i<lociPaths.size(); i++) {

			ifstream lociFile;
			lociFile.open(lociPaths[i].c_str());
			string loci;

			while(lociFile >> loci) { //Check if loci is the largest found so far

				int lociInt = atoi(loci.c_str());
				if(lociInt > lociMax) { lociMax = lociInt; }
			}
		}

		//Read in physical positions, only in cases where relevant
		//Assumes positions are ordered by size
		if(path->getPhysicalPosMatter() == true) {

			string physicalPosPath = path->getPhysicalPosPath();

			if(fileExists(physicalPosPath) == true) {
				ifstream posFile;
				posFile.open(physicalPosPath.c_str());
				string pos;

				while(posFile >> pos) { //Add pos to string and keep going until all pos have been imported

					int posInt = atoi(pos.c_str());
					physicalPos.push_back(posInt);

				}	

			} else {

				cout << "No physical pos file for current gene. Aborting.\n";
				exit(1);
			}
		}


		//Read in data
		for(unsigned int i=0;i<phPaths.size();i++) { //Loop over sets of partial haps


			setOfPartialHaps soph;

			ifstream phFile;
			ifstream lociFile;
			phFile.open(phPaths[i].c_str());
			cout << "Current path: " << phPaths[i] << "\n";
			string word1,word2,word3;

			//Import partial haplotype sequences, observations and contributions
			while(phFile >> word1) {

				//If word1 (column1) is an X, stop importing data
				if(word1.compare("X")==0) { break; }

				//Else add data
				partialHap ph;
				Sequence seq = Sequence(word1);  //Word1 correspond to the partial hap sequence

				//Correct for loci positions, i.e. add "-"s from front and back to get the sequence length to correspond to the FH
				//First from the front
				lociFile.open(lociPaths[i].c_str());
				string loci;
				lociFile >> loci;
				int lociInt = atoi(loci.c_str()); //This is the first loci of the partial haplotype
				for(int j=0;j<(lociMax+1);j++) {
					if(lociInt != j) {
						seq.addBaseFront('-');
					} else {
						break;
					}
				}
				//Then from the back
				while(lociFile >> loci) { }
				lociInt = atoi(loci.c_str()); //Last loci of the partial haplotype
				for(int j=lociMax; j>=0; j--) {
					if(lociInt != j) {
						seq.addBase('-');
					} else {
						break;
					}
				}
				lociFile.close();


				//Add the corected sequence to the partial haplotype
				ph.setSeq(seq);

				phFile >> word2; //word in column2, i.e. number of obs at time T1
				phFile >> word3; //word in column3, i.e. number of obs at time T2

				int obsBefore = atoi(word2.c_str());
				int obsAfter = atoi(word3.c_str());

				ph.setNb(obsBefore); //Obs before
				ph.setNa(obsAfter); //Obs after


				//Add haplotype to set
				soph.addPartialHap(ph);
			}


			//Count NbTot and NaTot for soph + update loci covered
			soph.countNbTot();
			soph.countNaTot();
			vector<int> lociCovered;
			Sequence ph0Seq = soph.getSeq(0);
			for(int l=0; l<ph0Seq.getLength(); l++) {
				if(ph0Seq.getBase(l) != '-') {

					lociCovered.push_back(l);
				}
			}
			cout << "lociCovered (temp print): "; printIntVector(lociCovered);
			soph.setLociCovered(lociCovered);


			//Add set of partial haplotypes
			data.addSetOfPartialHaps(soph);
		}

		if(path->getFilterData() == true) {

			//To mimick samfire, we remove partial haplotype observations with less than 10 observations at both time points
			removeLowCountObservations();


			//Fix any possible monomorphic sites arising from low count removal
			if(data.getNumDatasets() > 0 ) { //Only if any sophs remaining
				removeMonomorphicSites();
			}

			//Read in single locus read depths and remove loci for which read depth is too low
			ifstream single_trajs_file;
			vector<int> posWithSufficientReadDepth;
			single_trajs_file.open(slPath.c_str());
			string line;
			while(getline(single_trajs_file,line)) {

				//Read traj into string of "words"
				//Words is of the form: pos, al1, al2, #times #time1 #A #C #G #T #1tot #time2 #A #C #G #T #2tot
				//e.g. 860 G T 2 5 1 1 3464 923 4389 9 0 0 0 0 0
				stringstream ss(line);
				istream_iterator<std::string> begin(ss);
				istream_iterator<std::string> end;
				vector<string> words(begin, end);

				//Get the loci position
				int pos = stoi(words[0]);

				//Get before and after read depth
				int depthBefore = stoi(words[9]);
				int depthAfter = stoi(words[15]);

				int minReadDepth = path->getMinReadDepth();
				if(depthBefore >= minReadDepth && depthAfter >= minReadDepth) {

					posWithSufficientReadDepth.push_back(pos);
				}

			}


			if(path->getPhysicalPosMatter() == true) {

				removeLowReadDepth(posWithSufficientReadDepth);
			} else {

				cout << "Cannot remove loci with low total read depth if physical positions have not been loaded. Exiting.\n";
				exit(1);
			}

		}

	}


	//Import full haplotypes if requested   
	if(path->getImportHaps() == true) {

		importFullHaps = true;
		string importHapsFile = path->getImportHapsFile();
		ifstream importHapsStream;
		importHapsStream.open(importHapsFile.c_str());
		string fhString;

		importedFullHaps.clear();
		while(importHapsStream >> fhString) { //Read full haplotype into fh;

			Sequence fh = Sequence(fhString);
			importedFullHaps.push_back(fh);
		}
		importHapsStream.close();
	}

	//Import haplotypes from Mahan
	if(path->getImportHapsMahan() == true) {

		importHapsMahan = true;
		//1) Open file
		//2) Read line at a time
		//3) Append haps to importedFullHaps (same container as for normal import of full haps)
		//4) Append freqB to importedqB;
		//5) Append freqA to importedqA;
		
		ifstream importFile;
		importFile.open(path->getImportHapsMahanFile());
		string hapNum, hap, qB, qA;

		//Import data until none left
		while(importFile >> hapNum) {

			importFile >> hap;
			importFile >> qB;
			importFile >> qA;

			//Get data
			Sequence hapSeq = Sequence(hap);
			double qBdouble = atof(qB.c_str());
			double qAdouble = atof(qA.c_str());

			//Add to containers
			importedFullHaps.push_back(hapSeq);
			importedqB.push_back(qBdouble);
			importedqA.push_back(qAdouble);
		}

	}

	//Check if Times.in file present, and if so, read in deltaDays
	if(path->getTimesFile().compare("") != 0) {

		ifstream timesFile;
		timesFile.open(path->getTimesFile());
		string tBeforeString, tAfterString;
		timesFile >> tBeforeString;
		timesFile >> tAfterString;
		int tBefore = atoi(tBeforeString.c_str());
		int tAfter = atoi(tAfterString.c_str());
		deltaDays = tAfter - tBefore;
		cout << "DeltaDays: " << deltaDays << "\n";
		
	} else {
		
		deltaDays = 1; //Default
	}
}


//Simulate data (takes into account selection and epistasis)
void DataPH::simulateData(SimParam *sp) {

	cout << "In dataPH\n";
	spPH  = dynamic_cast<SimParamPH *>(sp);

	//Output file
	ofstream* outputFile = (spPH->getOutputFile());

	int Nt = spPH->getNt();
	int numGenes= spPH->getNumGenes();
	int dim =spPH->getDim();
	if(numGenes > 1) { cout << "WRONG USE OF DATA CLASS. USE DataPHMG instead. Aborting.."; return; }
	Sequence selSeq = *(spPH->getSel());
	bool selPresent = false;
	for(int i=0;i<selSeq.getLength();i++) {
		if(selSeq.getBase(i) != '-') { selPresent = true; break; }
	}


	//Set up RNG
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	unsigned long int seed = spPH->getSeed();
	gsl_rng_set(r,seed);


	//If simulating data based on real data, we need to make a copy of the original data.
	//This is necessary, as the data may be altered, e.g. sophs or even  monomorphic loci may
	//be removed. These alterations are not allowed, as we need the read depth for each soph 
	//to match that of the original dataset. We also cannot have sophs that cover fewer loci
	//than the real dataset. However, in order to realise whether we need to alter the data, 
	//we essentially first have to alter the data, if that makes sense. As such, we make a
	//copy of the initial state of the data. In the case where alterations have been made, 
	//we rerun the simulations based off the data saved in the copy.
	DataPH::phDataSingleSim dataCopy;
	if(spPH->getUseExternalReadStatistic() == true) {	
		dataCopy = data;
	}


	bool simulationError = true;
	int simulationErrorCounter = 0;
	while(simulationError == true) { //Rerun simulation of data until error free
		simulationError = false; //Assume no errors, then check during generation of simulated data

		vector<Sequence> fullHaps;
		//Construct full haplotypes differently depending on whether dealing with pure simulated data or pseudo real data
		if(spPH->getUseExternalReadStatistic()==false) {

			//Get a diploid sequence with random alleles, but with some alleles constrained
			DiploidSequence fullHapAlleles = generateRandomSemiConstrainedDiploidSequence(selSeq,r);


			//Obtain all possible full haplotypes with these allele combinations
			vector<Sequence> fullHapsAll = generateFullHapsFromAlleles(fullHapAlleles);


			//Pick dim full haplotypes for continued use
			int nPick = dim;


			bool allSitesAreVariants = false; //I.e. if have three haplotypes ACT, AGT, ACC then the first position is monomorphic which is not allowed, so check for this
			while(allSitesAreVariants == false) {
				cout << "Finding full haplotypes:\n";
				fullHaps.clear();
				vector<int> picks = pickRandomInts(r, (int) fullHapsAll.size(), nPick);
				for(int i=0; i<nPick;i++) {
					fullHaps.push_back(fullHapsAll[picks[i]]);
					fullHaps[i].print();
				}

				allSitesAreVariants = true; //Assume true, then check
				for(int i=0; i<fullHaps[0].getLength(); i++) { //Loop over loci
					bool currentSiteIsVariant = false;
					for(unsigned int j=1; j<fullHaps.size(); j++) {
						if(fullHaps[j].getBase(i) != fullHaps[0].getBase(i)) { //WLOG compare all haplotypes to the zeroth haplotype
							currentSiteIsVariant = true;
							break;
						}
					}
					if(currentSiteIsVariant == false) {
						allSitesAreVariants = false;
						break;
					}
				}

			}

			//If based on real sampling depth information, construct full haplotypes from real data
			//or get from imported full haplotypes.
		} else {

// Most recent concensus regarding full haplotypes for simulated real data: We have decided to usee all the exhaustively
// generated haplotypes in all cases. Old code retained in case it proves useful to revert later.
//
//			if(importFullHaps == true) { //Import of full haplotypes triggered
//	
//				fullHaps = importedFullHaps;
//
//			} else if(spPH->getFilterHaplotypesMethod()==3) { //Filter haplotypes using method 3
//		
//				cout << "Using filter method 3.\n";
//				AnalyserPH aPH;
//				aPH.loadData(this);
//				fullHaps = aPH.filterHaps3(r, spPH->getC(), spPH->getGeneIndex());
//
//			} else {

			if(importHapsMahan == true) { //Import of full haplotypes AND frequencies triggered
			
				fullHaps = importedFullHaps;	

			} else {

				//Standard approach: Generate all possible haplotypes that match the data
				AnalyserPH aPH;
				vector<Sequence> phs = data.getSequences(); //Get all (real data) partial haplotype sequences
				aPH.constructFullHaps(&phs,&fullHaps);
			}
			
		}


		//Assign (random) frequencies to the full haplotypes
		vector<double> pFHbefore;

		if(importHapsMahan == true) { //Use imported frequencies

			cout << "Setting pFHbefore = importedqB as import of Mahan haplotypes and frequencies were triggered.\n";
			pFHbefore = importedqB;
		
		} else { //Default approach
			
			//If there are less than 100 haplotypes we need to ensure that we don't get haplotypes with very small
			//frequencies. Why? Because with a small number of haplotypes, chances are that only a single or a few
			//haplotypes cover a specific allele. If this or these haplotypes have very low frequencies, they may
			//never be sampled and as such the loci in question becomes monomorphic during sampling. We want to
			//avoid this. For large number of haplotypes it is statistically unlikely that this will be a problem.
			if(fullHaps.size() <= 100) {
				bool freqTooSmall = true;
				double minFreq = 1/((double)fullHaps.size() + 12); //I.e. minFreq=0.05 when fullHaps.size() = 8
				cout << "minFreq: " << minFreq << "\n";
				if(minFreq < 0.03) { minFreq = 0.0001; }
				while(freqTooSmall==true) {
					pFHbefore.clear();
					for(unsigned int i=0;i<fullHaps.size();i++) {
						pFHbefore.push_back(gsl_rng_uniform_pos(r)); //Uniform on (0,1]
					}
					rescale(pFHbefore);

					freqTooSmall = false;
					for(unsigned int i=0; i<pFHbefore.size(); i++) {
						if(pFHbefore[i] < minFreq) { //Ensure no haplotypes are at very low frequencies
							freqTooSmall = true;
						}
					}
				}
			} else {

				pFHbefore.clear();
				for(unsigned int i=0;i<fullHaps.size();i++) {
					pFHbefore.push_back(gsl_rng_uniform_pos(r)); //Uniform on (0,1]
				}
				rescale(pFHbefore);
			}

		}
		

		//print full haps
		cout << "Full haplotypes prior to selection and transmission:\n";
		for(unsigned int i=0; i<fullHaps.size(); i++) {
			cout << i <<" Freq: " << pFHbefore[i] << ". Sequence: ";
			fullHaps[i].print();

		}

		//Save full haplotype variables - only for debugging use
		fullHapsTrue = fullHaps;
		qBtrue = pFHbefore;


		//Decompose full haplotypes into all possible partial haplotypes
		decomposeFullHapsFixedSize(fullHaps, data,fullHaps[0].getLength(),spPH->getUseExternalReadStatistic());

		//Draw partial data with a noise parameter of C
		drawObservationsBefore(spPH->getNb(), pFHbefore, spPH->getC(),r);


		vector<double> pFHafter;
		if(importHapsMahan == true) { //Use imported frequencies

			cout << "Setting pFHafter = importedqA as import of Mahan haplotypes and frequencies were triggered.\n";
			pFHafter = importedqA;

			//print full haps
			cout << "Full haplotypes after:\n";
			(* outputFile) << "Full haplotypes:\n";
			for(unsigned int i=0; i<fullHaps.size(); i++) {
				cout << i <<" Freq: " << pFHafter[i] << ". Sequence: ";
				fullHaps[i].print();

				//Write to file
				for(int j=0; j<fullHaps[i].getLength(); j++) {
					(*outputFile) << fullHaps[i].getBase(j);
				}
				(*outputFile) << "\t" << pFHbefore[i] << "\t" << pFHafter[i] << "\n";
			}

		
		} else { //Default approach
	

			/*
			 /
			 / Transmission and selection
			 /
			 */
			//Selective pressure, if any
			vector<double> pFHselection = pFHbefore;
			if(selPresent == true) {

				//Get fitnesses
				vector<double> hapFit = computeSimulationHapFitT(fullHaps,spPH);

				//Act with selection
				double tot=0;
				for(unsigned int i=0;i<fullHaps.size();i++) {

					pFHselection[i] = pFHselection[i]*exp(hapFit[i]);
					tot += pFHselection[i]; //This line is correct. Think about it: pFH[i] is assigned to new value above
				}

				for(unsigned int i=0;i<fullHaps.size();i++) {
					pFHselection[i]=pFHselection[i]/tot;
				}
			}


			//Perform transmission
			vector<int> transmittedCounts = multinomialSampling(Nt,pFHselection,r);
			vector<double> pFHtransmission(transmittedCounts.begin(), transmittedCounts.end()); //Parse to double
			rescaleZeroAccepted(pFHtransmission); //Turn into frequencies

			//print full haps
			cout << "Full haplotypes after selection and transmission:\n";
			for(unsigned int i=0; i<fullHaps.size(); i++) {
				cout << i <<" Freq: " << pFHtransmission[i] << ". Sequence: ";
				fullHaps[i].print();
			}


			//Perform within-host growth consisting of numGenerations of replication
			int numGenerations = spPH->getNumGenerations();
			deltaDays = spPH->getDeltaDays(); //Number of days between donor and recipient sampling
			int N = Nt; //Current population size
			int growthFactor = spPH->getGrowthFactor();
			pFHafter = pFHtransmission; //Current population frequency
			for(int i=1; i<=numGenerations; i++) { //Loop over generations

				//Drift step
				//Leads to growthFactor fold increase in population size
				//The multinomial sampling process below only works for N < 2147483647 (integer range).
				//As such, if N, prior to the x fold increase, is > 2147483647/x, then
				//set the new N to integer max value 2147483647.
				if(N>2147483647/growthFactor) { N=2147483647; }
				else { N=growthFactor*N; }
				vector<int> counts = multinomialSampling(N,pFHafter,r);
				vector<double> countsDouble(counts.begin(), counts.end()); //Parse to double
				rescaleZeroAccepted(countsDouble);
				pFHafter = countsDouble;

				//Within-host selection step
				if(spPH->getWithinHostSelectionPresent() == true) {

					//Get fitnesses
					vector<double> hapFitG = computeSimulationHapFitG(fullHaps,spPH);

					//Act with selection
					double tot=0;
					for(unsigned int i=0;i<fullHaps.size();i++) {

						pFHafter[i] = pFHafter[i]*exp(hapFitG[i]);
						tot += pFHafter[i]; //This line is correct. Think about it: pFH[i] is assigned to new value above
					}

					for(unsigned int i=0;i<fullHaps.size();i++) {
						pFHafter[i]=pFHafter[i]/tot;
					}
				}
			}

			//print full haps
			cout << "Full haplotypes after growth and within-host selection:\n";
			(* outputFile) << "Full haplotypes:\n";
			for(unsigned int i=0; i<fullHaps.size(); i++) {
				cout << i <<" Freq: " << pFHafter[i] << ". Sequence: ";
				fullHaps[i].print();

				//Write to file
				for(int j=0; j<fullHaps[i].getLength(); j++) {
					(*outputFile) << fullHaps[i].getBase(j);
				}
				(*outputFile) << "\t" << pFHbefore[i] << "\t" << pFHselection[i] << "\t" << pFHtransmission[i] << "\t" << pFHafter[i] << "\n";
			}

		}

		//Perform after observation
		//Population is now of order 10^8 (after 24 hours)
		drawObservationsAfter(spPH->getNa(), pFHafter, spPH->getC(),r);


		
		if(spPH->getFilterData().compare("Transmission")==0) {
	
			//Remove low count observations consistent with our transmission framework
			removeLowCountObservations();
			
		} else if(spPH->getFilterData().compare("Samfire")==0)  { //Samfire method

			//The samfire method removes only trajectories, i.e. partial haplotypes, not entire sets of partial haplotypes
			removeLowCountSamfire();
		
		} else { //No filtering whatsover
		
		}

		//Remove monomorphic sites and check that length of data still matches requirement from user
		if(spPH->getFilterData().compare("Transmission")==0 || spPH->getFilterData().compare("Samfire")==0) {
	
			if(spPH->getRemoveMonomorphicSim() == true) {
			//Fix any monomorphic SNPs arising from removing low count observations
			//Need to be updated, so if any monomorphic sites are found, it returns an error flag (and doesn't alter the data)
				removeMonomorphicSites();
			}


			//Check if outcome has required number of sites. Sometimes, removal of monomorphic sites leads to final outcome with fewer sites than intended.
			//Compare length of sCDS with length of data
			if(data.getNumberOfSites() < selSeq.getLength()) {

				simulationErrorCounter++;
				if(simulationErrorCounter > 10) {
					cout << "Outcome has fewer sites than required and it hasn't been fixed in 10 resimulations. Exiting.\n";
					exit(1);
				}



				cout << "Outcome has fewer sites than required. Rerunning simulations:\n";
				(* outputFile) << "Outcome has fewer sites than required. Rerunning simulations:\n";
				simulationError = true;

				if(spPH->getUseExternalReadStatistic()==false) {
					data.clear();
	
				}  else {

					//Revert to the copy of the original data to ensure that any alterations in "removeLowCountObservations"
					//and "removeMonomorphicSites" are undone. Note that we only care about alterations leading to a removal
					//of monomorhic sites. It doesn't matter if the final dataset has fewer sophs than the original - afterall
					//we did draw data for every soph using the original total read depth, i.e. it is just a coincidence that
					//some/all partial haplotypes in that soph end up being removed due to low counts.
					data = dataCopy;
				}
			}
		}
	}



	//Print to screen and file
	data.print();
	(*outputFile) << "\nPartial haplotypes:\n";
	data.writeToFile(*outputFile);


	if(spPH->getFormat().compare("Mahan")==0 || spPH->getFormat().compare("MahanAndHaps")==0) { //Export in Mahan format
	
		/*
		 * Print in Mahan format
		 */
		//Set up output folder
		ofstream outputFileMahan;
		string filePathFolder = spPH->getPathToFolder();	
		int geneIndex = spPH->getGeneIndex();
		stringstream ss;
		ss << filePathFolder << "SimulatedData_Mahan_Gene_" << geneIndex << ".dat";
		outputFileMahan.open((ss.str()).c_str()); 

		//Define physical pos as [1,2,3..,n] where n=num loci in full haps
		vector<int> physicalPos;
		//for(int i =1; i<=fullHaps[0].getLength(); i++) {
		for(int i =1; i<=data.getNumberOfSites(); i++) {
			physicalPos.push_back(i);
		}

		//Write to file
		data.writeToFileMahan(outputFileMahan, physicalPos);

	}


	if(spPH->getFormat().compare("Haps")==0 || spPH->getFormat().compare("MahanAndHaps")==0) { //Export in Hap_dataX.dat format
	
		/*
		 * Print in Haps format
		 */
		//Define output folder
		ofstream outputFileMahan;
		string filePathFolder = spPH->getPathToFolder();	
		int geneIndex = spPH->getGeneIndex();
		vector<string> genes {"HA", "MP", "NA", "NP", "NS", "PA", "PB1", "PB2"};
		string geneString = genes[geneIndex]; 
		stringstream ss;
		ss << filePathFolder << geneString << "/"; //Folder for gene
		string filePathGeneFolder = ss.str().c_str();

		//Check if folder exists. If not create it. Note, the code does not
		//remove existing data in folder. I have decided that starting to 
		//removing files through C++ is a potential hazard. As such, the user
		//must ensure that there is no existing data in the folder.
		if(fileExists(filePathGeneFolder)!= true) {
			cout << "Ouput folder for gene " << genes[geneIndex] << " created. Path: " << filePathGeneFolder << ".\n";
			mkdir(filePathGeneFolder.c_str(),0775);
		} else {

			cout << "Output folder for gene " << genes[geneIndex] << " already exists. Path: " << filePathGeneFolder << ". Ensure that there is no pre-existing data in the folder!!!\n";
		}		


		//Define physical pos as [1,2,3..,n] where n=num loci in full haps
		vector<int> physicalPos;
		//for(int i =1; i<=fullHaps[0].getLength(); i++) {
		for(int i =1; i<=data.getNumberOfSites(); i++) {
			physicalPos.push_back(i);
		}

		//Write to file
		data.writeToFileHaps(filePathGeneFolder, physicalPos);


		//Write physical pos to file
		ss.str("");
		ss << filePathGeneFolder << "Positions.dat";
		string filePath = ss.str();
		ofstream outputFile;
		outputFile.open(filePath.c_str());

		for(unsigned int i=0; i<physicalPos.size(); i++) {

			if(i!=0) { outputFile << "\n"; }
			outputFile << physicalPos[i];
		}
		outputFile.close();

	}



	//Compute degree of overlap (e.g. 5000 reads of length 1 + 4000*L2 + 3000*L3 + 2000*L2 + 1000*L1 / (5000+4000+3000+2000+1000) = 2.333)
	double dooBefore = data.computeDegreeOfOverlapBefore();
	double dooAfter = data.computeDegreeOfOverlapAfter();
	cout << "Degree of overlap: " << dooBefore << "\t" << dooAfter << "\n\n";
	(*outputFile) << "Degree of overlap: " << dooBefore << "\t" << dooAfter << "\n";


	//Clean up
	gsl_rng_free(r);

}

void DataPH::simulateDataRealistic(SimParam *sp, vector<int>& physicalPos, DiploidSequence& sCDS, int geneLength) {

	spPH  = dynamic_cast<SimParamPH *>(sp);

	//Output file
	ofstream* outputFile = (spPH->getOutputFile());


	//Length of full haplotypes
	int lengthFH = spPH->getNumLoci();
	if(lengthFH<2) {
		cout << "Error, number of loci must be above 2. Aborting simulation.\n";
		return;
	}

	int Nt = spPH->getNt();
	int numGenes= spPH->getNumGenes();
	if(numGenes > 1) { cout << "WRONG USE OF DATA CLASS. USE DataPHMG instead. Aborting.."; return; }
	Sequence selSeq = *(spPH->getSel());
	vector<double> selMagVec = *(spPH->getSelMagVec());
	bool selPresent = false;
	for(int i=0;i<selSeq.getLength();i++) {
		if(selSeq.getBase(i) != '-') { selPresent = true; break; }
	}


	//Set up RNG
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	unsigned long int seed = spPH->getSeed();
	gsl_rng_set(r,seed);

	bool simulationError = true;
	int simulationErrorCounter = 0;
	while(simulationError == true) { //Rerun simulation of data until error free
		simulationError = false; //Assume no errors, then check during generation of simulated data


		//Get a diploid sequence with random alleles, but with some alleles constrained
		//Note: selSeq is a sequence whereas sCDS is a diploidSequence. sCDS is needed because we want
		//replicates to have the same major AND minor alleles. For single replicate simulations, we only
		//care about the allele upon which selection acts. The other allele may be chosen at random.
		//Hence we can use a sequence in these (simpler) cases.
		DiploidSequence fullHapAlleles = generateRandomSemiConstrainedDiploidSequence(sCDS,r);

		//Obtain all possible full haplotypes with these allele combinations
		vector<Sequence> fullHapsAll = generateFullHapsFromAlleles(fullHapAlleles);


		//Pick dim full haplotypes for continued use
		int nPick = spPH->getDim(); //Dim == required number of full haplotypes
		vector<Sequence> fullHaps;
		bool allSitesAreVariants = false; //I.e. if have three haplotypes ACT, AGT, ACC then the first position is monomorphic which is not allowed, so check for this
		while(allSitesAreVariants == false) {
			cout << "Finding full haplotypes:\n";
			fullHaps.clear();
			vector<int> picks = pickRandomInts(r, (int) fullHapsAll.size(), nPick);
			for(int i=0; i<nPick;i++) {
				fullHaps.push_back(fullHapsAll[picks[i]]);
				fullHaps[i].print();
			}

			allSitesAreVariants = true; //Assume true, then check
			for(int i=0; i<fullHaps[0].getLength(); i++) { //Loop over loci
				bool currentSiteIsVariant = false;
				for(unsigned int j=1; j<fullHaps.size(); j++) {
					if(fullHaps[j].getBase(i) != fullHaps[0].getBase(i)) { //WLOG compare all haplotypes to the zeroth haplotype
						currentSiteIsVariant = true;
						break;
					}
				}
				if(currentSiteIsVariant == false) {
					allSitesAreVariants = false;
					break;
				}
			}
		}

		
		//Assign (random) frequencies to the full haplotypes
		vector<double> pFHbefore;
		bool freqTooSmall = true;
		double minFreq = 1/((double)fullHaps.size() + 12); //I.e. minFreq=0.05 when fullHaps.size() = 8
		if(minFreq < 0.03) { minFreq = 0.0001; }
		while(freqTooSmall==true) {
			pFHbefore.clear();
			for(unsigned int i=0;i<fullHaps.size();i++) {
				pFHbefore.push_back(gsl_rng_uniform_pos(r));
			}
			rescale(pFHbefore);

			freqTooSmall = false;
			for(unsigned int i=0; i<pFHbefore.size(); i++) {
				if(pFHbefore[i] < minFreq) { //Ensure no haplotypes are at very low frequencies
					freqTooSmall = true;
				}
			}
		}

		//print full haps
		cout << "Full haplotypes prior to selection and transmission:\n";
		for(unsigned int i=0; i<fullHaps.size(); i++) {
			cout << i <<" Freq: " << pFHbefore[i] << ". Sequence: ";
			fullHaps[i].print();
		}

		//Save full haplotype variables - only for debugging use
		fullHapsTrue = fullHaps;
		qBtrue = pFHbefore;


		//Decompose full haplotypes into all possible partial haplotypes of length L
		decomposeFullHapsFixedSize(fullHaps, data,fullHaps[0].getLength(), false);

		//Perform partial haplotype observations
		//Parameters based on H5N1 WithinhostF13
		//Might need to update these to reflect a more recent dataset
		int depth = 102825; //Set these from spPH later
		double meanReadLength = 119.68; //100
		double stDevReadLength = 136.88; //10
		double meanGapLength = 61.96;
		double stDevGapLength = 104.4783231;
		drawObservationsBeforeRealistic(depth, geneLength, physicalPos, meanReadLength, stDevReadLength, meanGapLength, stDevGapLength, fullHaps, pFHbefore, spPH->getC(),r);


		/*
		   /
		   / Transmission and selection
		   /
		   */

		//Selective pressure, if any
		vector<double> pFHselection = pFHbefore;
		if(selPresent == true) {

			//Get fitnesses
			vector<double> hapFit = computeSimulationHapFitT(fullHaps,spPH);

			//Act with selection
			double tot=0;
			for(unsigned int i=0;i<fullHaps.size();i++) {

				pFHselection[i] = pFHselection[i]*exp(hapFit[i]);
				tot += pFHselection[i]; //This line is correct. Think about it: pFH[i] is assigned to new value above
			}

			for(unsigned int i=0;i<fullHaps.size();i++) {
				pFHselection[i]=pFHselection[i]/tot;
			}
		}


		//Perform transmission
		vector<int> transmittedCounts = multinomialSampling(Nt,pFHselection,r);
		vector<double> pFHtransmission(transmittedCounts.begin(), transmittedCounts.end()); //Parse to double
		rescaleZeroAccepted(pFHtransmission); //Turn into frequencies

		//print full haps
		cout << "Full haplotypes after selection and transmission:\n";
		for(unsigned int i=0; i<fullHaps.size(); i++) {
			cout << i <<" Freq: " << pFHtransmission[i] << ". Sequence: ";
			fullHaps[i].print();
		}


		
		//Perform within host growth for numGenerations of replication
		int numGenerations = spPH->getNumGenerations();
		deltaDays = spPH->getDeltaDays(); //Number of days between donor and recipient sampling
		int N = Nt; //Current population size
		int growthFactor = spPH->getGrowthFactor();
		cout << "growthFactor: " << growthFactor << "\n";
		vector<double> pFHafter = pFHtransmission; //Current population frequency
		for(int i=1; i<=numGenerations; i++) { //Loop over generations

			//Drift step
			//Leads to growthFactor fold increase in population size
			//The multinomial sampling process below only works for N < 2147483647 (integer range).
			//As such, if N, prior to the x fold increase, is > 2147483647/x, then
			//set the new N to integer max value 2147483647.
			if(N>2147483647/growthFactor) { N=2147483647; }
			else { N=growthFactor*N; }
			vector<int> counts = multinomialSampling(N,pFHafter,r);
			vector<double> countsDouble(counts.begin(), counts.end()); //Parse to double
			rescaleZeroAccepted(countsDouble);
			pFHafter = countsDouble;

			cout << "After drift step " << i << "\n";
			for(unsigned int j=0; j<fullHaps.size(); j++) {
				cout << j <<" Freq: " << pFHafter[j] << ". Sequence: ";
				fullHaps[j].print();
			}

			//Selection step
			//Requires within host selection input
			//Ignore for now
			//Growth selection step
			if(spPH->getWithinHostSelectionPresent() == true) {

				//Get fitnesses
				vector<double> hapFitG = computeSimulationHapFitG(fullHaps,spPH);

				//Act with selection
				double tot=0;
				for(unsigned int i=0;i<fullHaps.size();i++) {

					pFHafter[i] = pFHafter[i]*exp(hapFitG[i]);
					tot += pFHafter[i]; //This line is correct. Think about it: pFH[i] is assigned to new value above
				}

				for(unsigned int i=0;i<fullHaps.size();i++) {
					pFHafter[i]=pFHafter[i]/tot;
				}
			}

			cout << "After within-host selection step " << i << "\n";
			for(unsigned int j=0; j<fullHaps.size(); j++) {
				cout << j <<" Freq: " << pFHafter[j] << ". Sequence: ";
				fullHaps[j].print();
			}

		}

		//print full haps
		cout << "Full haplotypes after growth and within-host selection:\n";
		(* outputFile) << "Full haplotypes:\n";
		for(unsigned int i=0; i<fullHaps.size(); i++) {
			cout << i <<" Freq: " << pFHafter[i] << ". Sequence: ";
			fullHaps[i].print();

			//Write to file
			for(int j=0; j<fullHaps[i].getLength(); j++) {
				(*outputFile) << fullHaps[i].getBase(j);
			}
			(*outputFile) << "\t" << pFHbefore[i] << "\t" << pFHselection[i] << "\t" << pFHtransmission[i] << "\t" << pFHafter[i] << "\n";
		}

		//Perform after observation
		drawObservationsAfterRealistic(depth, geneLength, physicalPos, meanReadLength, stDevReadLength, meanGapLength, stDevGapLength, fullHaps, pFHafter, spPH->getC(),r);

		cout << "Printing data prior to removal of low count observations and monomorphic sites:\n";
		data.print();

		if(spPH->getFilterData().compare("Transmission")==0) {
	
			//Remove low count observations consistent with our transmission framework
			removeLowCountObservations();
			
		} else if(spPH->getFilterData().compare("Samfire")==0)  { //Samfire method

			//The samfire method removes only trajectories, i.e. partial haplotypes, not entire sets of partial haplotypes
			removeLowCountSamfire();
		
		} else { //No filtering whatsover
		
		}

		//Remove monomorphic sites and check that length of data still matches requirement from user
		if(spPH->getFilterData().compare("Transmission")==0 || spPH->getFilterData().compare("Samfire")==0) {
	
			if(spPH->getRemoveMonomorphicSim() == true) {
			//Fix any monomorphic SNPs arising from removing low count observations
			//Need to be updated, so if any monomorphic sites are found, it returns an error flag (and doesn't alter the data)
				removeMonomorphicSites();
			}


			//Check if outcome has required number of sites. Sometimes, removal of monomorphic sites leads to final outcome with fewer sites than intended.
			//Compare length of sCDS with length of data
			if(data.getNumberOfSites() < sCDS.getLength()) {

				simulationErrorCounter++;
				if(simulationErrorCounter > 10) {
					cout << "Outcome has fewer sites than required and it hasn't been fixed in 10 resimulations. Exiting.\n";
					exit(1);
				}

				cout << "Outcome has fewer sites than required. Rerunning simulations:\n";
				(* outputFile) << "Outcome has fewer sites than required. Rerunning simulations:\n";
				simulationError = true;
				data.clear();
	
			}
		}


	}


	//Print to screen and file
	data.print();
	(*outputFile) << "\nPartial haplotypes:\n";
	data.writeToFile(*outputFile);

	if(spPH->getFormat().compare("Mahan")==0 || spPH->getFormat().compare("MahanAndHaps")==0) { //Export in Mahan format
	
		/*
		 * Print in Mahan format
		 */
		//Set up output folder
		ofstream outputFileMahan;
		string filePathFolder = spPH->getPathToFolder();	
		int geneIndex = spPH->getGeneIndex();
		stringstream ss;
		ss << filePathFolder << "SimulatedData_Mahan_Gene_" << geneIndex << ".dat";
		outputFileMahan.open((ss.str()).c_str()); 

		//Define physical pos as [1,2,3..,n] where n=num loci in full haps
		vector<int> physicalPos;
		//for(int i =1; i<=fullHaps[0].getLength(); i++) {
		for(int i =1; i<=data.getNumberOfSites(); i++) {
			physicalPos.push_back(i);
		}

		//Write to file
		data.writeToFileMahan(outputFileMahan, physicalPos);

	}


	if(spPH->getFormat().compare("Haps")==0 || spPH->getFormat().compare("MahanAndHaps")==0) { //Export in Hap_dataX.dat format
	
		/*
		 * Print in Haps format
		 */
		//Define output folder
		ofstream outputFileMahan;
		string filePathFolder = spPH->getPathToFolder();	
		int geneIndex = spPH->getGeneIndex();
		vector<string> genes {"HA", "MP", "NA", "NP", "NS", "PA", "PB1", "PB2"};
		string geneString = genes[geneIndex]; 
		stringstream ss;
		ss << filePathFolder << geneString << "/"; //Folder for gene
		string filePathGeneFolder = ss.str().c_str();

		//Check if folder exists. If not create it. Note, the code does not
		//remove existing data in folder. I have decided that starting to 
		//removing files through C++ is a potential hazard. As such, the user
		//must ensure that there is no existing data in the folder.
		if(fileExists(filePathGeneFolder)!= true) {
			cout << "Ouput folder for gene " << genes[geneIndex] << " created. Path: " << filePathGeneFolder << ".\n";
			mkdir(filePathGeneFolder.c_str(),0775);
		} else {

			cout << "Output folder for gene " << genes[geneIndex] << " already exists. Path: " << filePathGeneFolder << ". Ensure that there is no pre-existing data in the folder!!!\n";
		}		


		//Define physical pos as [1,2,3..,n] where n=num loci in full haps
		vector<int> physicalPos;
		//for(int i =1; i<=fullHaps[0].getLength(); i++) {
		for(int i =1; i<=data.getNumberOfSites(); i++) {
			physicalPos.push_back(i);
		}

		//Write to file
		data.writeToFileHaps(filePathGeneFolder, physicalPos);


		//Write physical pos to file
		ss.str("");
		ss << filePathGeneFolder << "Positions.dat";
		string filePath = ss.str();
		ofstream outputFile;
		outputFile.open(filePath.c_str());

		for(unsigned int i=0; i<physicalPos.size(); i++) {

			if(i!=0) { outputFile << "\n"; }
			outputFile << physicalPos[i];
		}
		outputFile.close();

	}


	//Compute degree of overlap (e.g. 5000 reads of length 1 + 4000*L2 + 3000*L3 + 2000*L2 + 1000*L1 / (5000+4000+3000+2000+1000) = 2.333)
	double dooBefore = data.computeDegreeOfOverlapBefore();
	double dooAfter = data.computeDegreeOfOverlapAfter();
	cout << "Degree of overlap: " << dooBefore << "\t" << dooAfter << "\n";
	(*outputFile) << "Degree of overlap: " << dooBefore << "\t" << dooAfter << "\n";

	//Clean up
	gsl_rng_free(r);


}


void DataPH::decomposeFullHaps(vector<Sequence> & fullHaps, DataPH::phDataSingleSim & phdss, vector<Sequence> &phSeqs) {

	//Remove any existing data
	phdss.clear();

	//Find the length of the longest haplotype
	int longestHap = 0;
	for(unsigned int i=0; i<phSeqs.size();i++) {

		int start = 0;
		for(int j=0; j<phSeqs[i].getLength(); j++) {
			if(phSeqs[i].getBase(j)!='-') { break; }
			start++;
		}
		int end = start+1;
		for(int j=end;j<phSeqs[i].getLength();j++) {
			if(phSeqs[i].getBase(j)=='-') { break; }
			end++;
		}

		int diff = end - start;
		if(diff > longestHap) { longestHap = diff; }
	}	

	decomposeFullHapsFixedSize(fullHaps,phdss,longestHap, false);

}

void DataPH::decomposeFullHapsFixedSize(vector<Sequence> &fullHaps, DataPH::phDataSingleSim &phdss, int size, bool useExternalReadStatistic) {

	for(int i=size; i>0;i--) {
		decomposeFullHapsLevel(fullHaps,phdss,i,useExternalReadStatistic);
	}
}

void DataPH::decomposeFullHapsLevel(vector<Sequence> & fullHaps, DataPH::phDataSingleSim & phdss, int level, bool useExternalReadStatistic) {

	int length = fullHaps[0].getLength();
	int numOfHaps = (int) fullHaps.size();

	int combinations = length + 1 -level; //e.g. for L=5, level=3, then combinations are 3: xxx--, -xxx-, --xxx.
	for(int i=0;i<combinations;i++) { //Loop over combinations, e.g. xxx--, -xxx-, --xxx
		setOfPartialHaps sph;

		for(int j=0;j<numOfHaps;j++) { //Loop over the full haplotypes to generate partial haplotypes from them

			Sequence currentSeq = fullHaps[j].subset(i,i+level); //subset: [i,i+level)

			//Add "-"s to Sequence if shorter than length
			if(i>0) {
				for(int k=0;k<i;k++) {
					currentSeq.addBaseFront('-');
				}
			}
			if((i+level) < length) {
				for(int k=i+level;k<length;k++) {
					currentSeq.addBase('-');
				}
			}

			bool partialHapUnique = true;

			for(int k=0;k<sph.getNumOfPartialHaps();k++) { //Loop over already found partial haps to check if new partial hap already exist

				Sequence sphSeq = sph.getSeq(k);

				if(identicalSeqs(currentSeq,sphSeq)==true) { //Partial hap already exists
					sph.addContrib(k, j); //add contribution to partial hap k
					partialHapUnique = false;
					break;
				}
			}

			//If partial hap has not been added before, add it to sph
			if(partialHapUnique==true) {
				partialHap phCurrent;
				phCurrent.setSeq(currentSeq);
				phCurrent.addContrib(j);
				sph.addPartialHap(phCurrent);
			}
		}

		//Find the loci covered by this soph
		vector<int> lociCovered;
		Sequence ph0Seq = sph.getSeq(0); //Consider zeroth haplotyp WLOG
		for(int l=0; l<ph0Seq.getLength(); l++) {
			if(ph0Seq.getBase(l) != '-') {

				lociCovered.push_back(l);
			}
		}
		sph.setLociCovered(lociCovered);


		//Case where external read information used, only include sophs that match up with read data
		if(useExternalReadStatistic == true) {

			//Loop over all sophs in phdss to see if it matches lociCovered
			//If it does, merge into phdss (i.e. new haplotypes, but keeping NbTot and NaTot)
			for(int j=0; j<phdss.getNumDatasets();  j++) { //Loop over sophs


				DataPH::setOfPartialHaps* sophCurrent = phdss.getSetOfPartialHaps(j);

				vector<int> lociCoveredReal = sophCurrent->getLociCovered();

				if(lociCoveredReal.size() == lociCovered.size()) { //Prerequisite to comparing is that they have same length
					if(identicalIntVectors(lociCoveredReal, lociCovered) == true) {

						cout << "Priting soph current seqs before swap-over:\n";
						for(int k=0; k<sophCurrent->getNumOfPartialHaps(); k++) { 

							sophCurrent->getPartialHap(k).getSeq().print();
						}

						//First remove real partial haps
						for(int k=sophCurrent->getNumOfPartialHaps()-1; k>=0; k--) { //Go through backwards

							sophCurrent->removePH(k);
						}

						//Then add new simulated partial haps
						for(int k=0; k<sph.getNumOfPartialHaps(); k++) {

							DataPH::partialHap phSim = sph.getPartialHap(k);
							sophCurrent->addPartialHap(phSim);
						}

						cout << "Priting soph current seqs after swap-over:\n";
						for(int k=0; k<sophCurrent->getNumOfPartialHaps(); k++) { 

							sophCurrent->getPartialHap(k).getSeq().print();
						}
						cout << "----------------\n";

					}
				} 
			}

		} else {

			phdss.addSetOfPartialHaps(sph); //Add set of partial haps to data
		}
	}

}

void DataPH::drawObservationsBefore(int Nb, vector<double> freqs, double C, gsl_rng *r) {

	int depth = Nb;

	for(int i=0; i<data.getNumDatasets(); i++) { //Loop over sophs

		setOfPartialHaps* sPH = data.getSetOfPartialHaps(i);
		if(spPH->getUseExternalReadStatistic() == true) {

			depth = sPH->getNbTot(); //Get depth from real data 
		}

		//Perform full haplotype draws
		vector<int> draws = DirMultSampling(depth, freqs, C, r);


		vector<int> phDraws(sPH->getNumOfPartialHaps(), 0); //Set all draws to zero

		//Get partial haplotype draws
		for(int j=0; j<sPH->getNumOfPartialHaps();j++) {
			vector<int> contribs = * (sPH->getContribs(j));
			for(unsigned int k=0; k<contribs.size();k++) {
				phDraws[j] += draws[contribs[k]];
			}
		}

		//Add draws to data
		sPH->setNbTot(depth);
		for(int j=0; j<sPH->getNumOfPartialHaps();j++) {
			sPH->setNb(j, phDraws[j]);
		}

	}
}

void DataPH::drawObservationsAfter(int Na, vector<double> freqs, double C, gsl_rng *r) {

	int depth = Na; //Depth we aim to have


	for(int i=0; i<data.getNumDatasets(); i++) { //Loop over sophs

		setOfPartialHaps* sPH = data.getSetOfPartialHaps(i);
		if(spPH->getUseExternalReadStatistic() == true) {

			depth = sPH->getNaTot(); //Get depth from real data 
		}

		//Perform full haplotype draws
		vector<int> draws = DirMultSampling(depth, freqs, C, r);


		vector<int> phDraws(sPH->getNumOfPartialHaps(), 0); //Set all draws to zero

		//Get partial haplotype draws
		for(int j=0; j<sPH->getNumOfPartialHaps();j++) {
			vector<int> contribs = * (sPH->getContribs(j));
			for(unsigned int k=0; k<contribs.size();k++) {
				phDraws[j] += draws[contribs[k]];
			}
		}

		//Add draws to data
		sPH->setNaTot(depth);
		for(int j=0; j<sPH->getNumOfPartialHaps();j++) {
			sPH->setNa(j, phDraws[j]);
		}

	}

}


void DataPH::drawObservationsBeforeRealisticWrong(int depth, int length, vector<int>& physicalPos, double meanReadLength, double stDevReadLength, double meanGapLength, double stDevGapLength, vector<Sequence>& fullHaps, vector<double>& freqs, double C, gsl_rng *r) {

	//Generate depth # of reads
	vector<DataPH::pairedEndRead> reads;
	for(int i=0; i<depth; i++) {

		DataPH::pairedEndRead read;
		//Generates a single paired-end read with a haplotype ID (i.e. integer) obtained from DirMult(n=1,C,freqs) sampling
		read.generateRead(meanReadLength, stDevReadLength, meanGapLength, stDevGapLength, length, r, C, freqs);
		reads.push_back(read);
		//	cout << "Printing read " << i << ":\n";
		//	read.print();
	}	

	//Initiate single locus trajectories. Initially set major = allele in fullHaps[0]
	slTrajs.clear();
	for(unsigned int i=0; i<physicalPos.size(); i++) {

		slTraj traj;
		traj.setPos(physicalPos[i]);
		traj.setMajor(fullHaps[0].getBase(i));
		//Find minor
		for(unsigned int j=0; j<fullHaps.size(); j++) {

			if(fullHaps[j].getBase(i) != fullHaps[0].getBase(i)) {
			
				traj.setMinor(fullHaps[j].getBase(i)); break;
			}
		}
		slTrajs.push_back(traj);
	}

	//Get single locus data
	int readsNotCoveringVariants = 0;
	for(unsigned int i=0; i<reads.size(); i++) {

		//Loop over all positions to check if they are in read
		bool found = false;
		for(unsigned int j=0; j<physicalPos.size(); j++) {
		
			if(reads[i].posCovered(physicalPos[j]) == true) { //Physical pos is covered by read

				char readChar = fullHaps[reads[i].getHapID()].getBase(j);
				if(slTrajs[j].getMinor() == readChar) {

					slTrajs[j].incrementNbMinor();

				} else if(slTrajs[j].getMajor() == readChar) {

					slTrajs[j].incrementNbMajor();

				} else {
				
					cout << "Error in generating single locus trajectories: Multiallelic site.\n"; exit(1);
				}

				found = true;
			}
		}
		if(found == false) {
			readsNotCoveringVariants++;
		}
	}
	cout << "Number of reads not covering variants: " << readsNotCoveringVariants << ". Total read depth: " << depth << " \n";
	//Correct major/minor alleles
	for(unsigned int i=0; i<slTrajs.size(); i++) {

		slTrajs[i].countNbTot(); //Counts NbTot and updates major/minor
		slTrajs[i].print();
	}

	//Loop over sets of partial haplotypes
	//NOTE: it is important that these are ordered such that the sophs covering most loci comes first
	//e.g. starting with xxxxx, then xxxx- and -xxxx, and then xxx--, -xxx-, --xxx, ..., ----x.
	for(int i=0; i<data.getNumDatasets(); i++) {

		//Get positions covered by this set of partial haps	        
		setOfPartialHaps* soph = data.getSetOfPartialHaps(i);
		vector<int> posCoveredBySOPH = soph->getPosCovered(physicalPos);

		vector<int> fullHaplotypeCounts(fullHaps.size(),0);
		int depthSOPH = 0;
		for(int j=(int)reads.size()-1; j>=0; j--) { //Go through reads from the back, so we can remove them when matched to soph

			//Checks if read covers the positions in the soph
			//As we are looping over sophs, starting with the sophs covering the most loci
			//we are guaranteed that if a read covers the soph, it cannot also cover a soph
			//with more loci. The read can, however, cover sophs with fewer loci, so we must
			//remove the read, once found
			if(reads[j].posCovered(posCoveredBySOPH) == true) {

				//Next we increase the haplotype count for the haplotype correspond to this read
				fullHaplotypeCounts[reads[j].getHapID()]++;
				reads.erase(reads.begin()+j); //Exclude this read from now on	
				depthSOPH++;
			}
		}

		vector<int> phDraws(soph->getNumOfPartialHaps(), 0); //Set all ph draws to zero

		if(depthSOPH > 0) { //Soph was sampled, so distribute unto partial haplotypes according to frequencies


			//Get partial haplotype draws
			for(int j=0; j<soph->getNumOfPartialHaps();j++) {
				vector<int> contribs = * (soph->getContribs(j));
				for(unsigned int k=0; k<contribs.size();k++) {
					phDraws[j] += fullHaplotypeCounts[contribs[k]];
				}
			}
		}

		//Add draws to data
		soph->setNbTot(depthSOPH);
		for(int j=0; j<soph->getNumOfPartialHaps();j++) {
			soph->setNb(j, phDraws[j]);
		}

	}

	cout << "Number of reads that didn't cover any partial haplotypes: " << reads.size() << ". Total read depth: " << depth << "\n";
}


void DataPH::drawObservationsBeforeRealistic(int depth, int length, vector<int>& physicalPos, double meanReadLength, double stDevReadLength, double meanGapLength, double stDevGapLength, vector<Sequence>& fullHaps, vector<double>& freqs, double C, gsl_rng *r) {

	//Generate depth # of reads
	vector<DataPH::pairedEndRead> reads;
	for(int i=0; i<depth; i++) {

		DataPH::pairedEndRead read;
		//Generates a single paired-end read with a haplotype ID (i.e. integer) obtained from DirMult(n=1,C,freqs) sampling
		read.generateRead(meanReadLength, stDevReadLength, meanGapLength, stDevGapLength, length, r, C, freqs);
		reads.push_back(read);
		//	cout << "Printing read " << i << ":\n";
		//	read.print();
	}

	//Initiate single locus trajectories. Initially set major = allele in fullHaps[0]
	slTrajs.clear();
	for(unsigned int i=0; i<physicalPos.size(); i++) {

		slTraj traj;
		traj.setPos(physicalPos[i]);
		traj.setMajor(fullHaps[0].getBase(i));
		//Find minor
		for(unsigned int j=0; j<fullHaps.size(); j++) {

			if(fullHaps[j].getBase(i) != fullHaps[0].getBase(i)) {
			
				traj.setMinor(fullHaps[j].getBase(i)); break;
			}
		}
		slTrajs.push_back(traj);
	}


	//Loop over sets of partial haplotypes
	//NOTE: it is important that these are ordered such that the sophs covering most loci comes first
	//e.g. starting with xxxxx, then xxxx- and -xxxx, and then xxx--, -xxx-, --xxx, ..., ----x.
	for(int i=0; i<data.getNumDatasets(); i++) {

		//Get positions covered by this set of partial haps	        
		setOfPartialHaps* soph = data.getSetOfPartialHaps(i);
		vector<int> posCoveredBySOPH = soph->getPosCovered(physicalPos);
		vector<int> lociCoveredBySOPH = soph->getLociCovered();
		cout << "Pos covered by soph " << i << ": "; printIntVector(posCoveredBySOPH);

		//Get number of reads covering this set of partial haps
		int depthSOPH = 0;
		for(int j=(int)reads.size()-1; j>=0; j--) { //Go through reads from the back, so we can remove them when matched to soph

			//Checks if read covers the positions in the soph
			//As we are looping over sophs, starting with the sophs covering the most loci
			//we are guaranteed that if a read covers the soph, it cannot also cover a soph
			//with more loci. The read can, however, cover sophs with fewer loci, so we must
			//remove exclude the read, once found
			if(reads[j].posCovered(posCoveredBySOPH) == true) {

				depthSOPH++;
				reads.erase(reads.begin()+j); //Exclude this read from now on	
			}
		}

		vector<int> phDraws(soph->getNumOfPartialHaps(), 0); //Set all ph draws to zero

		if(depthSOPH > 0) { //Soph was sampled, so distribute unto partial haplotypes according to frequencies

			//Perform full haplotype draws
			vector<int> draws = DirMultSampling(depthSOPH, freqs, C, r);


			//Get partial haplotype draws
			for(int j=0; j<soph->getNumOfPartialHaps();j++) {
				vector<int> contribs = * (soph->getContribs(j));
				for(unsigned int k=0; k<contribs.size();k++) {
					phDraws[j] += draws[contribs[k]];
				}
			}


			//Get SL data
			for(unsigned int j=0; j<lociCoveredBySOPH.size(); j++) {

				int loci = lociCoveredBySOPH[j];
				for(unsigned int k=0; k<draws.size(); k++) {

					int count = draws[k];
					char allele = fullHaps[k].getBase(loci);
					if(slTrajs[loci].getMajor() == allele) {
					
						slTrajs[loci].incrementNbMajor(count);

					} else if(slTrajs[loci].getMinor() == allele) {

						slTrajs[loci].incrementNbMinor(count);

					} else {
				
						cout << "Error in generating single locus trajectories: Multiallelic site.\n"; exit(1);
					}

				}
			}

		}

		//Add draws to data
		soph->setNbTot(depthSOPH);
		for(int j=0; j<soph->getNumOfPartialHaps();j++) {
			soph->setNb(j, phDraws[j]);
		}
	}

	//Correct major/minor alleles
	for(unsigned int i=0; i<slTrajs.size(); i++) {

		slTrajs[i].countNbTot(); //Counts NbTot and updates major/minor
		slTrajs[i].print();
	}
}


void DataPH::drawObservationsAfterRealisticWrong(int depth, int length, vector<int>& physicalPos, double meanReadLength, double stDevReadLength, double meanGapLength, double stDevGapLength, vector<Sequence>& fullHaps, vector<double>& freqs, double C, gsl_rng *r) {

	//Generate depth reads
	vector<DataPH::pairedEndRead> reads;
	for(int i=0; i<depth; i++) {

		DataPH::pairedEndRead read;
		//Generates a single paired-end read with a haplotype ID (i.e. integer) obtained from DirMult(n=1,C,freqs) sampling
		read.generateRead(meanReadLength, stDevReadLength, meanGapLength, stDevGapLength, length, r, C, freqs);
		reads.push_back(read);
	}	

	//Get single locus data
	int readsNotCoveringVariants = 0;
	for(unsigned int i=0; i<reads.size(); i++) {

		//Loop over all positions to check if they are in read
		bool found = false;
		for(unsigned int j=0; j<physicalPos.size(); j++) {
		
			if(reads[i].posCovered(physicalPos[j]) == true) { //Physical pos is covered by read

				char readChar = fullHaps[reads[i].getHapID()].getBase(j);
				if(slTrajs[j].getMinor() == readChar) {

					slTrajs[j].incrementNaMinor();
					found = true;

				} else if(slTrajs[j].getMajor() == readChar) {

					slTrajs[j].incrementNaMajor();
					found = true;

				} else {
				
					cout << "Error in generating single locus trajectories: Multiallelic site.\n"; exit(1);
				}

			}
		}
		if(found == false) {
			readsNotCoveringVariants++;
		}
	}
	cout << "Number of reads not covering variants: " << readsNotCoveringVariants << ". Total read depth: " << depth << " \n";
	//Count NaTot
	for(unsigned int i=0; i<slTrajs.size(); i++) {

		slTrajs[i].countNaTot(); 
		slTrajs[i].print();
	}

	//Loop over sets of partial haplotypes
	//NOTE: it is important that these are ordered such that the sophs covering most loci comes first
	//e.g. starting with xxxxx, then xxxx- and -xxxx, and then xxx--, -xxx-, --xxx, ..., ----x.
	for(int i=0; i<data.getNumDatasets(); i++) {

		//Get positions covered by this set of partial haps	        
		setOfPartialHaps* soph = data.getSetOfPartialHaps(i);
		vector<int> posCoveredBySOPH = soph->getPosCovered(physicalPos);
		cout << "Pos covered by soph " << i << ": "; printIntVector(posCoveredBySOPH);

		vector<int> fullHaplotypeCounts(fullHaps.size(),0);
		int depthSOPH = 0;
		for(int j=(int)reads.size()-1; j>=0; j--) { //Go through reads from the back, so we can remove them when matched to soph

			//Checks if read covers the positions in the soph
			//As we are looping over sophs, starting with the sophs covering the most loci
			//we are guaranteed that if a read covers the soph, it cannot also cover a soph
			//with more loci. The read can, however, cover sophs with fewer loci, so we must
			//remove the read, once found
			if(reads[j].posCovered(posCoveredBySOPH) == true) {

				//Next we increase the haplotype count for the haplotype correspond to this read
				fullHaplotypeCounts[reads[j].getHapID()]++;
				reads.erase(reads.begin()+j); //Exclude this read from now on	
				depthSOPH++;
			}
		}

		vector<int> phDraws(soph->getNumOfPartialHaps(), 0); //Set all ph draws to zero

		if(depthSOPH > 0) { //Soph was sampled, so distribute unto partial haplotypes according to frequencies


			//Get partial haplotype draws
			for(int j=0; j<soph->getNumOfPartialHaps();j++) {
				vector<int> contribs = * (soph->getContribs(j));
				for(unsigned int k=0; k<contribs.size();k++) {
					phDraws[j] += fullHaplotypeCounts[contribs[k]];
				}
			}
		}

		//Add draws to data
		soph->setNaTot(depthSOPH);
		for(int j=0; j<soph->getNumOfPartialHaps();j++) {
			soph->setNa(j, phDraws[j]);
		}

	}

	cout << "Number of reads that didn't cover any partial haplotypes: " << reads.size() << ". Total read depth: " << depth << "\n";
}


void DataPH::drawObservationsAfterRealistic(int depth, int length, vector<int>& physicalPos, double meanReadLength, double stDevReadLength, double meanGapLength, double stDevGapLength, vector<Sequence>& fullHaps, vector<double>& freqs, double C, gsl_rng *r) {

	//Generate depth # of reads
	vector<DataPH::pairedEndRead> reads;
	for(int i=0; i<depth; i++) {

		DataPH::pairedEndRead read;
		//Generates a single paired-end read with a haplotype ID (i.e. integer) obtained from DirMult(n=1,C,freqs) sampling
		read.generateRead(meanReadLength, stDevReadLength, meanGapLength, stDevGapLength, length, r, C, freqs);
		reads.push_back(read);
		//	cout << "Printing read " << i << ":\n";
		//	read.print();
	}


	//Loop over sets of partial haplotypes
	//NOTE: it is important that these are ordered such that the sophs covering most loci comes first
	//e.g. starting with xxxxx, then xxxx- and -xxxx, and then xxx--, -xxx-, --xxx, ..., ----x.
	for(int i=0; i<data.getNumDatasets(); i++) {

		//Get positions covered by this set of partial haps	        
		setOfPartialHaps* soph = data.getSetOfPartialHaps(i);
		vector<int> posCoveredBySOPH = soph->getPosCovered(physicalPos);
		vector<int> lociCoveredBySOPH = soph->getLociCovered();
		cout << "Pos covered by soph " << i << ": "; printIntVector(posCoveredBySOPH);

		//Get number of reads covering this set of partial haps
		int depthSOPH = 0;
		for(int j=(int)reads.size()-1; j>=0; j--) { //Go through reads from the back, so we can remove them when matched to soph

			//Checks if read covers the positions in the soph
			//As we are looping over sophs, starting with the sophs covering the most loci
			//we are guaranteed that if a read covers the soph, it cannot also cover a soph
			//with more loci. The read can, however, cover sophs with fewer loci, so we must
			//remove exclude the read, once found
			if(reads[j].posCovered(posCoveredBySOPH) == true) {

				depthSOPH++;
				reads.erase(reads.begin()+j); //Exclude this read from now on	
			}
		}

		vector<int> phDraws(soph->getNumOfPartialHaps(), 0); //Set all ph draws to zero

		if(depthSOPH > 0) { //Soph was sampled, so distribute unto partial haplotypes according to frequencies

			//Perform full haplotype draws
			vector<int> draws = DirMultSampling(depthSOPH, freqs, C, r);


			//Get partial haplotype draws
			for(int j=0; j<soph->getNumOfPartialHaps();j++) {
				vector<int> contribs = * (soph->getContribs(j));
				for(unsigned int k=0; k<contribs.size();k++) {
					phDraws[j] += draws[contribs[k]];
				}
			}


			//Get SL data
			for(unsigned int j=0; j<lociCoveredBySOPH.size(); j++) {

				int loci = lociCoveredBySOPH[j];
				for(unsigned int k=0; k<draws.size(); k++) {

					int count = draws[k];
					char allele = fullHaps[k].getBase(loci);
					if(slTrajs[loci].getMajor() == allele) {
					
						slTrajs[loci].incrementNaMajor(count);

					} else if(slTrajs[loci].getMinor() == allele) {

						slTrajs[loci].incrementNaMinor(count);

					} else {
				
						cout << "Error in generating single locus trajectories: Multiallelic site.\n"; exit(1);
					}

				}
			}

		}

		//Add draws to data
		soph->setNaTot(depthSOPH);
		for(int j=0; j<soph->getNumOfPartialHaps();j++) {
			soph->setNa(j, phDraws[j]);
		}
	}

	//Print SL data
	for(unsigned int i=0; i<slTrajs.size(); i++) {

		slTrajs[i].countNaTot(); //Counts NaTot
		slTrajs[i].print();
	}
}



//Samfire ensures that partial haplotypes with observations of less than 10 at at least one time point are removed
//To mimick samfire, we too remove such observations.
//This method is actually more strict than that. See comments in method below.
void DataPH::removeLowCountObservations() {

	//Loop over sets of partial haplotypes
	bool sophRemoved = false;
	int initialLength = data.getNumberOfSites(); //To be used later
	for(int i=(data.getNumDatasets()-1); i>=0; i--) { //Go through in reverse order, so we can remove sophs if necessary

		setOfPartialHaps* soph = data.getSetOfPartialHaps(i);

		//Remove any partial haps that do not have either 10 counts before or after transmission
		//This is in agreement with what samfire does
		for(int j=(soph->getNumOfPartialHaps()-1); j>=0;j--) { //Go through in reverse order, so we can remove partial haps if necessary


			if(soph->getNb(j) < 10 && soph->getNa(j) < 10) {

				soph->removePH(j);

			}

		}


		//Remove partial haplotypes with before counts below a threshold of 1%
		//Threshold needs to be recomputed every time a partial haplotype is removed
		//Partial haplotypes with lowest count must be removed first
		bool changeOccurred = true;
		while(changeOccurred == true) {

			//Recount Nb and Na
			soph->countNbTot();
			soph->countNaTot();

			if(soph->getNumOfPartialHaps() == 0) {	
				break; //Leave while loop
			}

			changeOccurred = false;
			int NbTot = soph->getNbTot();
			int cutOff = ceil(NbTot/100.0); //Cutoff frequency for observations is 1%

			//Find minimum Nb value
			vector<int> Nbs = soph->getNb();

			int minValue = Nbs[0];
			int minIndex = 0;
			for(unsigned int j=0; j<Nbs.size(); j++) {

				if(Nbs[j] < minValue) {

					minValue = Nbs[j];
					minIndex = j;
				}
			}

			if(minValue < cutOff || minValue == 0) { //Always remove zero before values

				soph->removePH(minIndex);
				changeOccurred = true;
			}

		}


		//Remove potential singletons/possible issues in after data
		//E.g. sometimes would find that an allele is present at substantial
		//frequency pre transmission but only found at e.g. a single copy after
		//transmission despite the sampling depth being large. This is likely
		//to be an error and can lead to some strong inferences of selection
		//on this locus.
		//
		//We fix this by indentifying these observations and setting them to
		//zero, i.e. we don't remove the haps entirely.
		//
		//Mathematically: If a partial haplotype observation x>0 after transmission is less than or equal to y and x is less than a fraction z of the total number of observations after transmission, then x is likely to be a sequencing error and is set to zero (x=0).
		int NaTot = soph->getNaTot();
		int y=1;
		double z=0.02;
		for(int j=0; j<soph->getNumOfPartialHaps(); j++) { //Loop over phs

			int Na = soph->getNa(j);
			if(Na > 0 && Na <= y) { //ph has low observations

				if(Na < z*NaTot) { //NaTot is sufficiently large that ph outcome of Na<=y is unlikely

					soph->setNa(j,0); //Set to zero to avoid issues
				}
			}	
		}




		//Extra care must be taken if soph contains just a single partial haplotype (either originally or after the above removal process).
		//In this case, a situation where either Nb or Na is zero leads to singularities.
		//Hence these partial haplotypes need to be removed as well.
		if(soph->getNumOfPartialHaps() == 1) {

			if(soph->getNb(0) == 0 || soph->getNa(0) == 0) {

				soph->removePH(0);
			}
		}

		//Remove set of partial haps if empty upon removing PHs
		//OR remove set of partial haps if NaTot is zero.
		//When computing variance in compound solution, the variance has a factor
		//of NaTot. So if NaTot is zero, the variance is zero and things diverge.
		//Note: NbTot is already non zero due to the above.
		if(soph->getNumOfPartialHaps() == 0 || soph->getNaTot() == 0) {
			sophRemoved = true;
			data.removeSOPH(i);
			cout << "Removing soph " << i << "\n";
		}
	}

	//Occasionally, when removing a soph, we might remove a soph which uniquely covers one or more loci, i.e. there are no other sophs
	//left in the dataset that covers these loci. As a result, we need to update lociCovered and physicalPos accordingly. Note that this
	//happens very rarely, but nonetheless from time to time.
	//Example: physPos = 100 200 300 400 500.
	//Loci = 0 1 2 3 4
	//If sophs uniquely covering 0 and 4 are removed, then we have to decrease lociCovered by 1 in all remaining sophs.
	//Additionally, physPos becomes 200 300 400.
	if(sophRemoved == true) {

		vector<int> lociFound; //Loci still existing

		//Find the loci left in the dataset
		for(int i=0; i<data.getNumDatasets(); i++) {

			setOfPartialHaps* soph = data.getSetOfPartialHaps(i);
			vector<int> lociCovered = soph->getLociCovered();
			for(unsigned int j=0; j<lociCovered.size(); j++) {

				bool alreadyExists = false;
				for(unsigned int k=0; k<lociFound.size(); k++) {
				
					if(lociCovered[j] == lociFound[k]) { alreadyExists = true; }
				}

				if(alreadyExists == false) { lociFound.push_back(lociCovered[j]); }
			}
		}

		//Find possible missing loci
		vector<int> lociMissing;
		for(int i=0; i<initialLength; i++) {

			bool lociExists = false;
			for(unsigned int j=0; j<lociFound.size(); j++) {

				if(i == lociFound[j]) { lociExists = true; }
			}
	
			if(lociExists == false) { lociMissing.push_back(i); }
		}

		//Going through missing loci in reverse and fix sophs and physPos
		for(int i=((int) lociMissing.size()) -1; i>= 0; i--) {

			//This updates the partial haplotype data and the loci covered
			//There is no need to merge and order haps as, we are not removing the locus due to monomorphism
			data.removeLocus(lociMissing[i]);	

			///Then correct physical pos if existing
			if(physicalPos.size() > 0) {
				physicalPos.erase(physicalPos.begin()+lociMissing[i]);
			}

			//Also correct imported haplotypes if existing
			if(importedFullHaps.size() > 0 ) {

				for(unsigned int j=0; j<importedFullHaps.size(); j++) {

					importedFullHaps[j].removeBase(lociMissing[i]);
					if(j==0) {
						cout << "Removing loci " << lociMissing[i] << " from full haplotypes.\n";
					}
				}
	
			}
		}
	}
}


//This method mimics the data filtering taking place in Samfire
//In particular, haplotypes with counts below 10 for both time points are removed.
//Additionally, haplotypes with frequencies below 1% for both time points are removed. The frequency cut-off is not recomputed as the haplotypes are removed.
void DataPH::removeLowCountSamfire() {

	//Loop over sets of partial haplotypes
	bool sophRemoved = false;
	int initialLength = data.getNumberOfSites(); //To be used later
	for(int i=(data.getNumDatasets()-1); i>=0; i--) { //Go through in reverse order, so we can remove sophs if necessary

		setOfPartialHaps* soph = data.getSetOfPartialHaps(i);

		//Remove any partial haps that do not have either 10 counts before or after transmission
		//This is in agreement with what samfire does
		for(int j=(soph->getNumOfPartialHaps()-1); j>=0;j--) { //Go through in reverse order, so we can remove partial haps if necessary


			if(soph->getNb(j) < 10 && soph->getNa(j) < 10) {

				soph->removePH(j);

			}

		}


		//Remove partial haplotypes with both before and after counts below a threshold of 1%

		//Count NbTot ad NaTot
		soph->countNbTot();
		soph->countNaTot();
		int NbTot = soph->getNbTot();
		int NaTot = soph->getNaTot();
		int cutOffBefore = ceil(NbTot/100.0); //Cutoff frequency for observations is 1%
		int cutOffAfter = ceil(NaTot/100.0); //Cutoff frequency for observations is 1%


		vector<int> Nbs = soph->getNb();
		vector<int> Nas = soph->getNa();

		//Loop over haps and remove if below cutoff for both time points
		for(int j=((int)Nbs.size()-1); j>=0; j--) {

			if(Nbs[j] <= cutOffBefore && Nas[j] <= cutOffAfter) {
				soph->removePH(j);
			}
		}
		//Recount total observations
		soph->countNbTot();
		soph->countNaTot();

		//Remove set of partial haps if empty upon removing PHs
		if(soph->getNumOfPartialHaps() == 0) {
			sophRemoved = true;
			data.removeSOPH(i);
			cout << "Removing soph " << i << "\n";
		}
	}

	//Occasionally, when removing a soph, we might remove a soph which uniquely covers one or more loci, i.e. there are no other sophs
	//left in the dataset that covers these loci. As a result, we need to update lociCovered and physicalPos accordingly. Note that this
	//happens very rarely, but nonetheless from time to time.
	//Example: physPos = 100 200 300 400 500.
	//Loci = 0 1 2 3 4
	//If sophs uniquely covering 0 and 4 are removed, then we have to decrease lociCovered by 1 in all remaining sophs.
	//Additionally, physPos becomes 200 300 400.
	if(sophRemoved == true) {

		vector<int> lociFound; //Loci still existing

		//Find the loci left in the dataset
		for(int i=0; i<data.getNumDatasets(); i++) {

			setOfPartialHaps* soph = data.getSetOfPartialHaps(i);
			vector<int> lociCovered = soph->getLociCovered();
			for(unsigned int j=0; j<lociCovered.size(); j++) {

				bool alreadyExists = false;
				for(unsigned int k=0; k<lociFound.size(); k++) {
				
					if(lociCovered[j] == lociFound[k]) { alreadyExists = true; }
				}

				if(alreadyExists == false) { lociFound.push_back(lociCovered[j]); }
			}
		}

		//Find possible missing loci
		vector<int> lociMissing;
		for(int i=0; i<initialLength; i++) {

			bool lociExists = false;
			for(unsigned int j=0; j<lociFound.size(); j++) {

				if(i == lociFound[j]) { lociExists = true; }
			}
	
			if(lociExists == false) { lociMissing.push_back(i); }
		}

		//Going through missing loci in reverse and fix sophs and physPos
		for(int i=((int) lociMissing.size()) -1; i>= 0; i--) {

			//This updates the partial haplotype data and the loci covered
			//There is no need to merge and order haps as, we are not removing the locus due to monomorphism
			data.removeLocus(lociMissing[i]);	

			///Then correct physical pos if existing
			if(physicalPos.size() > 0) {
				physicalPos.erase(physicalPos.begin()+lociMissing[i]);
			}

			//Also correct imported haplotypes if existing
			if(importedFullHaps.size() > 0 ) {

				for(unsigned int j=0; j<importedFullHaps.size(); j++) {

					importedFullHaps[j].removeBase(lociMissing[i]);
					if(j==0) {
						cout << "Removing loci " << lociMissing[i] << " from full haplotypes.\n";
					}
				}
	
			}
		}
	}
}

//Enforce a cutoff for observed after data such that observations with a frequency of
//less than cutOff is set to 0.
void DataPH::applyObsFrequencyCutOffAfter(double cutOff) {

	//Loop over sets of partial haplotypes
	for(int i=0; i<data.getNumDatasets(); i++) { 

		setOfPartialHaps* soph = data.getSetOfPartialHaps(i);

		//Find integer value of cutoff from total number of observations, NaTot
		int NaTot = soph->getNaTot();
		int cutOffValue = floor(cutOff*NaTot);

		//Set after observations to zero if partial hap has after observation of <= cutOffValue
		bool allBelowCutOff = true;
		for(int j=0; j<soph->getNumOfPartialHaps(); j++) {

			if(soph->getNa(j) <= cutOffValue) { 

				cout << "SOPH " << i << " PH " << j << " had after observation of " << soph->getNa(j) << " which is <= " << cutOff << " times the total number of observations for this soph (" << NaTot << "). The after observation of this PH has been set to 0.\n";
				soph->setNa(j,0); 
			} else {
				//	cout << "SOPH " << i << " PH " << j << " remain unchanged.\n";
				allBelowCutOff = false;
			}
		}

		if(allBelowCutOff == true) { //This happens if cutoff is too large

			cout << "All after observations are below threshold of " << cutOff << ". This outcome is undefined. Exiting.\n";
			exit(1);
		}

		soph->countNaTot(); //Recount NaTot after changes
	}
}

void DataPH::removeMonomorphicSites() {

	cout << "Removing any monomorphic sites:\n";

	bool monomorphicFound = true;
	while(monomorphicFound==true) { //Keep looking for and removing monomorphic sites until all have been dealt with

		monomorphicFound = false;
		int posToBeRemoved = -1;
		vector<Sequence> seqs = data.getSequences();
		int numLoci = seqs[0].getLength();
		cout << "Num loci: " << numLoci << "\n";
		cout << "Printing current data:\n";
		data.print();
		for(int i=0; i<numLoci;i++) {

			bool monomorphic = true;
			char nucleotide1 = '-';
			bool first = true;
			//Loop over ph sets and check loci one at a time
			for(unsigned int j=0; j<seqs.size();j++) {
				char current = seqs[j].getBase(i);
				if(current != '-') {
					if(first == true) {
						nucleotide1 = current;
						first = false;
					} else {
						if(current != nucleotide1) {
							monomorphic = false;
							break;
						}
					}
				}
			}

			//If monomorphic add loci to posToBeRemoved
			if(monomorphic == true) {
				posToBeRemoved = i;
				monomorphicFound = true;
				break;
			}
		}

		//Remove monomorphic site, redistributing the counts unto shorter sequences
		if(monomorphicFound == true) {
			cout << "Removing pos " << posToBeRemoved << "\n";
			data.removeLocus(posToBeRemoved);
			data.mergeHaps();
			data.orderHaps();	

			//Remove physical pos if existing
			if(physicalPos.size() > 0) {
				physicalPos.erase(physicalPos.begin()+posToBeRemoved);
			}
		}

		if(data.getNumDatasets() == 0) { 
			break; //No more sophs to remove
		}
	}

}


void DataPH::removeLowReadDepth(vector<int> posWithSufficientReadDepth) {

	for(int i=((int) physicalPos.size()-1); i>=0; i--) { //Go through in reverse so we can remove loci with complications


		bool lociFound = false;
		for(unsigned int j=0; j<posWithSufficientReadDepth.size(); j++) {

			if(physicalPos[i] == posWithSufficientReadDepth[j]) {
				lociFound = true; break;
			}
		}

		if(lociFound == false) { //Low read depth for this loci, so remove

			cout << "Removing pos " << physicalPos[i] << " (loci " << i << " ) due to low read depth.\n";
			data.removeLocus(i);
			data.mergeHaps();
			data.orderHaps();	

			//Remove physical pos
			physicalPos.erase(physicalPos.begin()+i);

		}
		
	}

}



//Checks if any of the partial haplotype sets cover the same positions
//If they do, they get merged
void DataPH::phDataSingleSim::mergeHaps() {

	//Merge haps that cover same positions
	vector<int> sophsToBeRemoved;
	cout << "Dataset.size() in mergeHaps: " << dataset.size() << "\n";
	cout << "Printing dataset:\n";
	for(unsigned int i=0; i<dataset.size(); i++) {
		cout << "i: " << i << "\n";
		dataset[i].print();
	}
	for(unsigned int i=0; i<dataset.size(); i++) { //Loop over sets of partial haps and try and merge i with j

		setOfPartialHaps* seti = &dataset[i];
		Sequence seqi0 = seti->getSeq(0); //Get sequence of zeroth ph
		//Find the position covered by this seq
		vector<int> posCoveredi;
		for(int l=0; l<seqi0.getLength(); l++) {

			if(seqi0.getBase(l) != '-') {

				posCoveredi.push_back(l);
			}
		}

		if(posCoveredi.size() == 0) { //The current soph covers nothing, so remove!

			cout << "SOPH " << i << " to be removed because it covers no loci. Printing soph:\n";
			seti->print();
			sophsToBeRemoved.push_back(i);

		} else { //Else compare with the other sophs to see if they cover same set of loci and can be merged 

			for(unsigned int j=i+1; j<dataset.size(); j++) {

				setOfPartialHaps* setj = &dataset[j];
				Sequence seqj0 = setj->getSeq(0); //Get sequence of zeroth ph


				//Find the position covered by this seq
				vector<int> posCoveredj;
				for(int l=0; l<seqj0.getLength(); l++) {

					if(seqj0.getBase(l) != '-') {

						posCoveredj.push_back(l);
					}
				}


				if(posCoveredi.size() == posCoveredj.size()) { //Same size, so can compare content of vectors

					if(identicalIntVectors(posCoveredi, posCoveredj) == true) { //Same set of loci covered, so merge

						cout << "posCovered for soph " << i << " and " << j << " identical.\n";
						
						//Merge i into j
						for(int k=0; k<seti->getNumOfPartialHaps(); k++) { //Loop over phs in soph i

							Sequence seqik = seti->getSeq(k);
							int Nbik = seti->getNb(k);
							int Naik = seti->getNa(k);

							bool merged = false;
							for(int l=0; l<setj->getNumOfPartialHaps(); l++) { //Compare phs k from soph i with phs l from soph j

								Sequence seqjl = setj->getSeq(l);
								int Nbjl = setj->getNb(l);
								int Najl = setj->getNa(l);

								if(identicalSeqs(seqik, seqjl)) { //Identical partial haps, so merge into j

									Nbjl += Nbik;
									Najl += Naik;
									setj->setNb(l,Nbjl);
									setj->setNa(l,Najl);
									merged = true;
									cout << "SOPH " << i << " ph " << k << " and  SOPH " << j << " ph " << l << " have been merged!\n";
									//	cout << "New Nb: " << setj->getNb(l) << "\n";
									//	cout << "New Na: " << setj->getNa(l) << "\n";
									break;
								}
							}

							if(merged == false) { //New partial haplotype found

								cout << "New partial haplotype add to soph " << j << ".\n";
								partialHap newPH;
								newPH.setSeq(seqik);
								newPH.setNb(Nbik);
								newPH.setNa(Naik);
								setj->addPartialHap(newPH);
							}

							//Recount NbTot and NaTot for soph j
							setj->countNbTot();
							setj->countNaTot();
						}

						sophsToBeRemoved.push_back(i);
						cout << "SOPH " << i << " to be removed. Printing soph:\n";
						seti->print();
						break; //soph i taken care of, move on to i+1
					}
				}

			}
		}

	}


	//Remove merged sophs
	for(int i=((int)sophsToBeRemoved.size())-1; i>=0; i--) { //Go through backwards

		cout << "Removing soph: " << sophsToBeRemoved[i] << "\n";
		removeSOPH(sophsToBeRemoved[i]);		
	}

}

//Order sets of haplotypes so that the longest come first
//Within sophs of similar lengths, the ones with the lowest position come first
void DataPH::phDataSingleSim::orderHaps() {

	bool changeOccurred = true;
	while(changeOccurred == true) {

		changeOccurred = false;
		for(int i=0; i<((int) dataset.size())-1; i++) { //Loop over sets of partial haps and try and swap i with j=i+1

			setOfPartialHaps* seti = &dataset[i];
			Sequence seqi0 = seti->getSeq(0); //Get sequence of zeroth ph
//			cout << "Sequence " << i << " zero: "; seqi0.print();
			//Find the position covered by this seq
			vector<int> posCoveredi;
			for(int l=0; l<seqi0.getLength(); l++) {

				if(seqi0.getBase(l) != '-') {

					posCoveredi.push_back(l);
				}
			}

			setOfPartialHaps* setj = &dataset[i+1];
			Sequence seqj0 = setj->getSeq(0); //Get sequence of zeroth ph
//			cout << "Sequence " << i+1 << " zero: "; seqj0.print();
			vector<int> posCoveredj;
			for(int l=0; l<seqj0.getLength(); l++) {

				if(seqj0.getBase(l) != '-') {

					posCoveredj.push_back(l);
				}
			}

			//Swap sophs if |j|>|i|, e.g. if j covers loci 2,3,4,5 and i covers 1,2,3, then swap
			if(posCoveredj.size() > posCoveredi.size()) {

				changeOccurred = true;

				//Swap sophs if |j|=|i| but j[0] < i[0], e.g. if j covers loci 1,2,3,4 and i covers loci 3,4,5,6 then swap
			} else if(posCoveredj.size() == posCoveredi.size() && posCoveredj[0] < posCoveredi[0]) {

				changeOccurred = true;
			}


			if(changeOccurred == true) { //Swap i and j


				cout << "Swapping sophs " << i << " and " << i+1 << ".\n"; 
				iter_swap(dataset.begin()+i, dataset.begin()+i+1);
				break; //Start over from i=0 again
			}
		}
	}

//	cout << "Done ordering.\n";

}


vector<Sequence> DataPH::getSequences() { return data.getSequences(); }


void DataPH::computeContribs(vector<Sequence> &fhs) {

	data.clearContributions();

	for(int i=0; i<data.getNumDatasets(); i++) {

		setOfPartialHaps* soph = data.getSetOfPartialHaps(i);
		// soph->clearContributions();
		vector<Sequence> phs = soph->getSequences();

		for(unsigned int j=0; j<phs.size(); j++) {

			Sequence partialHap = phs[j];
			vector<int> contribs;

			for(unsigned int k=0; k<fhs.size();k++) {

				//Get non '-' scope of partialHap
				int start = -1; //Initialise to minus one to catch errors
				int end = -1;
				bool first = true;
				for(int l=0;l<partialHap.getLength();l++) {

					if(partialHap.getBase(l) != '-' && first==true) {
						start = l;
						end = l+1;
						first = false;
					} else if(partialHap.getBase(l) != '-') {
						end = l+1;
					}
				}

				Sequence partialHapSubset = partialHap.subset(start, end);
				Sequence fullHapSubset = fhs[k].subset(start, end);

				//Check if partialHap is a subset of full hap
				if(identicalSeqs(partialHapSubset, fullHapSubset) == true) { contribs.push_back(k); }
			}

			for(unsigned int k=0; k<contribs.size(); k++) {
				soph->addContrib(j, contribs[k]);
			}
		}
	}
}


void DataPH::phDataSingleSim::computeContribs(vector<Sequence> &fhs) {

	clearContributions();

	for(unsigned int i=0; i<dataset.size(); i++) {

		vector<Sequence> phs = dataset[i].getSequences();

		for(unsigned int j=0; j<phs.size(); j++) {

			Sequence partialHap = phs[j];
			vector<int> contribs;

			for(unsigned int k=0; k<fhs.size();k++) {

				//Get non '-' scope of partialHap
				int start = -1; //Initialise to minus one to catch errors
				int end = -1;
				bool first = true;
				for(int l=0;l<partialHap.getLength();l++) {

					if(partialHap.getBase(l) != '-' && first==true) {
						start = l;
						end = l+1;
						first = false;
					} else if(partialHap.getBase(l) != '-') {
						end = l+1;
					}
				}

				Sequence partialHapSubset = partialHap.subset(start, end);
				Sequence fullHapSubset = fhs[k].subset(start, end);

				//Check if partialHap is a subset of full hap
				if(identicalSeqs(partialHapSubset, fullHapSubset) == true) { contribs.push_back(k); }
			}

			for(unsigned int k=0; k<contribs.size(); k++) {
				dataset[i].addContrib(j, contribs[k]);
			}
		}	
	}
}

vector<int> DataPH::getBeforeCounts(int sphIndex) { return data.getBeforeCounts(sphIndex); }
vector<int> DataPH::getAfterCounts(int sphIndex) { return data.getAfterCounts(sphIndex); }
int DataPH::numOfPHSets() { return data.numOfPHSets(); }




//Methods related to structs
DataPH::slTraj::slTraj() : NbMinor(0), NbMajor(0), NaMinor(0), NaMajor(0), NbTot(0), NaTot(0), pos(-1) { //Constructor initialises parameters to zero

	minor = '-';
	major = '-';
}
void DataPH::slTraj::setMinor(char c) { minor = c; }
void DataPH::slTraj::setMajor(char c) { major = c; }
void DataPH::slTraj::incrementNbMinor(int x /* = 1 */) { NbMinor += x; } //Default: x=1;
void DataPH::slTraj::incrementNbMajor(int x /* = 1 */) { NbMajor += x; }
void DataPH::slTraj::incrementNaMinor(int x /* = 1 */) { NaMinor += x; }
void DataPH::slTraj::incrementNaMajor(int x /* = 1 */) { NaMajor += x; }
int DataPH::slTraj::getNbMinor() { return NbMinor; }
int DataPH::slTraj::getNbMajor() { return NbMajor; }
int  DataPH::slTraj::getNaMinor() { return NaMinor; }
int DataPH::slTraj::getNaMajor() { return NaMajor; }

//Counts total counts before transmission
//Also updates major and minor alleles based on before time point
void DataPH::slTraj::countNbTot() {

	NbTot = NbMinor + NbMajor;
	if(NbMinor > NbMajor) { //Swap major/minor
		
		char trueMinor = major;
		major = minor;
		minor = trueMinor;

		int trueNbMinor = NbMajor;
		NbMajor = NbMinor;
		NbMinor = trueNbMinor;
	}
}
void DataPH::slTraj::countNaTot() { NaTot = NaMinor + NaMajor; }
int DataPH::slTraj::getNbTot() { return NbTot; }
int DataPH::slTraj::getNaTot() { return NaTot; }
void DataPH::slTraj::setPos(int p) { pos = p; }
int DataPH::slTraj::getPos() { return pos; }
char DataPH::slTraj::getMinor() { return minor; }
char DataPH::slTraj::getMajor() { return major; }
void DataPH::slTraj::print() {

	cout << pos << " " << minor << " " << major << " " << NbMinor << " " << NbMajor << " " << NbTot << " " << NaMinor << " " << NaMajor << " " << NaTot << "\n";

}





void DataPH::partialHap::setSeq(Sequence &s) { seq = s; }
void DataPH::partialHap::setNb(int NB) { Nb = NB; }
void DataPH::partialHap::setNa(int NA) { Na = NA; }
void DataPH::partialHap::addContrib(int c) { contribs.push_back(c); }
Sequence DataPH::partialHap::getSeq() { return seq; }
int DataPH::partialHap::getNb() { return Nb; }
int DataPH::partialHap::getNa() { return Na; }
vector<int>* DataPH::partialHap::getContribs() { return &contribs; }
void DataPH::partialHap::clearContributions() { contribs.clear(); }
void DataPH::partialHap::print() {

	cout << "Sequence: ";
	seq.print();

	cout << "Observations before: " << Nb << ". Observations after: " << Na << "\n";

	//cout << "Contributions ";
	//printIntVector(contribs);

}

void DataPH::setOfPartialHaps::setNb(int index, int Nb) { phVec[index].setNb(Nb); }
void DataPH::setOfPartialHaps::setNa(int index, int Na) { phVec[index].setNa(Na); }
void DataPH::setOfPartialHaps::setNbTot(int N) { NbTot = N; }
void DataPH::setOfPartialHaps::setNaTot(int N) { NaTot = N; }
void DataPH::setOfPartialHaps::countNbTot() {
	NbTot = 0;
	for(unsigned int i=0; i<phVec.size(); i++) {
		NbTot += phVec[i].getNb();
	}
}
void DataPH::setOfPartialHaps::countNaTot() {
	NaTot = 0;
	for(unsigned int i=0; i<phVec.size(); i++) {
		NaTot += phVec[i].getNa();
	}
}

int DataPH::setOfPartialHaps::getNbTot() {

	return NbTot;
}

int DataPH::setOfPartialHaps::getNaTot() {

	return NaTot;
}
void DataPH::setOfPartialHaps::addPartialHap(partialHap & ph) { phVec.push_back(ph); }
DataPH::partialHap DataPH::setOfPartialHaps::getPartialHap(int phIndex) { return phVec[phIndex]; }
Sequence DataPH::setOfPartialHaps::getSeq(int phIndex) { return phVec[phIndex].getSeq(); }
int DataPH::setOfPartialHaps::getNumOfPartialHaps() { return (int) phVec.size(); }
void DataPH::setOfPartialHaps::addContrib(int index, int contrib) { phVec[index].addContrib(contrib); }
vector<int>* DataPH::setOfPartialHaps::getContribs(int index) { return phVec[index].getContribs(); }
void DataPH::setOfPartialHaps::clearContributions() {
	for(unsigned int i=0;i<phVec.size();i++) {
		phVec[i].clearContributions();
	}
}
void DataPH::setOfPartialHaps::print() {

	for(unsigned int i=0; i<phVec.size(); i++) {
		phVec[i].print();
	}

	cout << "Total observations before: " << NbTot << ". Total observations after: " << NaTot << "\n";
}
vector<Sequence> DataPH::setOfPartialHaps::getSequences() {

	vector<Sequence> allSeqs;
	for(unsigned int i=0;i<phVec.size();i++) {
		allSeqs.push_back(phVec[i].getSeq());
	}

	return allSeqs;
}

int DataPH::phDataSingleSim::numOfPHSets() { return (int) dataset.size(); }



void DataPH::phDataSingleSim::addSetOfPartialHaps(setOfPartialHaps & setPH) {
	dataset.push_back(setPH);
}

void DataPH::phDataSingleSim::clear() { dataset.clear(); }

int DataPH::phDataSingleSim::getNumDatasets() { return (int) dataset.size(); }
DataPH::setOfPartialHaps* DataPH::phDataSingleSim::getSetOfPartialHaps(int index) { return & (dataset[index]); }
void DataPH::phDataSingleSim::print() {
	cout << "Printing " << dataset.size() << "sophs:\n";
	for(unsigned int i=0;i<dataset.size();i++) {
		dataset[i].print();
	}
}
void DataPH::phDataSingleSim::setNb(int index1, int index2, int Nb) { dataset[index1].setNb(index2,Nb); }

vector<Sequence> DataPH::phDataSingleSim::getSequences() {

	vector<Sequence> allSeqs;
	for(unsigned int i=0;i<dataset.size();i++) {
		vector<Sequence> phsSeqs = dataset[i].getSequences();
		for(unsigned int j=0; j<phsSeqs.size();j++) {
			allSeqs.push_back(phsSeqs[j]);
		}
	}

	return allSeqs;
}


vector<int> DataPH::phDataSingleSim::getBeforeCounts(int sphIndex) {

	vector<int> counts;

	vector<int> Nb = dataset[sphIndex].getNb(); //Vector of partial haplotype observations

	for(unsigned int i=0; i<Nb.size();i++) {

		vector<int>* contribs = dataset[sphIndex].getContribs(i);

		for(unsigned int k=0; k<contribs->size(); k++) {
			int c = (*contribs)[k];
			int currentLength = (int) counts.size();

			if((c+1) > currentLength) { //Dynamically update the size of counts as more contributions are discovered

				for(int l=currentLength; l<(c+1); l++) {
					counts.push_back(0);
				}
			}

			counts[c] += Nb[i];
		}
	}

	return counts;
}


vector<int> DataPH::phDataSingleSim::getAfterCounts(int sphIndex) {

	vector<int> counts;

	vector<int> Na = dataset[sphIndex].getNa(); //Vector of partial haplotype observations

	for(unsigned int i=0; i<Na.size();i++) {

		vector<int>* contribs = dataset[sphIndex].getContribs(i);

		for(unsigned int k=0; k<contribs->size(); k++) {
			int c = (*contribs)[k];
			int currentLength = (int) counts.size();

			if((c+1) > currentLength) { //Dynamically update the size of counts as more contributions are discovered

				for(int l=currentLength; l<(c+1); l++) {
					counts.push_back(0);
				}
			}

			counts[c] += Na[i];
		}
	}

	return counts;
}

vector<int> DataPH::setOfPartialHaps::getNb() {

	vector<int> Nb;
	for(unsigned int i=0; i<phVec.size(); i++) {
		Nb.push_back(phVec[i].getNb());
	}
	return Nb;
}

vector<int> DataPH::setOfPartialHaps::getNa() {

	vector<int> Na;
	for(unsigned int i=0; i<phVec.size(); i++) {
		Na.push_back(phVec[i].getNa());
	}
	return Na;
}

void DataPH::print() { data.print(); }

void DataPH::phDataSingleSim::clearContributions() {
	for(unsigned int i=0;i <dataset.size();i++) { dataset[i].clearContributions(); }
}

void DataPH::clearContributions() { data.clearContributions(); }
void DataPH::clear() { data.clear(); }
int DataPH::getNumOfDatasets() { return data.getNumDatasets(); }
int DataPH::phDataSingleSim::getNumOfPartialHaps(int datasetIndex) { return dataset[datasetIndex].getNumOfPartialHaps(); }
int DataPH::getNumOfPartialHaps(int datasetIndex) { return data.getNumOfPartialHaps(datasetIndex); }
vector<int>* DataPH::phDataSingleSim::getContribs(int sphIndex, int phIndex) { return dataset[sphIndex].getContribs(phIndex); }
vector<int> * DataPH::getContribs(int sphIndex, int phIndex) {
	return data.getContribs(sphIndex, phIndex);
}

vector<int>  DataPH::getNb(int sphIndex) { return data.getNb(sphIndex); }
vector<int>  DataPH::getNa(int sphIndex) { return data.getNa(sphIndex); }
vector<int>  DataPH::phDataSingleSim::getNb(int sphIndex) { return dataset[sphIndex].getNb(); }
vector<int>  DataPH::phDataSingleSim::getNa(int sphIndex) { return dataset[sphIndex].getNa(); }

vector<double> DataPH::getPHFreqs(int sphIndex, vector<double> & fhFreqs) {
	return data.getPHFreqs(sphIndex, fhFreqs);
}

vector<double> DataPH::phDataSingleSim::getPHFreqs(int sphIndex, vector<double> & fhFreqs) {
	return dataset[sphIndex].getPHFreqs(fhFreqs);
}
vector<double> DataPH::setOfPartialHaps::getPHFreqs(vector<double> & fhFreqs) {

	vector<double> p;
	for(unsigned int i=0; i<phVec.size(); i++) {

		vector<int> contribs = *phVec[i].getContribs();
		double pcurrent = 0;
		for(unsigned int j=0; j<contribs.size(); j++) {
			pcurrent += fhFreqs[contribs[j]];
		}
		p.push_back(pcurrent);
	}
	rescale(p); //Probably unnecessary. Check?

	return p;
}

vector<double> DataPH::getPHhapFit(int sphIndex, vector<double> &fhFreqs, vector<double> & fhHapFit) {
	return data.getPHhapFit(sphIndex, fhFreqs, fhHapFit);
}
vector<double> DataPH::phDataSingleSim::getPHhapFit(int sphIndex, vector<double> &fhFreqs, vector<double> & fhHapFit) {
	return dataset[sphIndex].getPHhapFit(fhFreqs, fhHapFit);
}
vector<double> DataPH::setOfPartialHaps::getPHhapFit(vector<double> &fhFreqs, vector<double> & fhHapFit) {


	vector<double> phHapFit;
	for(unsigned int i=0; i<phVec.size();i++) { //Loop over haplotypes in set

		vector<int> contribs = *phVec[i].getContribs();
		double hapFitCurrent = 0;
		double sumFHfreqs = 0;
		for(unsigned int j=0; j<contribs.size();j++) {
			hapFitCurrent += fhFreqs[contribs[j]]*exp(fhHapFit[contribs[j]]);
			sumFHfreqs += fhFreqs[contribs[j]];
		}
		phHapFit.push_back(log(hapFitCurrent/sumFHfreqs));
	}

	return phHapFit;
}

//Removes the partial hap at the phIndex position
//E.g. if vec={a,b,c,d,e,f}, then removePH(3) removes d
void DataPH::setOfPartialHaps::removePH(int phIndex) {
	phVec.erase(phVec.begin()+phIndex);
}

int DataPH::setOfPartialHaps::getNb(int phIndex) { return phVec[phIndex].getNb(); }
int DataPH::setOfPartialHaps::getNa(int phIndex) { return phVec[phIndex].getNa(); }

vector<vector<Sequence> >  DataPH::findMMS(std::vector<Sequence> & MSfhs, int seed) {
	return data.findMMS(MSfhs, seed);
}


//Method for finding the set of minimal-minimal sets of full haplotypes, i.e. the set containing
//sets of the smallest number of full haplotypes possible to explain the partial haplotype data
//
//E.g. if the true full haplotypes are ACT, AGT, CCG
//the partial haplotypes will be something like
//AC-  -CT  A--  -C-  --G
//AG-  -GT  C--  -G-  --G
//CC-  -CG
//
//from which one may deduce a minimal set of full haplotypes:
//ACT, AGT, ACG, CCT, CCG
//
//This method then finds the set of  minimal-minimal sets of full haplotypes:
//{ACT,AGT,CCG} and {AGT,ACG,CCT} (in this case) which both explain the partial haplotypes,
//however, only the former is the correct set. The MMS will have length <= the length of the MS.
vector<vector<Sequence> >  DataPH::phDataSingleSim::findMMS(std::vector<Sequence> & MSfhsInput, int seed) {


	vector<Sequence> MSfhs = MSfhsInput;	
	//Randomise MSfhs to speed up algorithm
	srand(seed); //Use built in random_shuffle rather than gsl rng
	random_shuffle(MSfhs.begin(), MSfhs.end());	
	computeContribs(MSfhs); //Updates contribs wrt full haplotypes

	vector<vector<Sequence> > result;

	//Use greedy set cover method

	//First, assign contributions to all full haplotypes
	vector<vector<int> > fhContribs(MSfhs.size());
	int phCounter = 0; //Counts the number of partial haplotypes we loop over
	for(int i=0; i<numOfPHSets(); i++) {

		for(int j=0; j<getNumOfPartialHaps(i);j++) {
			vector<int> phContribs = *(getContribs(i,j)); //Get the full haplotypes the current ph corresponds to
			for(unsigned int k=0; k<phContribs.size(); k++) { //Loop over the full haplotypes
				fhContribs[phContribs[k]].push_back(phCounter);
			}
			phCounter ++; //Next partial haplotype
		}
	}
	cout << "full hap contribs assigned.\n";


	//Next, perform greedy algorithm, i.e. at each step, choose haplotype with most contributions in remainingPHs
	//Keep doing this until cover is achieved (i.e. all the number from 0 to phCounter-1)
	vector<Sequence> greedySet;
	vector<int> fullSetPHs; //The full set of partial haplotype IDs (e.g. numbers 0 to phCounter-1)
	for(int i=0; i<phCounter;i++) {
		fullSetPHs.push_back(i);
	}
	vector<int> pickedPHs; //The IDs of the partial haplotypes that have been picked (through picking full haplotypes)
	vector<int> remainingPHs = fullSetPHs; //The set of ph IDs we haven't yet covered
	bool cover = false; //Check whether we have full cover. Start with false
	while(cover == false) {

		int fullHapWithMostContribsInRemainingPHs = -1; //Need to find the ID of the full hap with most contribs in remaining
		int mostContribs = 0;
		for(unsigned int i=0; i<fhContribs.size(); i++) {
			vector<int> currentContribs = fhContribs[i];
			int currentContribsCount = 0;
			for(unsigned int j=0; j<currentContribs.size(); j++) { //Loop over contribution from current full hap
				for(unsigned int k=0; k<remainingPHs.size(); k++) { //Compare each of them to the ones in the remaining set
					if(currentContribs[j] == remainingPHs[k]) { //If they are the same
						currentContribsCount++; //Up the count
						break;
					}
				}
			}
			if(currentContribsCount > mostContribs) { //Check if this full haplotype has the most cover in remaining
				fullHapWithMostContribsInRemainingPHs = i;
				mostContribs = currentContribsCount;
			}
		}

		greedySet.push_back(MSfhs[fullHapWithMostContribsInRemainingPHs]);
		vector<int> pickedContribs = fhContribs[fullHapWithMostContribsInRemainingPHs];
		//Update remaining and picked
		cout << "New full hap pick with " << pickedContribs.size() << " contribs\n";
		for(unsigned int i=0; i<pickedContribs.size(); i++) {

			//Check if contrib i has already been covered
			bool alreadyPicked = false;
			for(unsigned int j=0; j<pickedPHs.size(); j++) {
				if(pickedPHs[j] == pickedContribs[i]) {
					alreadyPicked = true; break;
				}
			}
			if(alreadyPicked == false) { //add it to picked and remove it from remaining
				pickedPHs.push_back(pickedContribs[i]);
				//sort(pickedPHs.begin(), pickedPHs.end()); //Sort vector (ascending order)
				//reverse(pickedPHs.begin(),pickedPHs.end()); //Reverse to get descending order

				cout << "pickedPH = " << pickedContribs[i] << "\n";
				cout << "remainingPHs before: \n";
				printIntVector(remainingPHs);
				for(unsigned int k=remainingPHs.size()-1; k>=0; k--) {
					if(remainingPHs[k] == pickedContribs[i]) {
						remainingPHs.erase(remainingPHs.begin() + k); //This is fine as we are doing in descending order
						break;
					}
				}
				cout << "remainingPHs after: \n";
				printIntVector(remainingPHs);
			}
		}

		//Check cover
		cover = true;
		sort(pickedPHs.begin(),pickedPHs.end()); //Sort ascending order
		cout << "cover check:\n";
		for(int i=0; i<phCounter; i++) {
			if(pickedPHs[i] == i) {
				cout << "Match: " << i << "\n";
			}
			if(pickedPHs[i] != i) { //Can do this because we sorted the vector pickedPHs
				cout << "Mismatch: " << pickedPHs[i] << " " << i << "\n";
				cover = false; break;
			}	
		}

	}




	//if(MSfhs.size() < 30) { //Do exhaustive method if not too involved
	if(1==2) { // Don't do this bit
		//if(1==1) { // Always do this bit }


		//Vector containing the fhs removed at the various levels, e.g. {3,5,4} means
		//FH3 has been removed at level 0, FH5 at level 1 and FH4 at level 2
		vector<int> indicesRemoved = {-1}; 
		int level = 0; //Indicating which level we are currently at
		vector<Sequence> removedSeqs = {MSfhs[0]}; //Set it to a Sequence, which is irrelevant as it will be overwritten
		bool newLevel = false; //Indicated whether a shift in level has just occurred
		bool endOfLevel = false;
		int shortestMMS = (int) greedySet.size(); 

		//Keep finding sets as long as we haven't checked the last index on the lowest level
		//(We start at the first index and then work the way down the levels, then back up again and go to the second index, etc.
		while(indicesRemoved[0] != (int) MSfhs.size()-1) {

			cout << "In while loop...\n";
			printIntVector(indicesRemoved);
			if(newLevel == true) { //The value of level has just been adjusted

				if((int) indicesRemoved.size() >= level +2) { //Remove the last n levels

					//cout << "A\n";
					while((int) indicesRemoved.size() > level +1) {
						indicesRemoved.pop_back();
						removedSeqs.pop_back();
					}
					//cout << "B\n";	
					//Get the next index - returns -1 if none left
					int nextLeft = getNextLeft(indicesRemoved, (int) MSfhs.size() -1);
					//cout << "C\n";
					if(nextLeft == -1) {
						endOfLevel = true; 
					} else {
						indicesRemoved[level] = nextLeft;
						removedSeqs[level] = MSfhs[nextLeft];
					}


					//if(level <=3) { printIntVector(indicesRemoved); }
				}

				if((int) indicesRemoved.size() == level) { //Add a new level
					//cout << "Going down a level: \n";	
					indicesRemoved.push_back(indicesRemoved[level-1]);
					int nextLeft = getNextLeft(indicesRemoved, (int) MSfhs.size() -1);
					//cout << "nextLeft: " << nextLeft << "\n";

					if(nextLeft == -1) {
						endOfLevel = true;
					} else {				
						indicesRemoved[level] = nextLeft;
						removedSeqs.push_back(MSfhs[nextLeft]);
						//MSfhs[nextLeft].print();
					}
				}

				newLevel = false;

			} else { //Stay at the same level

				int nextLeft = getNextLeft(indicesRemoved, (int) MSfhs.size() -1);
				//cout << "NextLeft in stay at same level: " << nextLeft << "\n";
				if(nextLeft == -1) {
					endOfLevel = true;
				} else {
					indicesRemoved[level] = nextLeft;
					removedSeqs[level] = MSfhs[nextLeft];
				}

			}

			//We have reached the end of a certain level
			if(endOfLevel == true) {

				level--;
				newLevel = true;
				endOfLevel = false;

			} else {
				//Check if this set is a MMS
				vector<Sequence> possibleMMS;
				for(unsigned int i=0; i<MSfhs.size(); i++) {

					bool includeHap = false;
					for(unsigned int j=0; j<indicesRemoved.size(); j++) {
						if((int) i==indicesRemoved[j]) { //Don't include this fh
							includeHap = true;
							break;
						}
					}

					if(includeHap == true) { possibleMMS.push_back(MSfhs[i]); }
				}	

				phDataSingleSim possibleMMSphdss; //Empty data holder
				vector<Sequence> truePHSeqs = getSequences(); //The actual partial haplotypes from dataset
				decomposeFullHaps(possibleMMS,possibleMMSphdss, truePHSeqs);

				//Check if the true partial haplotypes are contained within the decomposed set
				vector<Sequence> possibleMMSphSeqs = possibleMMSphdss.getSequences();
				bool MMS = true;
				for(unsigned int i=0; i<truePHSeqs.size(); i++) {

					bool hapFound = false;
					for(unsigned int j=0; j<possibleMMSphSeqs.size(); j++) {
						hapFound = identicalSeqs(possibleMMSphSeqs[j], truePHSeqs[i]);
						if(hapFound == true) {
							break;
						}
					}
					if(hapFound == false) { 
						MMS = false;
						break;
					}
				}

				if(MMS == true) { //Add the MMS and stay at currennt level

					result.push_back(possibleMMS);
					//cout << "Added length: " << possibleMMS.size() << "\n";
					int min = possibleMMS.size();

					//If the new MMS is shorter than the previous, remove the previous entries
					if(min < shortestMMS) {
						result.clear();
						result.push_back(possibleMMS);
						shortestMMS=min;
						cout << "New value for shortest MMS: " << min << "\n";
						newLevel = false; //Stay at same level and check for other solutions of same length
					} else if(indicesRemoved[level] == (int) MSfhs.size()-1) { //Cannot go down any further
						//So go up
						level--;
						newLevel = true;
					} else {//Else stay where we are
						newLevel = false;
					}
				} else {

					if(indicesRemoved[level] == (int) MSfhs.size()-1) { //Cannot go down any further
						//So go up
						level--;
						newLevel = true;
					} else if (level+1 == shortestMMS) { //We only want solutions of the same length or shorter than the shortestMMS
						//Hence, if level+1 is of the same length as the best solution, stay at the same level
						newLevel = false; 
					} else{ //Go a level down to search for more advane solutions
						level++;
						newLevel = true;
					}
				}
			}
		}

	} else { //Else just go with greedy result

		result.push_back(greedySet);

	}

	return result;
}

int DataPH::getNbTot() { return data.getNbTot(); }
int DataPH::getNaTot() { return data.getNaTot(); }
int DataPH::getNtot() { return getNbTot() + getNaTot(); }
int DataPH::phDataSingleSim::getNbTot() { 

	int result = 0;
	for(unsigned int i=0; i<dataset.size(); i++) {
		result += dataset[i].getNbTot();
	}
	return result;
}
int DataPH::phDataSingleSim::getNaTot() { 

	int result = 0;
	for(unsigned int i=0; i<dataset.size(); i++) {
		result += dataset[i].getNaTot();
	}
	return result;
}


bool DataPH::decomposeAndCheck(vector<Sequence> &fhs) {
	return data.decomposeAndCheck(fhs);
}

bool DataPH::phDataSingleSim::decomposeAndCheck(vector<Sequence> &fhs) {

	phDataSingleSim possibleMMSphdss; //Empty data holder
	vector<Sequence> truePHSeqs = getSequences(); //Actual partial haplotypes from dataset
	decomposeFullHaps(fhs,possibleMMSphdss,truePHSeqs);

	//Check if the true partial haplotypes are contained within the decomposed set
	vector<Sequence> possibleMMSphSeqs = possibleMMSphdss.getSequences();
	bool MMS = true;
	for(unsigned int i=0; i<truePHSeqs.size(); i++) {

		bool hapFound = false;
		for(unsigned int j=0; j<possibleMMSphSeqs.size(); j++) {
			hapFound = identicalSeqs(possibleMMSphSeqs[j], truePHSeqs[i]);
			if(hapFound == true) {
				break;
			}
		}
		if(hapFound == false) { 
			MMS = false;
			break;
		}
	}
	return MMS;
}

//Method for writing the partial haplotype data to a file
//Output of form:
//[Haplotype BeforeCount AfterCount]
//e.g.
//ATA 10 5
//ACG 15 20
//AT- 50 45
//AC- 60 65
//...
void DataPH::writeToFile(ofstream &outputFile) { data.writeToFile(outputFile); }
void DataPH::writeToFileMahan(ofstream &outputFile, vector<int>& physicalPos) { data.writeToFileMahan(outputFile, physicalPos); }
void DataPH::writeToFileHaps(string &filePathGeneFolder, vector<int>& physicalPos) { data.writeToFileHaps(filePathGeneFolder, physicalPos); }
void DataPH::phDataSingleSim::writeToFile(ofstream &outputFile) {
	for(unsigned int i=0;i<dataset.size();i++) {
		dataset[i].writeToFile(outputFile);
	}
}
void DataPH::phDataSingleSim::writeToFileMahan(ofstream &outputFile, vector<int>& physicalPos) {
	for(unsigned int i=0;i<dataset.size();i++) { //Loop over sophs
		dataset[i].writeToFileMahan(outputFile, physicalPos);
	}
}
void DataPH::phDataSingleSim::writeToFileHaps(string &filePathGeneFolder, vector<int>& physicalPos) {
	for(unsigned int i=0;i<dataset.size();i++) { //Loop over sophs. Assumes they are ordered by number of loci covered, large to small.
		dataset[i].writeToFileHaps(filePathGeneFolder, physicalPos, i);
	}
}

void DataPH::setOfPartialHaps::writeToFile(ofstream &outputFile) {

	for(unsigned int i=0; i<phVec.size(); i++) {
		phVec[i].writeToFile(outputFile);
	}

	outputFile << "Total\t" <<  NbTot << "\t" << NaTot << "\n";
}
void DataPH::setOfPartialHaps::writeToFileMahan(ofstream &outputFile, vector<int>& physicalPos) {

	//Find positions for this specific soph
	vector<int> physPosSoph;
	for(unsigned int i=0; i<lociCovered.size(); i++) {
		physPosSoph.push_back(physicalPos[lociCovered[i]]);
	}

	for(unsigned int i=0; i<phVec.size(); i++) {
		phVec[i].writeToFileMahan(outputFile, physPosSoph);
	}

}
void DataPH::setOfPartialHaps::writeToFileHaps(string &filePathGeneFolder, vector<int>& physicalPos, int sophIndex) {

	//Create Hap_dataX.dat file
	stringstream ss;
	ss << filePathGeneFolder << "Hap_data" << sophIndex << ".dat";
	string filePath = ss.str();
	ofstream outputFile;
	outputFile.open(filePath.c_str());

	for(unsigned int i=0; i<phVec.size(); i++) {

		if(i!=0) { outputFile << "\n"; }
		outputFile << phVec[i].getSeq().printToStringNoHyphens() << " " << phVec[i].getNb() << " " << phVec[i].getNa(); 
	}
	outputFile.close();

	//Create LociX.dat file
	ss.str("");
	ss << filePathGeneFolder << "Loci" << sophIndex << ".dat";
	filePath = ss.str();
	outputFile.open(filePath.c_str());

	for(unsigned int i=0; i<lociCovered.size(); i++) {

		outputFile << lociCovered[i];
		if(i!=lociCovered.size()-1) {
			outputFile << " ";
		}
	}
	outputFile.close();

//// Not sure what I was trying to do with the below??
//
//	//Find positions for this specific soph
//	vector<int> physPosSoph;
//	for(unsigned int i=0; i<lociCovered.size(); i++) {
//		physPosSoph.push_back(physicalPos[lociCovered[i]]);
//	}
//
//
//	//Print physical position for this soph
//	ss.str("");
//	ss << filePathGeneFolder << "Positions.dat";
//	filePath = ss.str();
//	outputFile.open(filePath.c_str());
//	for(unsigned int i=0; i<physPosSoph.size(); i++) {
//		
//		outputFile << physPosSoph[i];
//		if(i!=physPosSoph.size()-1) {
//			outputFile << " ";
//		}
//	}
	outputFile.close();

}
void DataPH::partialHap::writeToFile(ofstream &outputFile) {

	for(int i=0; i<seq.getLength(); i++) {
		outputFile << seq.getBase(i);
	}

	outputFile << "\t" <<  Nb << "\t" << Na << "\n";
}
void DataPH::partialHap::writeToFileMahan(ofstream &outputFile, vector<int>& physPosPH) {

	//Find number of non-hyphen loci in ph
	int numNonHyphenLoci = 0;
	for(int i=0; i<seq.getLength(); i++) {
		if(seq.getBase(i) != '-') { numNonHyphenLoci++; }
	}
	outputFile << numNonHyphenLoci << "\t";

	//Write physical pos for this PH
	for(unsigned int i=0; i<physPosPH.size(); i++) {
		outputFile << physPosPH[i] << "\t";
	}

	//Print partial hap without hyphens
	for(int i=0; i<seq.getLength(); i++) {
		if(seq.getBase(i) != '-') {
			outputFile << seq.getBase(i);
		}
	}
	
	//Print number of time points (2 - before and after)
	outputFile << "\t2\t";

	//Print timepoint and number of reads
	outputFile << "0\t" <<  Nb << "\t1\t" << Na << "\n";
}
double DataPH::phDataSingleSim::computeDegreeOfOverlapBefore() {

	int weightedTotal = 0;
	int totalNumberOfReads = 0;

	for(unsigned int i=0; i<dataset.size(); i++) {

		dataset[i].computeDegreeOfOverlapBefore(weightedTotal, totalNumberOfReads);
	}


	double dooBefore = weightedTotal/((double) totalNumberOfReads);

	return dooBefore;

}

void DataPH::setOfPartialHaps::computeDegreeOfOverlapBefore(int& weightedTotal, int& totalNumberOfReads) {

	for(unsigned int i=0; i<phVec.size(); i++) {

		//Get current partial hap
		DataPH::partialHap ph = phVec[i];
		Sequence seq = ph.getSeq();
		int Nb = ph.getNb();

		//Get the number of loci covered by partial hap
		int phCoverage = 0;
		for(int j=0; j<seq.getLength(); j++) {

			if(seq.getBase(j) != '-') { phCoverage++; }
		}

		weightedTotal += Nb*phCoverage;
		totalNumberOfReads += Nb;

	}

}

double DataPH::phDataSingleSim::computeDegreeOfOverlapAfter() {

	int weightedTotal = 0;
	int totalNumberOfReads = 0;

	for(unsigned int i=0; i<dataset.size(); i++) {

		dataset[i].computeDegreeOfOverlapAfter(weightedTotal, totalNumberOfReads);
	}

	double dooAfter = weightedTotal/((double) totalNumberOfReads);

	return dooAfter;

}

void DataPH::setOfPartialHaps::computeDegreeOfOverlapAfter(int& weightedTotal, int& totalNumberOfReads) {

	for(unsigned int i=0; i<phVec.size(); i++) {

		//Get current partial hap
		DataPH::partialHap ph = phVec[i];
		Sequence seq = ph.getSeq();
		int Na = ph.getNa();

		//Get the number of loci covered by partial hap
		int phCoverage = 0;
		for(int j=0; j<seq.getLength(); j++) {

			if(seq.getBase(j) != '-') { phCoverage++; }
		}

		weightedTotal += Na*phCoverage;
		totalNumberOfReads += Na;
	}
}

void DataPH::setOfPartialHaps::setLociCovered(vector<int>& lc) { lociCovered = lc; }
vector<int> DataPH::setOfPartialHaps::getLociCovered() { return lociCovered; }

void DataPH::phDataSingleSim::removeLocus(int locus) {

	for(unsigned int i=0; i<dataset.size(); i++) {
		dataset[i].removeLocus(locus);
	}
}

void DataPH::setOfPartialHaps::removeLocus(int locus) {

	for(unsigned int i=0; i<phVec.size(); i++) {
		phVec[i].removeLocus(locus);
	}

	//Also remove locus from lociCovered and adjust accordingly,
	//i.e. if there are 6 loci (0,1,2,3,4,5) and loci (1,2,3,4) are covered
	//and locus 3 is removed, then loci covered becomes (1,2,3)
	int toBeErased = -1;
	for(unsigned int i=0; i<lociCovered.size(); i++) {

		if(lociCovered[i] == locus) {

			toBeErased = i;

		} else if(lociCovered[i] > locus) {

			lociCovered[i]--; //Adjust loci positions larger than removed locus
		}
	}
	if(toBeErased != -1) {
		lociCovered.erase(lociCovered.begin()+toBeErased); //Remove locus from lociCovered
	}
	cout << "Printing loci covered temporarily: "; printIntVector(lociCovered);
}

void DataPH::partialHap::removeLocus(int locus) {

	seq.removeBase(locus);
}

vector<int> DataPH::getPhysicalPos() { return physicalPos; }
void DataPH::clearPhysicalPos() { physicalPos.clear(); }

vector<Sequence>& DataPH::getImportedFullHaps() { return importedFullHaps; }
void DataPH::setImportedFullHaps(vector<Sequence>& ihf) { importedFullHaps = ihf; }

//Given imported haplotypes, we need to update the data to avoid any issues. We have to look out for two things:
//1) Some loci in the full haplotypes are monomorphic. This means that this locus needs 
//deleting from all partial haps. Need to update physical pos as well, if present.
//2) Some partial haps may no longer be represented by a full haplotype. In this case we
//need to shorten the haplotype until it can be merged with a shorter haplotype.
void DataPH::updateDataWithRespectToImportedFullHaps() {


	//First, we check if there are monomorphic loci.
	vector<int> monoLoci;
	for(int i=0; i<importedFullHaps[0].getLength(); i++) { //Loop over loci

		char refAllele = importedFullHaps[0].getBase(i);
		bool monomorphic = true;
		for(unsigned int j=0; j<importedFullHaps.size(); j++) { //Loop over full haplotypes
		
			if(importedFullHaps[j].getBase(i) != refAllele) {
				monomorphic = false; break;
			}
		}

		if(monomorphic == true) { monoLoci.push_back(i); }

	}

	//Then we remove the monomorphic loci one at a time, correcting the data as we go along
	for(int i= ((int) (monoLoci.size())) -1; i>=0; i--) { //Go through loci in reverse to avoid fuck-ups

		cout << "Index i in monoLoci loop: " << i <<"\n";

		//First correct the full haps
		for(unsigned int j=0; j<importedFullHaps.size(); j++) {

			importedFullHaps[j].removeBase(monoLoci[i]);
			if(j==0) {
				cout << "Removing loci " << monoLoci[i] << " from full haplotypes.\n";
			}
		}


		//Then correct partial haplotype data
		cout << "Removing loci " << monoLoci[i] << " from partial haps.\n";
		data.removeLocus(monoLoci[i]);
		cout << "Merging partial haplotype sets:\n";
		data.mergeHaps();
		cout << "Ordering partial haplotype sets:\n";
		data.orderHaps();


		//Then correct physical pos if existing
		if(physicalPos.size() > 0) {
//			cout << "Correcting physical pos.\n";
			physicalPos.erase(physicalPos.begin()+monoLoci[i]);
		}

	}

	cout << "Done with monomorphic loci.\n";

	
	//Next, we check that all partial haplotypes are represented by some full haplotype
	//The following command loops through all partial haplotypes and checks that they are
	//represented by a full haplotype. If this is not the case, it is shortened and placed
	//in new ph set.
	if(data.getNumDatasets()>0) { //Only if any data left
		cout << "Check whether all PHs are represented by imported full haps. If not, they get subsetted until they do.\n";
		data.ensurePHsInImportedHaps(importedFullHaps);
		cout << "Merging partial haplotype sets to make up for any new sophs created by the previous process.\n";
		data.mergeHaps();
		cout << "Ordering partial haplotype sets to make up for any disorder created by previous processes.\n";
		data.orderHaps();


		//Update contribs
		cout << "Computing contribs for partial haps with respect to imported full haps.\n";
		data.computeContribs(importedFullHaps);
	}
}

//Checks that all partial haplotypes are represented by some full haplotype. This is sometimes not the case.
//If this is not the case, then the partial haplotype is shortened to the largest possible partial haplotype
//that is represented by a full haplotype. In the worst case, the partial haplotype is entirely discarded.
void DataPH::phDataSingleSim::ensurePHsInImportedHaps(vector<Sequence>& fhs) {

	//Loop over sets of partial haps
	for(int i=((int)dataset.size())-1; i>=0; i--) { //Go through in reverse so we can remove SOPHS without messing up index i

		cout << "i: " << i << "\n";	
		//Get the data
		vector<Sequence> phs = dataset[i].getSequences();
		vector<int> Nb = dataset[i].getNb();
		vector<int> Na = dataset[i].getNa();
		vector<int> lociCovered = dataset[i].getLociCovered();


		//Go through each ph in turn
		for(int j=((int) phs.size())-1; j>=0; j--) {

			cout << "j: " << j << "\n";	
			Sequence ph = phs[j];
			cout << "Sequence j: "; ph.print();
			cout << "lociCovered: "; printIntVector(lociCovered);
			Sequence subsetPH = ph.subset(lociCovered[0],lociCovered[0]+lociCovered.size()); //Subset gives [start,finish)

		
			bool phRepresentedByFH = false;

			//Check if represented by a full haplotype
			for(unsigned int k=0; k<fhs.size(); k++) { //Loop over imported full haps

				//E.g. if loci covered are 2,3,4, then get subset [2,4]
				Sequence subsetFH = fhs[k].subset(lociCovered[0],lociCovered[0]+lociCovered.size()); //Subset gives [start,finish)

				bool compare = identicalSeqs(subsetFH,subsetPH);
				if(compare == true) {

					phRepresentedByFH = true;
					break;
				}
			}

			if(phRepresentedByFH == false) {

				cout << "SOPH " << i << " partial haplotype " << j << " is not represented by a full haplotype:\n";
				ph.print();

				//Next we wish to find all the possible partial haplotype subsets of ph.
				//A partial haplotype subset is an unbroken subset, i.e. --XXX--- is a
				//subset of --XXXX--, but --X-XX-- is not.
				vector<vector<int> > subsets;
				for(unsigned int removedFromLeft=0; removedFromLeft<lociCovered.size(); removedFromLeft++) {
					for(unsigned int removedFromRight=0; removedFromRight<lociCovered.size(); removedFromRight++) {

						//Try removing loci from left and right, but only as long as we don't remove them all.
						//There needs to be at least one loci remaining.
						if(removedFromLeft + removedFromRight < lociCovered.size() && removedFromLeft + removedFromRight > 0) {
							vector<int> subset = lociCovered; //Set subset equal to original list of covered loci
							subset.erase(subset.end()-removedFromRight,subset.end()); //Remove from right
							subset.erase(subset.begin(),subset.begin()+removedFromLeft); //Remove from left

							subsets.push_back(subset);
						}
					}
				}

				//Order subsets by size, largest first
				sort(subsets.begin(), subsets.end(), compareBySize);

				//Check if any of the new subsets are working, i.e. that the subsetted partial haplotype can be represented by
				//some subsetted full haplotype.
				int workingSubset = -1;
				for(unsigned int k=0; k<subsets.size(); k++) {
		
					Sequence subsetPH = ph.subset(subsets[k][0],subsets[k][0]+subsets[k].size()); //Subset gives [start,finish)
					for(unsigned int l=0; l<fhs.size(); l++) {
						Sequence subsetFH = fhs[l].subset(subsets[k][0],subsets[k][0]+subsets[k].size()); // Subset gives [start,finish)

						//Compare the subsetted partial and full haplotype
						//If identical, choose this as the working subset.
						//This might create a slight bias in that it chooses subsets
						//of the form XXX-- before e.g. --XXX due to the way the subsets
						//have been generated.
						bool compare = identicalSeqs(subsetFH,subsetPH);
						if(compare == true) {
				
							workingSubset = k;
							break;
						}
					}

					if(workingSubset != -1) { break; } //The first subset we find will be the longest, so we stick with that one
				}

				//If a working subset is found, transfer haplotype j to a new soph covering the loci in the working subset.
				if(workingSubset != -1) {

					cout << "A working subset was found: "; printIntVector(subsets[workingSubset]);
					cout << "The haplotype has been transferred.\n";

					//Get the (working) subset of the PH
					Sequence newPHseq = ph.subset(subsets[workingSubset][0],subsets[workingSubset][0]+subsets[workingSubset].size()); //Subset gives [start,finsh)
					
					//Add '-' before and after accordingly, i.e. make into the length of a full haplotype
					for(int k=0; k<subsets[workingSubset][0]; k++) {

						newPHseq.addBaseFront('-');
					}
					for(int k=subsets[workingSubset][0]+subsets[workingSubset].size(); k<fhs[0].getLength(); k++) {

						newPHseq.addBase('-');
					}
					
					//Create the new partial haplotype with number of reads from the old
					DataPH::partialHap newPH;
					newPH.setSeq(newPHseq);	
					newPH.setNb(Nb[j]);
					newPH.setNa(Na[j]);
					
					//Add to set of partial haps
					DataPH::setOfPartialHaps newSOPH;
					newSOPH.addPartialHap(newPH);
					newSOPH.countNbTot();
					newSOPH.countNaTot();
					newSOPH.setLociCovered(subsets[workingSubset]);

					//Add to data - to be merged with possible identical sophs later
					dataset.push_back(newSOPH);
					
				} else {

					cout << "No working subset found. Haplotype has been deleted.\n";
				}	
				


				//Remove the old partial haplotype. Also the case if no working subset is found
				dataset[i].removePH(j);

				//Update the old SOPH
				if(dataset[i].getNumOfPartialHaps() > 0) {
					
					//Update Nb and Na
					dataset[i].countNbTot();
					dataset[i].countNaTot();

				} else { //No more haplotypes left, so need to delete soph

					dataset.erase(dataset.begin() + i);
						
				}
					
			}
		}
	}
}



//Full haplotype variables -- for debugging use only
vector<Sequence> DataPH::getFullHapsTrue() { return fullHapsTrue; }
vector<double> DataPH::getqBtrue() { return qBtrue; }


//Reads
void DataPH::pairedEndRead::generateRead(double meanReadLength, double stDevReadLength, double meanGapLength, double stDevGapLength, int length, gsl_rng* r, double C, vector<double>& freqs) {

	//Length is the length of the gene, e.g. [1,1270]
	bool withinGene = false;

	while(withinGene == false) { //Keep generating reads until one is within the gene limits

		start1 = gsl_rng_uniform_int(r,length) + 1; //[1,length]
		end1 = -1;

		//		while(end1 <= start1) { //normal variate can be < 0
		end1 = start1 + gsl_ran_gaussian(r, stDevReadLength) + meanReadLength;
		//		}

		start2 = -1;
		//		while(start2 <= end1) {
		start2 = end1 + gsl_ran_gaussian(r, stDevGapLength) + meanGapLength;
		//		}

		end2 = -1;
		//		while(end2 <= start2) {
		end2 = start2 + gsl_ran_gaussian(r, stDevReadLength) + meanReadLength;
		//		}


		//		if(end2 <= length) { withinGene = true; }


		//Check if start1, end1 and start2, end2 forms a proper read, i.e. a read where
		// 0 < start1 < end1 < start2 < end2 <= length
		if(start1 < end1 && end1 < start2 && start2 < end2 && end2 <= length) {

			withinGene = true;
		}
	}

	//Generate haplotype index, i.e. integer on [0,freqs.size()) from an overdispersed sampling method
	vector<int> draw = DirMultSampling(1,freqs,C,r); //Get DirMult sample with N=1
	for(unsigned int i=0; i<draw.size(); i++) {
	
		if(draw[i] == 1) { 
			
			hapID = i;
			break;
		}
	}

}


//Check if a position is covered by read of the form
//-------|start1 xxxxxxxxxxxxxx end1|----(gap)-----|start2 xxxxxxxxxxxxxx end2|-----------
bool DataPH::pairedEndRead::posCovered(int& pos) {


	if(pos < start1) {
		return false;
	} else if(pos > end1 && pos < start2) {
		return false;
	} else if(pos > end2) {
		return false;
	}

	//If position is inside read, return true
	return true;
}


//Check if a vector of positions are covered by read of the form
//-------|start1 xxxxxxxxxxxxxx end1|----(gap)-----|start2 xxxxxxxxxxxxxx end2|-----------
bool DataPH::pairedEndRead::posCovered(vector<int>& pos) {

	//For each position in vector, check if it is covered by read, if not return false
	for(unsigned int i=0; i<pos.size(); i++) {

		if(posCovered(pos[i]) == false) { return false; }

	}

	//If no positions are outside read, return true
	return true;
}

void DataPH::pairedEndRead::print() {

	cout << "Start1: " << start1 << "\tEnd1: " << end1 << "\tStart2: " << start2 << "\tEnd2: " << end2 <<"\n";
}

int DataPH::pairedEndRead::getHapID() { return hapID; }

vector<int> DataPH::setOfPartialHaps::getPosCovered(vector<int>& allPos) {

	Sequence seq0 = phVec[0].getSeq();
	int seqLength = seq0.getLength();
	vector<int> result;

	if(seqLength == (int) allPos.size()) {

		for(unsigned int i=0; i<allPos.size(); i++) {

			if(seq0.getBase(i) != '-') { //Locus covered by this set of partial haps

				result.push_back(allPos[i]);
			}
		}

	} else {

		cout << "Error in DataPH::setOfPartialHaps::getPosCovered: input vector not same size as partial haplotype sequence.\n";
	}

	return result;
}

void DataPH::phDataSingleSim::removeSOPH(int index) {

	dataset.erase(dataset.begin() + index);
}

int DataPH::phDataSingleSim::getNumberOfSites() {

	return dataset[0].getNumberOfSites();
}

int DataPH::setOfPartialHaps::getNumberOfSites() {

	return phVec[0].getNumberOfSites();
}

int DataPH::partialHap::getNumberOfSites() {

	return seq.getLength();

}

bool DataPH::getImportFullHaps() { return importFullHaps; }

int DataPH::getDeltaDays() { return deltaDays; }

//Minor freqs
vector<double> DataPH::getDonorSLFreqs() {

	vector<double> donorFreqs;
	for(unsigned int i=0; i<slTrajs.size(); i++) { //Loop over trajectories

		donorFreqs.push_back(slTrajs[i].getNbMinor()/((double) slTrajs[i].getNbTot()));
	}
	return donorFreqs;
}

//Minor freqs
vector<double> DataPH::getRecipientSLFreqs() {

	vector<double> recipientFreqs;
	for(unsigned int i=0; i<slTrajs.size(); i++) { //Loop over trajectories

		recipientFreqs.push_back(slTrajs[i].getNaMinor()/((double) slTrajs[i].getNaTot()));
	}
	return recipientFreqs;
}

//Minor obs
vector<int> DataPH::getRecipientSLObs() {

	vector<int> recipientObs;
	for(unsigned int i=0; i<slTrajs.size(); i++) { //Loop over trajectories

		recipientObs.push_back(slTrajs[i].getNaMinor());
	}
	return recipientObs;
}

//NaTot
vector<int> DataPH::getRecipientSLObsTot() { 

	vector<int> recipientObsTot;
	for(unsigned int i=0; i<slTrajs.size(); i++) { //Loop over trajectories

		recipientObsTot.push_back(slTrajs[i].getNaTot());
	}
	return recipientObsTot;
}

