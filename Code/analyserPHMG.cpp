//
//  analyserPHMG.cpp
//  FluTransmissionProject
//
//  Created by Casper Lumby on 08/12/2015.
//  Copyright © 2015 Casper Lumby. All rights reserved.
//

#include "analyserPHMG.hpp"
#include <iostream>
#include "misc.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <limits>
#include <math.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "model.hpp"

#ifdef _OPENMP
	#include <omp.h>
#else
	#define omp_get_thread_num() 0
#endif

using namespace std;

AnalyserPHMG::AnalyserPHMG() { } //Constructor does nothing
AnalyserPHMG::~AnalyserPHMG() { } //Destructor does nothing

//Data *d pointer passed by copy, not by reference. Can't change original data.
void AnalyserPHMG::loadData(Data *d) {
    
	cout << "Loading data for full haplotypes. \n";
	data = dynamic_cast<DataPHMG *>(d);
}
 
void AnalyserPHMG::runAnalysis(AnaParam *ap) {

	//Check if any data all
	if(data->getDataPresent() == false) {

		cout << "No data is present for this transmission event at the time of analysis. This may be anticipated. If not, investigate. Exiting.\n";
		ofstream outputFile;
		string filePathFolder = ap->getOutputFolder();
		stringstream ss;
		ss << filePathFolder << "noDataError.dat";
                string filePath = ss.str();
                outputFile.open(filePath.c_str());
		outputFile << "No data was present for this transmission event at the time of analysis. This may be anticipated. If not, investigate.\n";
		outputFile.close();
		exit(1);
	}
    
	cout << "Running analysis for multiple genes partial haplotypes.\n";

	//Set up RNG
	gsl_rng * r;
	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_set(r, ap->getSeed());


	//Setup output folder
	ofstream outputFile;
	string filePathFolder = ap->getOutputFolder();
	stringstream ss;
	string filePath;
						
	
	/*
	 * Case where we filter haplotypes within the C++ code (in contrast to A) using all the possible haplotypes that could explain the dataset
	 * or B) importing pre-filtered haplotypes). 
	 *
	 * The below process is for the filter 1 scenario:
	 *
	 * 	- Generate all possible haplotypes
	 * 	- Infer before frequencies only (i.e. we don't infer a covariance matrix!)
	 *	- Remove haplotypes with frequencies <= 2*10^-10
	 *	- If there are more than 100 haplotypes remaining, increase cutoff by factor of 10
	 *	- Keep doing this until < 100 haplotypes 
	 *	- Pass these along as "imported haplotypes" further downstream
	 *
	 * NB: This is an early version of this method. It works, but code may need to be cleaned up slightly
	 * and outcomes verified before official usage. Filter 3 is the currently used method for filtering.
	 */
	vector<vector<Sequence> > filteredHaplotypes; //May never be initiated!
	double cutOff = 2e-10;
	double C = ap->getC();
	if(ap->getFilterHaplotypesMethod() == 1) { //Filter method 1 triggered

		for(int g=0; g<getNumGenes(); g++) {

			vector<Sequence> fullHapsGene;
			if(data->getGenePresent(g) == true) {
				DataPH* d = data->getGene(g);

				//Generate full haps from partial haps
				vector<Sequence> phs = d->getSequences(); //Get partial haplotypes
			
				AnalyserPH aPH;
				aPH.loadData(d);
				aPH.constructFullHaps(&phs, &fullHapsGene);
				vector<Sequence> fullHapsGeneOriginal = fullHapsGene; //Make a copy
				d->computeContribs(fullHapsGene); //Updates contribs wrt full haplotypes

				cout << "Number of haplotypes before filtering for gene " << g << " : " << fullHapsGene.size() << "\n";


				//Infer frequencies for haplotypes
				bool hapsLimitReached = false;
				int hapsLimit = 100; //Less than 100 haplotypes needed
				while(hapsLimitReached == false) {

					vector<double> qBmean;
					double LqBgene;
					fullHapsGene = fullHapsGeneOriginal;
					aPH.optimiseqBmeanFromPre((int)fullHapsGene.size(), r, C, &qBmean, &LqBgene); 


					//Sort haplotypes by size of frequencies
					typedef pair<double, int> hapFreq; //Temporary container for frequency and haplotype index
					vector<hapFreq> hapFreqs;
					for(unsigned int i=0; i<fullHapsGene.size(); i++) {
	
						hapFreqs.push_back(hapFreq(qBmean[i],i));
					}
					sort(hapFreqs.begin(), hapFreqs.end()); //Sorts by frequency, small to large

					int numBelowThreshold = 0;
					for(unsigned int i=0; i<hapFreqs.size(); i++) {
			
						cout << "Index " << i << " freq: " << hapFreqs[i].first << "\n";
						if(hapFreqs[i].first < cutOff) { numBelowThreshold++; }
					}
					cout << "Num below threshold: " << numBelowThreshold << "\n";


					
					//Remove haplotypes if frequency less than cutoff
					for(int i = ((int) qBmean.size())-1; i >= 0; i--) { //Loop through backwards so we can remove haplotypes
	
						if(qBmean[i] <= cutOff) {
		
							qBmean.erase(qBmean.begin()+i);
							fullHapsGene.erase(fullHapsGene.begin()+i);
						}
					} 

					cout << "Number of haplotypes after filtering for gene " << g << " : " << fullHapsGene.size() << "\n";

					if((int) fullHapsGene.size() <= hapsLimit) {

						hapsLimitReached = true;
						cout << "Haplotype limit reached. Number of haplotypes: " << fullHapsGene.size() << ". Continuing with analysis using these haplotypes.\n";
					} else {

						hapsLimitReached = false;
						if(cutOff * 10 < 0.001) {  //Check that cutoff doesn't get unrealistically high
							cutOff = cutOff * 10;
							cout << "Required number of haplotypes not reached. Increasing filtering cutoff to " << cutOff << "\n";
						} else {
							cout << "Unable to reach the required number of haplotypes. Exiting.\n";
							exit(1);
						}
					}
				}
			}

			filteredHaplotypes.push_back(fullHapsGene); //Save for later usage
			exit(1); //Temporary
		}
	}



	/*
	 * Case where we filter haplotypes within the C++ code (in contrast to A) using all the possible haplotypes that could explain the dataset
	 * or B) importing pre-filtered haplotypes). 
	 *
	 * The below process is for the filter 2 scenario:
	 *
	 * 	- Generate all possible haplotypes
	 * 	- Infer before frequencies only (i.e. we don't infer a covariance matrix!)
	 *	- Remove haplotypes with frequencies <= 2*10^-10
	 *	- Update data with respect to these haplotypes (might remove entire sets of partial haplotypes and result in removal of monomorphic loci)
	 *	- If there are more than 100 haplotypes remaining, try inferring frequencies again and remove frequencies <2*10^-10
	 *	- Keep doing this until < 100 haplotypes 
	 *	- Pass these along as "imported haplotypes" further downstream
	 *
	 * NB: This is an early version of this method. It works, but code may need to be cleaned up slightly
	 * and outcomes verified before official usage. Filter 3 is the currently used method for filtering.
	 */
	cutOff = 2e-10;
	if(ap->getFilterHaplotypesMethod() == 2) { //Filter method 2 triggered

		for(int g=0; g<getNumGenes(); g++) {

			vector<Sequence> fullHapsGene;
			if(data->getGenePresent(g) == true) {
				DataPH* d = data->getGene(g);

				//Generate full haps from partial haps
				vector<Sequence> phs = d->getSequences(); //Get partial haplotypes
				AnalyserPH aPH;
				aPH.loadData(d);
				aPH.constructFullHaps(&phs, &fullHapsGene);
				//d->computeContribs(fullHapsGene); //Updates contribs wrt full haplotypes

				cout << "Number of haplotypes before filtering for gene " << g << " : " << fullHapsGene.size() << "\n";


				//Infer frequencies for haplotypes
				int hapsLimit = 100; //Less than 100 haplotypes needed
				for(int reductionCounter=1; reductionCounter < 999; reductionCounter++) { //Keep reduction set of haplotypes until hapsLimit reached

					bool reductionOccurred = false;
					vector<double> qBmean;
					double LqBgene;
		
					/*  Need to update the data with respect to fullHapsGene (as changes to
					 *  fullHapsGene may render some partial haplotype set invalid which in
					 *  turn may lead to monomorphic. Thus we need to update accordingly.					
					 *  Using filtered haplotypes
					 */
					d->setImportedFullHaps(fullHapsGene); //Need to load haplotypes into DataPH object


					//Need to do further processing here
					//Essentially, two things can happen: 
					//1) Some loci in the full haplotypes are monomorphic. This means that this locus needs 
					//deleting from all partial haps. Need to update physical pos as well, if present.
					//2) Some partial haps may no longer be represented by a full haplotype. In this case we
					//need to shorten the haplotype until it can be merged with a shorter haplotype.
					d->updateDataWithRespectToImportedFullHaps(); //Also computes contribs.

					//The above might move partial haplotypes from one set to a shorter subset.
					//In this process, some sophs become illegal. This is the case if NaTot = 0,
					//which is not defined for our inference process.
					//To this end, we rerun the "removeLowCountObservations" method. This takes care of that.
					//It also updates lociCovered and physicalPos.
					//Note that removing possible monomorphic SNPs should not be necessary, as the imported
					//full haplotypes should already be polymorphic.
					d->removeLowCountObservations();

					//Get possibly truncated full haplotypes back
					fullHapsGene = d->getImportedFullHaps();

					d->computeContribs(fullHapsGene); //Updates contribs wrt full haplotypes
					aPH.optimiseqBmeanFromPre((int)fullHapsGene.size(), r, C, &qBmean, &LqBgene); 

					cout << "Printing haps and freqs:\n";
					for(unsigned int i=0; i<fullHapsGene.size(); i++) {

						cout << qBmean[i] << "\t"; fullHapsGene[i].print();
					}

					//Remove haplotypes if frequency less than cutoff
					for(int i = ((int) qBmean.size())-1; i >= 0; i--) { //Loop through backwards so we can remove haplotypes
	
						if(qBmean[i] <= cutOff) {
		
							qBmean.erase(qBmean.begin()+i);
							fullHapsGene.erase(fullHapsGene.begin()+i);
							reductionOccurred = true;
						}
					} 

					cout << "Number of haplotypes after filtering for gene " << g << " for reduction round " << reductionCounter <<  " : " << fullHapsGene.size() << "\n";

					if((int) fullHapsGene.size() <= hapsLimit) {

						cout << "Haplotype limit reached. Number of haplotypes: " << fullHapsGene.size() << ". Continuing with analysis using these haplotypes.\n";
						break;
					} else {

						if(reductionOccurred == false) {

							cout << "No further filtering reductions possible. Were unable to reach the required number of haplotypes. Exiting.\n";
							exit(1);
						}
					}

				}
			}

			filteredHaplotypes.push_back(fullHapsGene); //Save for later usage
		}
	}

	if(ap->getFilterHaplotypesMethod() == 3) { //Filter method 3 triggered

		filteredHaplotypes = filterHaps3(r, C, filePathFolder);
	}


	
	//Obtain full haplotypes
	vector<vector<Sequence> > fullHaps;
	cout << "Number of genes: " << getNumGenes() << "\n";
	for(int g=0;g<getNumGenes();g++) {
		if(data->getGenePresent(g) == true) {
			DataPH* d = data->getGene(g);

			cout << "Printing data for gene " << g << ":\n";
			d->print();
	
			//If either imported full haps have beeen loaded or if haplotypes have been filtered above
			if(d->getImportFullHaps() == true || ap->getFilterHaplotypesMethod() > 0) { 

				vector<Sequence> importedFullHaps;
				if(d->getImportFullHaps() == true) {
					cout << "Imported full haps have been loaded. Using these for analysis. Haps:\n";
					importedFullHaps = d->getImportedFullHaps();
					for(unsigned int i=0; i<importedFullHaps.size(); i++) {
						importedFullHaps[i].print();
					}
				} else {
			
					//Using filtered haplotypes
					d->setImportedFullHaps(filteredHaplotypes[g]); //Need to load haplotypes into DataPH object
					importedFullHaps = filteredHaplotypes[g];

				}

				//Need to do further processing here
				//Essentially, two things can happen: 
				//1) Some loci in the full haplotypes are monomorphic. This means that this locus needs 
				//deleting from all partial haps. Need to update physical pos as well, if present.
				//2) Some partial haps may no longer be represented by a full haplotype. In this case we
				//need to shorten the haplotype until it can be merged with a shorter haplotype.

				if(importedFullHaps.size() > 0) { //Only update data if there actually are imported full haplotypes 
				
					d->updateDataWithRespectToImportedFullHaps(); //Also computes contribs.

					//The above might move partial haplotypes from one set to a shorter subset.
					//In this process, some sophs become illegal. This is the case if NaTot = 0,
					//which is not defined for our inference process.
					//To this end, we rerun the "removeLowCountObservations" method. This takes care of that.
					//It also updates lociCovered and physicalPos.
					//Note that removing possible monomorphic SNPs should not be necessary, as the imported
					//full haplotypes should already be polymorphic.
					d->removeLowCountObservations();

				} else { //Remove all data, as no full haplotypes to match

					cout << "No imported full haps, so removing all partial haplotype for gene " << g << " if present.\n";
					d->clear();
					d->clearPhysicalPos(); //If no data, can't have physical pos either
					data->setGenePresent(g,false);	
				}


				//Recall full haps from data (in case it has changed in the above processing) and add to container
				importedFullHaps = d->getImportedFullHaps();
				fullHaps.push_back(importedFullHaps);


				//Physical positions is implemented somewhat badly. In particular, each data entity has its own copy, but the upper
				//levels do not make use of the lower levels, e.g. if the lower levels are changed, the upper levels don't.
				//At some future restructuring of the code, this should be attended to. For the time being, here is a simple function
				//that copies data from lower structures up to higher ones.
				if(data->physicalPosPresent()==true) {
					data->updatePhysicalPos(g, d);
				}

				d->computeContribs(importedFullHaps);

			} else {

				//Generate full haps from partial haps
				vector<Sequence> phs = d->getSequences(); //Get partial haplotypes
				vector<Sequence> fhs; //Full haplotypes for current gene
            
				cout << "Printing partial haplotypes before constructing full haplotypes:\n";
				for(unsigned int j=0; j<phs.size(); j++) {
					phs[j].print();
				}
			
				AnalyserPH aPH;
				aPH.loadData(d);
				aPH.constructFullHaps(&phs, &fhs);
//Add to use true full haps			fhs = d->getFullHapsTrue(); //TEMPORARILY ADDED TO CHECK IMPACT ON USING TRUE HAPS
				cout << "Minimum set of full haplotypes: \n";
				cout << "fhs.size(): " << fhs.size() << "\n";
				for(unsigned int j=0; j<fhs.size();j++) {
					fhs[j].print();
					cout << j+1 << " of " << fhs.size() << "\n";
				}

				d->computeContribs(fhs); //Updates contribs wrt full haplotypes
				fullHaps.push_back(fhs); //Add to full haps for all genes
			}

		} else { //I.e. no data for current gene, so add empty variables

			vector<Sequence> fhs;
			fullHaps.push_back(fhs);
		}
	}

	//Print inferred full haplotypes
	for(unsigned int i=0; i<fullHaps.size(); i++) {

		ss.str(""); //Clear stringstream
		ss << filePathFolder << "InferredFullHaps_Gene_" << i << ".dat";
		filePath = ss.str();
		outputFile.open(filePath.c_str());

		for(unsigned int j=0; j<fullHaps[i].size(); j++) {

			for(int k=0; k<fullHaps[i][j].getLength(); k++) {

				outputFile << fullHaps[i][j].getBase(k);
			}
			outputFile << "\n";
		}
		outputFile.close();
			
	}

	//Print physical positions
	vector<vector<int> > physPos = data->getPhysicalPos();
	ss.str(""); //Clear stringstream
	ss << filePathFolder << "physPos.dat";
	filePath = ss.str();
	outputFile.open(filePath.c_str());
	for(unsigned int i=0; i<physPos.size(); i++) { //Loop over genes

		for(unsigned int j=0; j<physPos[i].size(); j++) {

			outputFile << physPos[i][j] << " ";
		}
		outputFile << "\n";
	}
	outputFile.close();	

	//Flag of 1 indicates to terminate analysis after inference of full haplotypes (without inferring frequencies)
	if(ap->getTerminateEarly() == 1) {

		cout << "Terminating after haplotype inference as requested.\n";
		return;
	}

	//Infer frequencies qB as a mean vector and covariance matrix. Infer from xB only.
	vector<vector<double> > qBmean; //Outer vector level: gene
	vector<gsl_matrix*> qBvar;
	vector<double> LqB; //Likelihoods for qB optimisation

	for(int g=0;g<getNumGenes();g++) {
		if(data->getGenePresent(g) == true) {
			DataPH* d = data->getGene(g);

			cout << "Inferring full haplotype frequencies for gene " << g << ":\n";
			AnalyserPH aPH;
			aPH.loadData(d);
            
			//Compute qB with mean and var from pre-transmission only
			vector<double> qBmeanGene;
			gsl_matrix* qBvarGene = gsl_matrix_alloc(fullHaps[g].size(), fullHaps[g].size()); //Allocate here, otherwise can't be passed back
			double LqBgene = 0;
			cout << "List of full haps for gene " << g << ":\n";
			for(unsigned int i=0; i<fullHaps[g].size(); i++) {

				fullHaps[g][i].print();
			}

			if(ap->getMeanOnly() == true) { //Infer qBmean only, not qBvar

				aPH.optimiseqBmeanFromPre((int)fullHaps[g].size(), r, C, &qBmeanGene, &LqBgene);
				cout << "Printing optimised qB mean for gene " << g << ":\n";
				printDoubleVector(qBmeanGene);
		
				//Even though we don't use a matrix, we define it as a 1x1 matrix with 0 entries.
				//Safer than just using an unassigned matrix. Easier to spot if errors due occur.
				gsl_matrix* qBvarGene = gsl_matrix_alloc(1,1);
				gsl_matrix_set(qBvarGene, 0,0, 0); //Define a 1x1 matrix with element 0 as an "empty" matrix - no other way to do this, as can't have it unassigned or 0 size
		

			} else {
	
				aPH.optimiseqBmeanAndVarFromPre((int)fullHaps[g].size(), r, C, &qBmeanGene, qBvarGene, &LqBgene); 
				//aPH.optimiseqBmeanAndVarSA((int)fhs.size(), r, C, &qBmeanGene, qBvarGene, &LqBgene); 
	

				cout << "Printing optimised qB mean for gene " << g << ":\n";
				printDoubleVector(qBmeanGene);
				cout << "Printing optimised qB variance for gene " << g << ":\n";
				printMatrixMathematica(qBvarGene);
			}
			
			qBmean.push_back(qBmeanGene);
			qBvar.push_back(qBvarGene);
			LqB.push_back(LqBgene);
	
		} else { //I.e. no data for current gene, so add empty variables

			vector<double> qBmeanGene;
			qBmean.push_back(qBmeanGene);
			gsl_matrix* qBvarGene = gsl_matrix_alloc(1,1);
			gsl_matrix_set(qBvarGene, 0,0, 0); //Define a 1x1 matrix with element 0 as an "empty" matrix - no other way to do this, as can't have it unassigned or 0 size
			qBvar.push_back(qBvarGene);
			LqB.push_back(0);
		}
	}


	//Print qBmean
	ss.str(""); //Clear stringstream
	ss << filePathFolder << "qBmeanPre.dat";
	filePath = ss.str();
	outputFile.open(filePath.c_str());
	for(unsigned int i=0; i<qBmean.size(); i++) { //Loop over genes

		for(unsigned int j=0; j<qBmean[i].size(); j++) {

			outputFile << qBmean[i][j] << " ";
		}
		outputFile << "\n";
	}
	outputFile.close();
	

	//Print qBvar
	ss.str(""); //Clear stringstream
	ss << filePathFolder << "qBvarPre.dat";
	filePath = ss.str();
	outputFile.open(filePath.c_str());
	for(unsigned int i=0; i<qBvar.size(); i++) { //Loop over genes
		
		if(data->getGenePresent(i) == true) {
			string qBvarGeneString = printMatrixMathematicaToString(qBvar[i]);
		
			outputFile << qBvarGeneString << "\n";
		
		} else {
			outputFile << "\n";
		}
	}
	outputFile.close();

	//Print LqB
	ss.str(""); //Clear stringstream
	ss << filePathFolder << "LqB.dat";
	filePath = ss.str();
	outputFile.open(filePath.c_str());
	for(unsigned int i=0; i<LqB.size(); i++) { //Loop over genes
		
		if(data->getGenePresent(i) == true) {
		
			outputFile << LqB[i] << "\n";
		
		} else {
			outputFile << "\n";
		}
	}
	outputFile.close();

	//Flag of 2 indicates to terminate analysis after inference of full haplotypes and associated frequencies
	if(ap->getTerminateEarly() == 2) {

		cout << "Terminating after haplotype frequency inference as requested.\n";
		return;
	}

	//Read in within-selection and compute within-host fitness if present.
	//Only needs to be computed once, as depends only on full haplotypes.
	//Note: within-host selection data is not read in until now. This is to
	//ensure that physical pos has been updated accordingly (they may have 
	//been altered earlier in this file).
	string whFolder = ap->getWithinHostSelectionFolder();
	if(!whFolder.empty()) {

		data->loadWithinHostSelection(whFolder);
		DataPHMG::WHSelMG whSelMG = data->getWHSelMG();

		if(whSelMG.getSelPresent() == true) {
		
			hapFitGMG = computeHapFitGrowth(fullHaps, data->getDeltaDays());
			cout << "Printing out hapFitG:\n";
			for(unsigned int g=0; g<hapFitGMG.size(); g++) {
				cout << "HapFitGgene " << g << ": "; printDoubleVector(hapFitGMG[g]);
				cout << "Full haps with associated with fitnesses:\n";
				for(unsigned int i=0; i<fullHaps[g].size(); i++) {
					cout << hapFitGMG[g][i] << "\t"; fullHaps[g][i].print();
				}
			}
		}
	}

	/*
	/ Analyse full model of selectiom
	*/
	cout << "Before analysis.\n";
	vector<Model> bestModels = analyseSelection(ap,fullHaps, r, qBmean, qBvar);



	for(unsigned int numFittedParams=0; numFittedParams<bestModels.size(); numFittedParams++) {
            
		cout << "numFittedParams: " << numFittedParams << "\n";

		//Most probably scenario for numFittedParameters
		Model currentModel = bestModels[numFittedParams]; 

		//Check that BIC value is not infinity. If infinity, the model is invalid and shouldn't be printed
		if(currentModel.getBICbest() < numeric_limits<double>::max()) {
	
			//PRETTY SURE THESE ARE ALL THE SAME AND CAN BE DELETED!
			//CHECK!
			//Write qB to file
			ss.str(""); //Clear stringstream
			ss << filePathFolder << "qB_" << numFittedParams << "_Fitted_Params.dat";
			filePath = ss.str();
			outputFile.open(filePath.c_str());
			for(int i=0; i<getNumGenes(); i++) {
		
				vector<double> qB = currentModel.getqBmeanBest(i);
				for(unsigned int j=0; j<qB.size(); j++) {
					outputFile << qB[j] << " ";
				}
				outputFile << "\n";
			}
			outputFile.close();

		    
			//Write selection result to file
			ss.str(""); //Clear stringstream
			ss << filePathFolder << "Selection_" << numFittedParams << "_Fitted_Params.dat";
			filePath = ss.str();
			outputFile.open(filePath.c_str());

		
		    
			//Print BIC
			outputFile.setf(ios_base::fixed);
			outputFile << currentModel.getBICbest() << "\t";
			outputFile.unsetf(ios_base::fixed);
		    
			//Print other stats
			outputFile << currentModel.getNtBest() << "\t";
		    
			//Print one gene at a time, all on one tab delimited line
			outputFile << getNumGenes() << "\t";
			for(int j=0; j<getNumGenes(); j++) {
			
				Model::modelSingleGene MSGcurrent = currentModel.getModelSingleGene(j);
				DiploidSequence seq = MSGcurrent.getModel();
		     
				int lengthOfCurrentGene = seq.getMajor().getLength();
				outputFile << lengthOfCurrentGene << "\t";

				//Print major and minor sequence of selection model
				outputFile << seq.getMajor().printToString() << "\t" << seq.getMinor().printToString() << "\t";
			
				vector<double> selCoefs = MSGcurrent.getSelCoefsBest();

				//Print selection coefficients (for major)
				for(unsigned int k=0; k<selCoefs.size();k++) {
					outputFile << selCoefs[k] << "\t";
					
				}
			
			
				//Print epi model
				outputFile << MSGcurrent.getEpiModelLength() << "\t";
				for(int k=0;k<MSGcurrent.getEpiModelLength(); k++) {
			    
					//Print positions
					vector<int> positions = MSGcurrent.getEpiModel()[k].getPositions();
					outputFile << positions.size() << "\t"; //Print number of positions
					for(unsigned int l=0; l<positions.size(); l++) {
						outputFile << positions[l] << "\t"; //Print actual positions
					}
			    
					//Print corresponding coefficient
					outputFile << MSGcurrent.getEpiCoefsBest(k) << "\t";
			
				}
		       
			}

			outputFile << "\n";
			outputFile.close();

		}

	}

	/*
 	/ Use single locus models of analysis
	*/
	if(ap->getAnalyseSL() == true) {
		int NtSLApprox = analyseLeonard("Approximate");
		ss.str(""); //Clear stringstream
		ss << filePathFolder << "Nt_SL_LeonardApproximate.dat";
		filePath = ss.str();
		outputFile.open(filePath.c_str());
		outputFile << NtSLApprox;
		outputFile.close();
	
		int NtSLExact = analyseLeonard("Exact");
		ss.str(""); //Clear stringstream
		ss << filePathFolder << "Nt_SL_LeonardExact.dat";
		filePath = ss.str();
		outputFile.open(filePath.c_str());
		outputFile << NtSLExact;
		outputFile.close();

	}
	
}

int AnalyserPHMG::getNumGenes() { return AnalyserPHMG::data->getNumGenes(); }


vector<Model> AnalyserPHMG::analyseSelection(AnaParam *ap, vector<vector<Sequence> > &fullHaps, gsl_rng *r, vector<vector<double> >& qBmean, vector<gsl_matrix*>& qBvar) {
    

	//Approach as of 09/06
	//Outline
	//Construct collapsed set (still relevant)
	//Construct infinite for loop with break at relevant time (convert to a more correct while loop at a later stage)
	//For each numParam in for loop, find next set of models based on previous best
	//If no new models found (almost impossible, break out)
	//Compute L for these and find new best (add to bestModels)
	//If no new best found, break out


	//First construct the "collapsed set", i.e. the diploidSequence describing the full haps
	vector<DiploidSequence> collapsedFullHaps;
	for(unsigned int i=0; i<fullHaps.size(); i++) {
		collapsedFullHaps.push_back(constructDiploidFromFullHap(fullHaps[i]));
	}
	
	vector<Model> bestModels;
    
	for(int i=0; i<=ap->getMaxNumFittedParams(); i++) { //Loop over number of parameters to be fitted
		Model empty; //Empty gets initiated with infinite BIC and no selection
		bestModels.push_back(empty);
        
		if(i==0) { //Neutral case
			for(unsigned int j=0; j<collapsedFullHaps.size(); j++) { //Loop over genes
				
				Model::modelSingleGene neutralMSG = Model::modelSingleGene();
				DiploidSequence emptyDS = DiploidSequence(collapsedFullHaps[j].getLength());
				neutralMSG.setModel(emptyDS);
				neutralMSG.setSelectionPresent(false); //Neutral gene, no selection present
				
				//Set qBmean and qBvar to previously optimised values
				neutralMSG.setqBmeanBest(qBmean[j]);
				neutralMSG.setqBmeanNew(qBmean[j]); 

				neutralMSG.setqBvarBest(qBvar[j]);
				neutralMSG.setqBvarNew(qBvar[j]);

				//No epistasis, so add straight away
				bestModels[i].addModelSingleGene(neutralMSG);
			}

			//Run analysis for neutral scenario
			analyseSelectionCombination(ap, fullHaps, r, bestModels[i]); //TODO: Is it fine to use r and not rPar here? I guess so, as not parallel region. Should it be parallelised?

			//Print
			cout << "Best neutral outcome:\n";
			cout << "BIC: \t" << bestModels[i].getBICbest() << "\n";
			bestModels[i].print();


			//Check if Nt is equal to maxNt, this could be indicative of a problem
			if(bestModels[i].getNtBest() == ap->getMaxNt()) {

				cout << "Inferred Nt for model " << i << " was " << bestModels[i].getNtBest() << " which is equal to the maximum inferrable Nt. Whilst this isn't universally true, generally this indicates a potential problem with the data, e.g. if inferred qBmean is << inferred qBvar. Try rerunning the code with the flag '-noVar true'.\n";
				
				ofstream outputFile;
				string filePathFolder = ap->getOutputFolder();
				stringstream ss;
				ss << filePathFolder << "maxNtError.dat";
				string filePath = ss.str();
				outputFile.open(filePath.c_str());
				outputFile << "Inferred Nt for model " << i << " was " << bestModels[i].getNtBest() << " which is equal to the maximum inferrable Nt. Whilst this isn't universally true, generally this indicates a potential problem with the data, e.g. if inferred qBmean is << inferred qBvar. Try rerunning the code with the flag '-noVar true'.";
				outputFile.close();
			}

	
		} else { //Selection case

			vector<Model> currentModels = generateModelsFromPreviousBest(collapsedFullHaps,bestModels[i-1]); //Vector of sel models with i number of parameters to be analysed
			if(currentModels.size() == 0) { //No new models found, should rarely happen, but could in theory
				cout << "Breaking out of analysis loop as no additional models were found for analysis.\n";
				break;
			}


			#pragma omp parallel default(none) shared(i, cout, ap, fullHaps, bestModels, gsl_rng_default, currentModels, r)
			{

				#pragma omp for schedule(runtime) //Main for loop to be parallised
				for(unsigned int j=0; j<currentModels.size(); j++) { //Loop over different models in this category
            
					//Set up the RNG (can't have a shared RNG - makes no sense in parallel).
					//We can no longer recreate results using seeds, as parallelisation is unpredictable
					gsl_rng * rPar;
					const gsl_rng_type * TPar;
					gsl_rng_env_setup();
					TPar = gsl_rng_default;
					rPar = gsl_rng_alloc (TPar);
					gsl_rng_set(rPar, time(NULL)); //Set seed as current time

					//Run analysis
					if(omp_get_num_threads()>1) {
						cout << "Analysing in parallel.\n";
						analyseSelectionCombination(ap, fullHaps, rPar, currentModels[j]);
					} else {
						cout << "Analysing in serial.\n";
						analyseSelectionCombination(ap, fullHaps, r, currentModels[j]);
					}
	

					//Check for issues, e.g. very large selection
					bool issueFound = false;
					for(int g=0; g<getNumGenes(); g++) {

						Model::modelSingleGene MSG = currentModels[j].getModelSingleGene(g);
						vector<double> selCoefs = MSG.getSelCoefsBest();
				
						//Check for issues in sel coefs
						for(unsigned int s=0; s<selCoefs.size(); s++) {

							///If any selection coefficients >+/-9, then flag as issue
							if(selCoefs[s] > 9 || selCoefs[s] < -9) {
								issueFound = true; break;
							}
						}
		
						//Check for issues in epi coefs
						for(int e=0; e<MSG.getEpiModelLength(); e++) {

							//If any epistasis coefficients >+/-9, then flag as issue
							if(MSG.getEpiCoefsBest(e) > 9 || MSG.getEpiCoefsBest(e) < -9) {
								issueFound = true; break;
							}
						}
					}
					if(issueFound == true) { //There is an issue with this model, so set BIC to infinity
						currentModels[j].setBICbest(numeric_limits<double>::max());
					
					}	


					//Print outcome
					#pragma omp critical (print1) //Ensure printing in correct order
					{
						cout << "Running (multi locus) analysis for " << i << " fitted parameters. Currently at combination " << j+1 << " of " << currentModels.size() << ". Thread used: " << omp_get_thread_num() << "\n";

						//Print current model
						currentModels[j].print();

						cout << "BIC: \t" << currentModels[j].getBICbest() << "\n";
					}
            
					#pragma omp critical (update_BIC) //Ensure no erros when writing to this object
					{
						if(currentModels[j].getBICbest() < bestModels[i].getBICbest()) {
							bestModels[i]=currentModels[j];
						}
					}


					//Check if Nt is equal to maxNt, this could be indicative of a problem
					if(bestModels[i].getNtBest() == ap->getMaxNt()) {

						cout << "Inferred Nt for model " << i << " was " << bestModels[i].getNtBest() << " which is equal to the maximum inferrable Nt. Whilst this isn't universally true, generally this indicates a potential problem with the data, e.g. if inferred qBmean is << inferred qBvar. Try rerunning the code with the flag '-noVar true'.\n";
				
						ofstream outputFile;
						string filePathFolder = ap->getOutputFolder();
						stringstream ss;
						ss << filePathFolder << "maxNtError.dat";
						string filePath = ss.str();
						outputFile.open(filePath.c_str());
						outputFile << "Inferred Nt for model " << i << " was " << bestModels[i].getNtBest() << " which is equal to the maximum inferrable Nt. Whilst this isn't universally true, generally this indicates a potential problem with the data, e.g. if inferred qBmean is << inferred qBvar. Try rerunning the code with the flag '-noVar true'.";
						outputFile.close();
					}

				}
			}
        
			cout << "-----------------------------------------------------\n";
			cout << "Best model for " << i << " parameters:\n";
			bestModels[i].print();
			cout << "BIC best: " << bestModels[i].getBICbest() << "\n";
			cout << "-----------------------------------------------------\n";


			//Compare with best model for previous number of fitted parameters
			//if(bestModels[i].getBICbest() > bestModels[i-1].getBICbest()) {
			if(bestModels[i].getBICbest() > bestModels[i-1].getBICbest()-ap->getBICpenalty()) {
			
				//Remove the most recently added model
				//bestModels.pop_back(); //Temporarily added
				
				cout << "Breaking out of analysis loop as no further improvements were found.\n";
				cout << "Best model found:\n";
				bestModels[i-1].print();
				cout << "BIC best: " << bestModels[i-1].getBICbest() << "\n";

				break; //Don't add any more models
			}
		}

		if(i==ap->getMaxNumFittedParams()) {

			cout << "Optimisation interrupted as maximum number of fitted parameters reached: " << ap->getMaxNumFittedParams() << ".\n";
		}
	}
	return bestModels;
}

int AnalyserPHMG::analyseLeonard(string method) {

	int maxNt = 500; //Get from ap later
	vector<double> logLNts(maxNt,0);

	for(int g=0; g<getNumGenes(); g++) {

		cout << "Gene " << g << "\n";
		DataPH* dPH = data->getGene(g);
		AnalyserPH aPH;
		aPH.loadData(dPH);
		if(method.compare("Approximate") == 0) {
			aPH.analyseLeonard(logLNts);
		
		} else if(method.compare("Exact") == 0) {
		
			aPH.analyseLeonardExact(logLNts);
		}

	}
	
	//Find optimal Nt
	int bestNt = -1;
	double bestLogL = -numeric_limits<double>::max();
	for(int i=1;i<=logLNts.size(); i++) {

		cout << i << " " << logLNts[i-1] << "\n";
		if(logLNts[i-1] > bestLogL) {

			bestNt = i;
			bestLogL = logLNts[i-1];
		}
	}

	cout << "Best Nt was " << bestNt << " with logL of " << bestLogL << "\n";
	
	return bestNt;
}

void AnalyserPHMG::analyseSelectionCombination(AnaParam *ap, vector<vector<Sequence> > &fullHaps, gsl_rng *r, Model &M) {

	M.initialiseAndRandomiseCoefs(r, fullHaps);
	if(ap->getSelectionCapPresent() == true) {
		M.setSelectionCap(ap->getSelectionCap());
	}

	//Parameters defining update loop
	int nUpdate = 100;
	double delta = 1;
	int maxNt = ap->getMaxNt(); //Upper limit
	int nRunsTotal = 0;
    
	int k = 0;
	bool discreteSolutionComputed = false;
	while(delta > 0.01 || discreteSolutionComputed == false) { //Original 0.01

		if(delta >0.01) {

			k++;
			if(k%nUpdate==0) {
				double acceptanceRate = M.updateAcceptanceRates();
				//cout << "Acceptance rate: " << acceptanceRate << "\tDelta: " << delta << "\tnUpdate: " << nUpdate << "\n";
				double tempDelta=delta*(0.90+acceptanceRate);
				if (tempDelta < 5) { //If delta > 5, computations get out of hand
					delta=tempDelta;//delta*(0.95+acceptanceRate);
				}
            
				//Update value for nUpdate, i.e. the number of rounds before the acceptance value is computed
				//For high acceptance values we want to spend more time with the current delta
				//For low acceptance value, we wish to progress towards smaller deltas faster
				if(acceptanceRate > 0.1) {
					if(nUpdate < 500) { //Limit
						nUpdate += 10;  //Increase nUpdate
					}
				} else {
					if(nUpdate >= 20) {
						nUpdate -= 10; //Decrease nUpdate only if resulting value is >= 10
					}
				}

				//Reset counter k
				k = 1;
		
				cout << "Best BIC: " << M.getBICbest() << "\tdelta: " << delta << "\tnUpdate: " << nUpdate <<  "\n";
				cout << "NtBest: " << M.getNtBest() << "\n";
				//cout << "Current virtual memory usage: " << getValue() << "\n";
			}
       		
	 
			//Make incremental changes
			M.update(r, delta,maxNt);
	
		} else { //Final round, compute discrete solution

			discreteSolutionComputed = true;
			ap->setUseIntegralApproach(true);
			cout << "Best BIC before computing discrete solution: " << M.getBICbest() << "\tdelta: " << delta << "\tnUpdate: " << nUpdate <<  "\n";
		}		



		double Ltot=0;
		int Ntot = 0;
		for(int l=0;l<M.getNumGenes();l++) { //Loop over genes

			if(data->getGenePresent(l) == true) { //Data avaiable for current gene
				DataPH* d = data->getGene(l);
				AnalyserPH aPH;
				aPH.loadData(d);
				//cout << "Gene " << l << "\n";
            
				//Compute likelihood
				double Lfreq = 0;
				
				//Case with selection for transmission
				if(M.getSelectionPresent(l) == true) {
					
					vector<double> hapFitT = M.computeHapFitTransmission(l, fullHaps[l]); //Compute hapfit for gene l only

					//There are 2 subcases: with and without within-host transmission
					if(data->getWHSelPresent(l) == true) {

						Lfreq = aPH.computeLSelTSelGVar(ap, M.getNtNew(), M.getqBmeanNew(l), M.getqBvarNew(l), hapFitT, hapFitGMG[l]);

					} else{

						Lfreq = aPH.computeLSelTVar(ap, M.getNtNew(), M.getqBmeanNew(l), M.getqBvarNew(l), hapFitT);
					}

				//Case without selection for transmission (neutral)
				} else { 

					//There are 2 subcases: with and without within-host transmission
					if(data->getWHSelPresent(l) == true) {
						
//						M.setNtNew(10);
						Lfreq = aPH.computeLNeutralSelGVar(ap,M.getNtNew(), M.getqBmeanNew(l), M.getqBvarNew(l), hapFitGMG[l]);
//						cout << "gene " << l << " selG likelihood: " << Lfreq << "\n";
//						exit(1);
//						if(M.getNtNew()>110) { exit(1); }
					
					} else {
	
//						M.setNtNew(10);
						Lfreq = aPH.computeLNeutralVar(ap,M.getNtNew(), M.getqBmeanNew(l), M.getqBvarNew(l));
//						cout << "gene " << l << " likelihood: " << Lfreq << "\n";
					}
//					if(M.getNtNew()>150 || M.getNtNew() < 10) { exit(1); }
//					cout << "Nt: " << M.getNtNew() << "\n";
				}
           			
//				cout << "Lfreq for gene " << l << ": " << Lfreq << "\n"; 
				Ltot += Lfreq;
				Ntot += d->getNtot();

			}
		}
//		cout << "Nt: " << M.getNtNew() << "\n";
//		cout << "Ltot: " << Ltot << "\n";

		//Check if Ltot is nan. From experience, nan errors normally results from rounding errors. The obvious example of rounding errors are cases of 
		//strong selection in which haplotype frequencies are pushed to such extents that the successive application of matrix multiplication etc. leads
		//to rounding errors. However, it may also happen in other cases, e.g. if the dimension of the full haplotypes are very large (more computations
		//needed, so rounding errors stack up) or if we (by chance) have filtered and optimised a set of full haplotypes badly. Overall, the rounding
		//errors leads to the resulting convariance matrix being slightly wrong. This in turn can lead to e.g. negative determinant, which makes the
		//log likelihood (in which we take log(det)) nan. In most cases the covariance matrix is only slightly off, i.e. by modifying e.g. the 10th
		//significant digit in one of the entries might change the determinant from very negative to very positive. In other words, the method is correct,
		//there are just some situations in which rounding errors stack up. In these cases, we just have to accept that we can't get a useful likelihood.
		//As a result, we terminate optimisation of the current selection model and return a BIC of infinity.
		if(::isnan(Ltot) != 0) { //Ltot is nan

			cout << "A likelihood value of nan was encountered. This is likely due to rounding errors stacking up. As therefore return the current selection model with infinity BIC and continue optimising subsequent models.\n";
			M.setBICbest(numeric_limits<double>::max());
			ap->setUseIntegralApproach(false); //Shouldn't be necessary, but reset nonetheless

			return; //Exit optimising this selection model

		}

		if(discreteSolutionComputed == false) { //Update BIC as normal

			//Updates the best BIC value of new BIC is better. If not, it reverts changes. Also updates number of acceptances
			M.updateBICbest(-2*Ltot+M.getNumParamsToBeFitted()*log(Ntot)); //Compute BIC from likelihood and number of parameters
//			cout << "Best BIC: " << M.getBICbest() << "\tdelta: " << delta <<  "\n";
//			M.print();

		} else { //When computing the discrete solution (final step before exiting), update only the BIC value, not the other parameters

			M.setBICbest(-2*Ltot+M.getNumParamsToBeFitted()*log(Ntot));
			cout << "Best BIC after computing discrete solution: " << M.getBICbest() << "\tdelta: " << delta << "\tnUpdate: " << nUpdate <<  "\n";
			break;
		}

		nRunsTotal++;
	}

	cout << "Best BIC at end of run: " << M.getBICbest() << "\n";
	ap->setUseIntegralApproach(false); //Reset
	cout << "Total number of runs: " << nRunsTotal << "\n";
}

//Assumes selGmacVec is in 12 hour units, such that 24 hour units = 2*12 hour units
//DeltaDays describe the difference between before and after, e.g. if before = day 3 and after = day 5, then deltaDays = 2.
//This effectively assumes that system can be modelled by an effective growth size followed by deltaDays worth of WH selection.
//This is a slight approximation. In theory it should be Ng1, whSel1, Ng2, whSel2 with whSel1=whSel2, but Ng2 >> Ng1.
vector<vector<double> > AnalyserPHMG::computeHapFitGrowth(vector<vector<Sequence> > &fullHaps, int deltaDays) {

	cout << "In compute hapFitGrowth\n";

	vector<vector<double> > hapFitG;
	DataPHMG::WHSelMG withinHostSelMG = data->getWHSelMG();
	
	for(unsigned int i=0; i<fullHaps.size(); i++) { //Loop over genes
		
		cout << "In computeHapFitGrowth gene " << i << "\n";
		vector<double> hapFitGgene(fullHaps[i].size(), 0);
		DataPHMG::WHSel withinHostSel = withinHostSelMG.getWHSel(i);
	
		if(withinHostSel.getSelPresent() == true) { //Only need to check for within host haplotype fitness modifications if wh sel is present

			Sequence selG = withinHostSel.getSelG();
			vector<double> selGmag = withinHostSel.getSelGmag();
			vector<Model::epistasis> epiG = withinHostSel.getEpiG();
			vector<double> epiGmag = withinHostSel.getEpiGmag();		

//			cout << "Printing selG for gene: "; selG.print();
//			cout << "Printing selGmag for gene: "; printDoubleVector(selGmag);
//			cout << "Printing haplotypes for gene: ";
//			for(unsigned int j=0; j<fullHaps[i].size(); j++) {
//
//				fullHaps[i][j].print();
//			}

		
			//Check for selection
			for(unsigned int j=0; j<fullHaps[i].size(); j++) { //Loop over haplotypes

//				cout << "Haplotype " << j << "\n";
				for(int k=0; k<fullHaps[i][j].getLength(); k++) { //Loop over each locus in haplotype

//					cout << "Locus " << k << "\n";
//					cout << "selGvec[i].getLength(): " << selGvec[i].getLength() << "\n";
					if(selG.getBase(k) != '-') { //Locus under selection

						if(fullHaps[i][j].getBase(k) == selG.getBase(k)) { //Haplotype under selection

							hapFitGgene[j] += deltaDays*2*selGmag[k]; //Multiply by 2 to get 24 hour units

							//Ensure that no haplotype fitness are larger than +/- 10
							if(hapFitGgene[j] > 10) {
								hapFitGgene[j] = 10;
							} else if(hapFitGgene[j] < -10) {
								hapFitGgene[j] = -10;
							}
						}
					}
				}
			}		

			//Check for epistasis
			if(epiG.size() > 0) { //Possibility of epistasis
	
				for(unsigned int j=0; j<fullHaps[i].size(); j++) { //Loop over haplotypes

					for(unsigned int k=0; k<epiG.size(); k++) { //Loop over epistatic models
	
						vector<int> positions = epiG[k].getPositions();
						bool epistasisPresent = true; //Assume it is true, then check if the haplotype fits epistatic model
					
						for(unsigned int l=0; l<positions.size(); l++) { //Loop over positions in epi model

							if(fullHaps[i][j].getBase(positions[l]) != selG.getBase(positions[l])) {

								epistasisPresent = false; break;
							}
						}

						if(epistasisPresent == true) {
	
							hapFitGgene[j] += deltaDays*2*epiGmag[k]; //Multiply by 2 to get 24 hour units

							//Ensure that no haplotype fitness are larger than +/- 10
							if(hapFitGgene[j] > 10) {
								hapFitGgene[j] = 10;
							} else if(hapFitGgene[j] < -10) {
								hapFitGgene[j] = -10;
							}
						}
					}
				}
			}

		}

		hapFitG.push_back(hapFitGgene);
	} 	

	return hapFitG;
}



/*
 * Case where we filter haplotypes within the C++ code (in contrast to A) using all the possible haplotypes that could explain the dataset
 * or B) importing pre-filtered haplotypes). 
 *
 * The below process is for the filter 3 scenario:
 *
 * 	- Generate all possible haplotypes
 * 	- Infer before frequencies only (i.e. we don't infer a covariance matrix!)
 * 	- Sort haplotypes by frequency magnitude
 * 	- Remove one haplotype at a  time, starting with the ones with lowest frequency
 * 	- Check that the remaining full haplotypes can still represent all the partial haplotype data
 * 	- If this is not the case, don't remove the most recently removed haplotype. Continue with the next haplotype.
 * 	- Keep removing until the required number of haplotypes have been reached or no further haplotypes can be removed (fail scenario)
 */

//FilePathFolder argument is optional. If not supplied, required frequency haps limit will not be printed.
//This is useful if used to generate haplotypes only.
//Assumes dataPHMG has been loaded in to this AnalyserPHMG object
//Repindex is optional
vector<vector<Sequence> > AnalyserPHMG::filterHaps3(gsl_rng* r, double C, string filePathFolder /* = ""*/, int repIndex /* = -1 */) {

//	double cutOff = 2e-10; //Minimum cut-off that should be enforced at all times
	vector<vector<Sequence> > fullHaps;
	for(int g=0; g<data->getNumGenes(); g++) {

//		vector<Sequence> fullHapsGene;
//		if(data->getGenePresent(g) == true) {
//			DataPH* dPH = data->getGene(g);
////
//			//Generate full haps from partial haps
//			vector<Sequence> phs = dPH->getSequences(); //Get partial haplotypes
////			AnalyserPH aPH;
//			aPH.loadData(dPH);
//			aPH.constructFullHaps(&phs, &fullHapsGene);
//			dPH->computeContribs(fullHapsGene); //Updates contribs wrt full haplotypes
//
//			//Infer frequencies
//			vector<double> qBmean;
//			double LqBgene;
//			aPH.optimiseqBmeanFromPre((int)fullHapsGene.size(), r, C, &qBmean, &LqBgene); 
//
//			cout << "Number of haplotypes before filtering for gene " << g << " : " << fullHapsGene.size() << "\n";
//
//			//Sort haplotypes by size of frequencies
//			typedef pair<double, int> hapFreq; //Temporary container for frequency and haplotype index
//			vector<hapFreq> hapFreqs;
//			for(unsigned int i=0; i<fullHapsGene.size(); i++) {
//
//				hapFreqs.push_back(hapFreq(qBmean[i],i));
//			}
//			sort(hapFreqs.begin(), hapFreqs.end()); //Sorts by frequency, small to large
//
//
//			vector<int> toBeRemoved;
//			int hapsLimit = 100; //Less than 100 haplotypes needed
//			bool hapsLimitReached = false;
//			int kept = 0;
////			double currentCutOff = cutOff; //Real currentCutOff may be lower in reality, but unimportant to us
//			vector<Sequence> fullHapsGeneOld = fullHapsGene;
//			for(unsigned int i=0; i<hapFreqs.size(); i++) {
//		
//				cout << "Index " << i << " freq: " << hapFreqs[i].first << "\n";
//				cout << "currentCutOff: " << currentCutOff << "\n";
//
////				//Check that sufficiently many haplotypes have been removed such that #intial - #removed < limit
//				if((int) (fullHapsGene.size() - toBeRemoved.size()) <= hapsLimit && hapsLimitReached == false) {
//
//					cout << "Haps limit was reached:\n";
//					cout << "Num removed: " << toBeRemoved.size() << "\n";
//					cout << "Num kept: " << kept << "\n";
//					cout << "Might continue removing haplotypes if cutOff limit hasn't been reached.\n";
//
//					hapsLimitReached = true;
//	
//					//Print the required frequency for reaching the haps limit
//					if(filePathFolder.compare("") == false) { //Only print if an output folder is supplied
//						stringstream ss;
//						ss << filePathFolder << "requiredFrequencyHapsLimit_" << hapsLimit << "_gene" << g << ".dat";
//						string filePath = ss.str();
//						ofstream outputFile;
//						outputFile.open(filePath.c_str());
//						outputFile << currentCutOff;
//						outputFile.close();
//					}
//
//				}
//
//				//If haps limit has been reached, also check that we have removed haplotypes up to the cut-off limit, 
//				//e.g. we may have met the hapsLimit but still have frequencies of < cutOff.
//				//In this case we want to keep removing haplotypes until we have removed all haplotypes below the cut-off.
//				currentCutOff = hapFreqs[i].first; //Get the current frequency cutoff for haplotype i
//				if(hapsLimitReached == true && currentCutOff > cutOff) {
//
//					//Save new haplotypes
//					fullHapsGene = fullHapsGeneOld;
//					break;
//				}
//
//
//				/*
//				 * If limits haven't been reached, continue/start removing haplotypes
//				 */ 
//
//				//Add the next haplotype to new list of removed haplotypes
//				vector<int> toBeRemovedNew = toBeRemoved;
//				toBeRemovedNew.push_back(hapFreqs[i].second);
//
//				//Sort in ascending order and remove from fullHapsGeneNew						
//				sort(toBeRemovedNew.begin(), toBeRemovedNew.end()); //Ascending order
//				vector<Sequence> fullHapsGeneNew = fullHapsGene;
//				for(int j=((int) toBeRemovedNew.size()) -1; j>=0; j--) {
//
//					fullHapsGeneNew.erase(fullHapsGeneNew.begin() + toBeRemovedNew[j]);	
//				}
//
//				//Check if fullHapsGeneNew can de be decomposed to represent all partial haplotypes
//				if(dPH->decomposeAndCheck(fullHapsGeneNew) == true) {
//
//					//Keep change
//					toBeRemoved = toBeRemovedNew;
//					cout << "Haplotype " << i << " removed.\n";
//				} else {
//
//					//Don't keep change
//					cout << "Haplotype " << i << " kept.\n";
//					kept++;
//				}
//
//				//Next round new = old
//				fullHapsGeneOld = fullHapsGeneNew;
//
//			}
//
//			if(hapsLimitReached == true) {
//
//				cout << "Haps and cut-off limit were reached. Continuing.\n";		
//				
//				//I don't think this is necessary as recomputed downstream.
//				dPH->computeContribs(fullHapsGene); //Updates contribs wrt full haplotypes
//
//			} else {
//
//				cout << "Haps limit wasn't reached, so exiting.\n";
//				exit(1);
//			}
//		}


		//New approach: put all code into analyserPH
		DataPH* dPH = data->getGene(g);
		AnalyserPH aPH;
		aPH.loadData(dPH);
		vector<Sequence> fullHapsGene;
		if(repIndex != -1) {
			fullHapsGene = aPH.filterHaps3(r,C,filePathFolder,g,repIndex);
		} else {
			fullHapsGene = aPH.filterHaps3(r,C,filePathFolder,g);
		}

		fullHaps.push_back(fullHapsGene); //Save for later usage
	}

	return fullHaps;
}












