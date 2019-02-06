//
//  analyserPHMR.cpp
//  FluTransmissionProject
//
//  Created by Casper Lumby on 08/12/2015.
//  Copyright © 2015 Casper Lumby. All rights reserved.
//

#include "analyserPHMR.hpp"
#include <iostream>
#include "misc.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <limits>
#include <math.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "modelMR.hpp"

#ifdef _OPENMP
	#include <omp.h>
#else
	#define omp_get_thread_num() 0
#endif

using namespace std;

AnalyserPHMR::AnalyserPHMR() { } //Constructor does nothing
AnalyserPHMR::~AnalyserPHMR() { } //Destructor does nothing

void AnalyserPHMR::loadData(Data *d) {
    
    cout << "Loading data for full haplotypes. \n";
    data = dynamic_cast<DataPHMR *>(d);
}


void AnalyserPHMR::runAnalysis(AnaParam *ap) {
    
	cout << "Running analysis for multiple replicates.\n";
    
    
	//Set up RNG
	gsl_rng * rng;
	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	rng = gsl_rng_alloc (T);
	gsl_rng_set(rng, ap->getSeed());
    
	//Get output folder
	ofstream outputFile;
	string filePathFolder = ap->getOutputFolder();
	double C = ap->getC();

	vector<vector<vector<Sequence> > > filteredHaplotypes;
	if(ap->getFilterHaplotypesMethod() == 3) { //Filter method 3 triggered

		filteredHaplotypes = filterHaps3(rng, C, filePathFolder);
	}


	//Obtain full haplotypes
	vector<vector<vector<Sequence> > > fullHaps;
	for(int r=0;r<data->getNumReplicates(); r++) { //Loop over replicates

		DataPHMG *dPHMG = data->getReplicate(r);
		vector<vector<Sequence> > fullHapsPHMG;

		for(int g=0; g<dPHMG->getNumGenes(); g++) { //Loop over genes within replicate
			if(dPHMG->getGenePresent(g) == true) {

				DataPH* dPH = dPHMG->getGene(g);
				cout << "Printing data for rep " << r << " gene " << g << ":\n";
				dPH->print();
	
				//If either imported full haps have beeen loaded or if haplotypes have been filtered above
				if(dPH->getImportFullHaps() == true || ap->getFilterHaplotypesMethod() > 0) { 

					vector<Sequence> importedFullHaps;
					if(dPH->getImportFullHaps() == true) {
						cout << "Imported full haps have been loaded. Using these for analysis. Haps:\n";
						importedFullHaps = dPH->getImportedFullHaps();
						for(unsigned int i=0; i<importedFullHaps.size(); i++) {
							importedFullHaps[i].print();
						}
					} else {
			
						//Using filtered haplotypes
						dPH->setImportedFullHaps(filteredHaplotypes[r][g]); //Need to load haplotypes into DataPH object
						importedFullHaps = filteredHaplotypes[r][g];

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
						dPH->clearPhysicalPos(); //If not data, can't have physical pos either
						dPHMG->setGenePresent(g,false);	
					}

					//Recall full haps from data (in case it has changed in the above processing) and add to container
					importedFullHaps = dPH->getImportedFullHaps();
					fullHapsPHMG.push_back(importedFullHaps);
	
					//Physical positions is implemented somewhat badly. In particular, each data entity has its own copy, but the upper
					//levels do not make use of the lower levels, e.g. if the lower levels are changed, the upper levels don't.
					//At some future restructuring of the code, this should be attended to. For the time being, here is a simple function
					//that copies data from lower structures up to higher ones.
					if(dPHMG->physicalPosPresent()==true) {
						dPHMG->updatePhysicalPos(g, dPH);
					}


				} else {

					//Generate full haps from partial haps
					vector<Sequence> phs = dPH->getSequences();
					vector<Sequence> fhs; //Full haps for current gene
				
					cout << "Printing partial haplotypes before constructing full haplotypes:\n";
						for(unsigned int j=0; j<phs.size(); j++) {
						phs[j].print();
					}

					AnalyserPH aPH;
					aPH.loadData(dPH);
					aPH.constructFullHaps(&phs, &fhs);
					cout << "Constructed full haplotypes:\n";
					for(unsigned int j=0; j<fhs.size();j++) {
						fhs[j].print();
						cout << j+1 << " of " << fhs.size() << "\n";
					}		


					dPH->computeContribs(fhs);
					fullHapsPHMG.push_back(fhs);
				}
			
			} else { //I.e. no data for current gene, so add empty variables

				vector<Sequence> fhs;
				fullHapsPHMG.push_back(fhs);	
			}

		}

		fullHaps.push_back(fullHapsPHMG);

		//Physical positions are implemented somewhat badly. See description above.
		if(data->physicalPosPresent()==true) {
			data->updatePhysicalPos(r,dPHMG);
		}
	}


	//Print inferred full haplotypes
	stringstream ss;
	string filePath;
	for(unsigned int r=0; r<fullHaps.size(); r++) {
	
		for(unsigned int g=0; g<fullHaps[r].size(); g++) {
	
			ss.str(""); //Clear stringstream
			ss << filePathFolder << "InferredFullHaps_Rep" << r << "_Gene" << g <<".dat";
			filePath = ss.str();
			outputFile.open(filePath.c_str());
			
			for(unsigned int i=0; i<fullHaps[r][g].size(); i++) {
	
				for(int j=0; j<fullHaps[r][g][i].getLength(); j++) {

					outputFile << fullHaps[r][g][i].getBase(j);
				}
				outputFile << "\n";
			}
			outputFile.close();
		}
	}


       	//Print physical positions
	vector<vector<vector<int> > > physPos = data->getPhysicalPos();
	for(unsigned int r=0; r<physPos.size(); r++) {
		ss.str(""); //Clear stringstream
		ss << filePathFolder << "physPos_rep" << r << ".dat";
		filePath = ss.str();
		outputFile.open(filePath.c_str());

		for(unsigned int g=0; g<physPos[r].size(); g++) { //Loop over genes

			for(unsigned int i=0; i<physPos[r][g].size(); i++) {
			
				outputFile << physPos[r][g][i] << " ";
			}
			outputFile << "\n";
		}
		outputFile.close();
	}




	//Infer frequencies qB as a mean and covariance matrix. Infer from xB only.
	vector<vector<vector<double> > > qBmean;
	vector<vector<gsl_matrix*> > qBvar;
	vector<vector<double> > LqB; //Likelihoods for qB optimisation

	for(int r=0;r<data->getNumReplicates(); r++) { //Loop over replicates

		DataPHMG *dPHMG = data->getReplicate(r);
		vector<vector<double> > qBmeanPHMG;
		vector<gsl_matrix*> qBvarPHMG;
		vector<double> LqBPHMG;


		for(int g=0; g<dPHMG->getNumGenes(); g++) { //Loop over genes within replicate
			if(dPHMG->getGenePresent(g) == true) {

				DataPH* dPH = dPHMG->getGene(g);
				cout << "Inferring full haplotype frequencies for rep " << r << " gene " << g << ":\n";
				AnalyserPH aPH;
				aPH.loadData(dPH);

				cout << "full haplotypes: \n";
				for(unsigned int i=0; i<fullHaps[r][g].size(); i++) {
					fullHaps[r][g][i].print();
				}
	
				//Compute qB with mean and var from pre-transmission only
				vector<double> qBmeanPH;
				gsl_matrix* qBvarPH = gsl_matrix_alloc(fullHaps[r][g].size(), fullHaps[r][g].size()); //Allocate here, otherwise can't be passed back
				double LqBPH = 0;
                   	  	aPH.optimiseqBmeanAndVarFromPre((int)fullHaps[r][g].size(), rng, C, &qBmeanPH, qBvarPH, &LqBPH); 


				cout << "Printing optimised qB mean for replicate " << r << " gene " << g << ":\n";
				printDoubleVector(qBmeanPH);
				cout << "Printing optimised qB variance for replicate " << r << " gene " << g << ":\n";
				printMatrixMathematica(qBvarPH);

				qBmeanPHMG.push_back(qBmeanPH);
				qBvarPHMG.push_back(qBvarPH);
				LqBPHMG.push_back(LqBPH);

			} else { //I.e. no data for current gene, so add empty variables

				cout << "Not optisiming mean and var for rep " << r << " gene " << g << "\n";

				vector<double> qBmeanPH;
				qBmeanPHMG.push_back(qBmeanPH);
				gsl_matrix* qBvarPH = gsl_matrix_alloc(1,1);
				gsl_matrix_set(qBvarPH, 0,0, 0); //Define a 1x1 matrix with element 0 as an "empty" matrix - no other way to do this, as can't have it unassigned or 0 size
				qBvarPHMG.push_back(qBvarPH);
				LqBPHMG.push_back(0);

			}

		}

		qBmean.push_back(qBmeanPHMG);
		qBvar.push_back(qBvarPHMG);
		LqB.push_back(LqBPHMG);

	}

	//Print qBmean
	ss.str(""); //Clear stringstream
	for(unsigned int r=0; r<qBmean.size(); r++) { //Loop over reps
		ss.str(""); //Clear stringstream
		ss << filePathFolder << "qBmean_Rep" << r <<".dat";
		filePath = ss.str();
		outputFile.open(filePath.c_str());
	
		for(unsigned int g=0; g<qBmean[r].size(); g++) { //Loop over genes

			for(unsigned int h=0; h<qBmean[r][g].size(); h++) { //Loop over haplotypes

				outputFile << qBmean[r][g][h] << " ";
			}
			outputFile << "\n";
		}
		outputFile.close();
	}
	
	//Print qBvar
	for(unsigned int r=0; r<qBvar.size(); r++) { //Loop over reps
		ss.str(""); //Clear stringstream
		ss << filePathFolder << "qBvar_Rep" << r << ".dat";
		filePath = ss.str();
		outputFile.open(filePath.c_str());
		for(unsigned int g=0; g<qBvar[r].size(); g++) { //Loop over genes
		
			string qBvarGeneString = printMatrixMathematicaToString(qBvar[r][g]);
		
			outputFile << qBvarGeneString << "\n";
		}
		outputFile.close();
	}


	//Print LqB
	for(unsigned int r=0; r<LqB.size(); r++) {

		ss.str(""); //Clear stringstream
		ss << filePathFolder << "LqB_Rep" << r << ".dat";
		filePath = ss.str();
		outputFile.open(filePath.c_str());
		for(unsigned int g=0; g<LqB[r].size(); g++) {

			outputFile << LqB[r][g] << "\n";
		}
		outputFile.close();
	}


	//Read in within-selection and compute within-host fitness if present.
	//Only needs to be computed once, as depends only on full haplotypes.
	//Note: within-host selection data is not read in until now. This is to
	//ensure that physical pos has been updated accordingly (they may have 
	//been altered earlier in this file).
	string whFolder = ap->getWithinHostSelectionFolder();
	if(!whFolder.empty()) {

		data->loadWithinHostSelection(whFolder);

		for(int r=0; r<data->getNumReplicates(); r++) {

			DataPHMG* dataPHMG = data->getReplicate(r);
			DataPHMG::WHSelMG whSelMG = dataPHMG->getWHSelMG();
			vector<vector<double> > hapFitGMG;
		
			if(whSelMG.getSelPresent() == true) {

				AnalyserPHMG aPHMG;
				aPHMG.loadData(dataPHMG);
				hapFitGMG = aPHMG.computeHapFitGrowth(fullHaps[r], data->getDeltaDays(r));
				cout << "Pritning out hapFitG for replicate " << r << ":\n";
				for(unsigned int g=0; g<hapFitGMG.size(); g++) {

					cout << "HapFitGgene " << g << ": "; printDoubleVector(hapFitGMG[g]);
					cout << "Full haps with associated fitnesses:\n";
					for(unsigned int i=0; i<fullHaps[r][g].size(); i++) {
						cout << hapFitGMG[g][i] << "\t"; fullHaps[r][g][i].print();
					}
				}
			}

			hapFitGMR.push_back(hapFitGMG);
		}
	}


    /*
     / Analyse full model of selectiom
     */
	cout << "Before analysis.\n";
	vector<ModelMR> bestModels = analyseSelection(ap,fullHaps, rng, qBmean, qBvar);

        
	/*
 	/Print replicate information, e.g. positions, alleles, maps from replicate to all
	*/
	vector<vector<vector<int> > > map = bestModels[0].getMap(); //Get values from neutral model - will always be added even if no selection was found
	vector<vector<int> > allPos = bestModels[0].getAllPos();
	vector<DiploidSequence> collapsedFullHaps = bestModels[0].getCollapsedFullHaps();


	//Print all (physical) pos
	ss.str(""); //Clear stringstream
	ss << filePathFolder << "physPos.dat";
	filePath = ss.str();
	outputFile.open(filePath.c_str());
	for(unsigned int i=0; i<allPos.size(); i++) { //Loop over genes

		for(unsigned int j=0; j<allPos[i].size(); j++) {
			outputFile << allPos[i][j] << " ";
		}
		outputFile << "\n";
	}
	outputFile.close();	

	//Print alleles
	for(unsigned int i=0; i<collapsedFullHaps.size(); i++) { //Loop over genes

		ss.str(""); //Clear stringstream
		ss << filePathFolder << "Alleles_Gene_" << i << ".dat";
		filePath = ss.str();
		outputFile.open(filePath.c_str());
		outputFile << collapsedFullHaps[i].printToString();
		outputFile.close();	
	}

	//Print replicate maps
	for(unsigned int i=0; i<map.size(); i++) { //Loop over genes

		for(unsigned int j=0; j<map[i].size(); j++) { //Loop over replicates

			ss.str(""); //Clear stringstream
			ss << filePathFolder << "Map_Replicate_" << j << "_Gene_" << i << ".dat";
			filePath = ss.str();
			outputFile.open(filePath.c_str());

			for(unsigned int k=0; k<map[i][j].size(); k++) { //Loop over positions
				outputFile << map[i][j][k] << " ";
			}
			outputFile.close();
		}
	}


	/*
 	/Print results   
	*/
	for(unsigned int numFittedParams=0; numFittedParams<bestModels.size(); numFittedParams++) {
            
		cout << "numFittedParams: " << numFittedParams << "\n";

		//Most probably scenario for numFittedParameters
		ModelMR currentModel = bestModels[numFittedParams]; 

		//Check that BIC value is not infinity. If infinity, the model is invalid and shouldn't be printed
		if(currentModel.getBICbest() < numeric_limits<double>::max()) {

			/*
			* Write combined model to file
			*/
			ss.str(""); //Clear stringstream
			ss << filePathFolder << "Selection_" << numFittedParams << "_Fitted_Params.dat";
			filePath = ss.str();
			outputFile.open(filePath.c_str());
		    
			//Print BIC
			outputFile << currentModel.getBICbest() << "\t";
		    
			//Print other stats
			vector<int> currentNts = currentModel.getNtBest();
			outputFile << currentNts.size() << "\t"; //Number of bottlenecks to follow
			for(unsigned int i=0; i<currentNts.size(); i++) {
				outputFile << currentNts[i] << "\t";
			}
		    
			//Print one gene at a time, all on one tab delimited line
			outputFile << currentModel.getNumGenes() << "\t";
			for(int j=0; j<currentModel.getNumGenes(); j++) {
			
				Model::modelSingleGene MSGcurrent = currentModel.getModelSingleGene(j);
				DiploidSequence seq = MSGcurrent.getModel();
		     
				int lengthOfCurrentGene = seq.getMajor().getLength();
				outputFile << lengthOfCurrentGene << "\t";

				//Print major and minor sequence of model
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
 	* Write selection for individual replicates to files
 	*/
        
	for(unsigned int numFittedParams=0; numFittedParams<bestModels.size(); numFittedParams++) {
            
		cout << "numFittedParams: " << numFittedParams << "\n";

		ModelMR currentModelMR = bestModels[numFittedParams]; 

		//Check that BIC value is not infinity. If infinity, the model is invalid and shouldn't be printed
		if(currentModelMR.getBICbest() < numeric_limits<double>::max()) {

			for(int r=0;r<data->getNumReplicates(); r++) { //Loop over replicates

				//Most probable scenario for numFittedParameters
				Model currentModel = currentModelMR.getModelRep(r); //Get the model for replicate r
		    
				//Write result to file
				ss.str(""); //Clear stringstream
				ss << filePathFolder << "Rep_" << r << "_Selection_" << numFittedParams << "_Fitted_Params.dat";
				filePath = ss.str();
				outputFile.open(filePath.c_str());
		    
				//Print BIC
				outputFile << currentModel.getBICbest() << "\t";
		    
				//Print other stats
				outputFile << currentModel.getNtBest() << "\t";
		    
				//Print one gene at a time, all on one tab delimited line
				outputFile << currentModel.getNumGenes() << "\t";
				for(int j=0; j<currentModel.getNumGenes(); j++) {
			
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
	}
}


vector<ModelMR> AnalyserPHMR::analyseSelection(AnaParam *ap, vector<vector<vector<Sequence> > > &fullHaps, gsl_rng* rng, vector<vector<vector<double> > >& qBmean, vector<vector<gsl_matrix*> >& qBvar) {

	//Find all loci from the replicates and return them in order
	vector<vector<int> > allPos = data->findAllPos(); //Of form: <gene<pos> >

	vector<vector<vector<int> > > map = data->findMapFromAllPosToReps(allPos); //Mapping of form vec<gene,rep,pos> with pos e.g. {-1,0,1,-1,2,3,-1,4} if zeroth index in all pos not in rep, first index in all pos corresponds to zeroth index in rep, etc..
	vector<DiploidSequence> collapsedFullHaps = createCollapsedFullHapsReplicates(fullHaps,map); //Ignores frequencies


	//From here the process is similar to that of analyserPHMG.cpp:
	vector<ModelMR> bestModels;


	for(int i=0; i<=ap->getMaxNumFittedParams(); i++) {
		ModelMR emptyMMR; //Empty gets initiated with infinite BIC and no selection
		emptyMMR.setMap(map);
		emptyMMR.setAllPos(allPos);
		emptyMMR.setCollapsedFullHaps(collapsedFullHaps);

//		if(i==1) {

//			cout << "NtBest[0]: " << bestModels[0].getNtBest(0) << "\n";
//			printMatrixMathematica(bestModels[0].qBvar[0][0]);
		//	printMatrixMathematica(emptyMMR.qBvar[0][0]);
			
//		}
		
		bestModels.push_back(emptyMMR);

		if(i==0) { //Neutral case

			for(unsigned int j=0; j<map.size(); j++) { //Loop over genes

				Model::modelSingleGene neutralMSG;
				cout << "collapsedFullHaps[j]: "; collapsedFullHaps[j].print();
				DiploidSequence emptyDS = DiploidSequence(collapsedFullHaps[j].getLength());
				cout << "Printing emptyDS: "; emptyDS.print();
				neutralMSG.setModel(emptyDS);
				neutralMSG.setSelectionPresent(false); //Neutral gene, no selection present

				bestModels[i].addModelSingleGene(neutralMSG); //Method defined in super	
			}

			//Set qBmean and var
			bestModels[i].setqBmean(qBmean); //Contains qBmeans for all replicates
			bestModels[i].setqBvar(qBvar);

			//Run analysis for neutral scenario
			analyseSelectionCombination(ap,  fullHaps, rng, bestModels[i]);

			//Print
			cout << "Best neutral outcome:\n";
			cout << "BIC: \t" << bestModels[i].getBICbest() << "\n";
			bestModels[i].print();

		} else { //Selection case

			vector<ModelMR> currentModels = generateModelsFromPreviousBest(collapsedFullHaps,bestModels[i-1]); //Vector of sel models with i number of parameters to be analysed

			if(currentModels.size() == 0) { //No new models found, should rarely happen, but could in theory
				cout << "Breaking out of analysis loop as no additional models were found for analysis.\n";
				break;
			}

			#pragma omp parallel default(none) shared(i, cout, ap, fullHaps, bestModels, gsl_rng_default, currentModels, rng)
			{

				#pragma omp for schedule(runtime) //Main for loop to be parallised
				for(unsigned int j=0; j<currentModels.size(); j++) { //Loop over different models in this category

					//Set up the RNG (can't have a shared RNG - makes no sense in parallel).
					//We can no longer recreate results using seeds, as parallelisation is unpredictable
					gsl_rng * rngPar;
					const gsl_rng_type * TPar;
					gsl_rng_env_setup();
					TPar = gsl_rng_default;
					rngPar = gsl_rng_alloc (TPar);
					gsl_rng_set(rngPar, time(NULL)); //Set seed as current time


					//Run analysis
					if(omp_get_num_threads()>1) {
						cout << "Analysing in parallel.\n";
						analyseSelectionCombination(ap, fullHaps, rngPar, currentModels[j]);
					} else {
						cout << "Analysing in serial.\n";
						analyseSelectionCombination(ap, fullHaps, rng, currentModels[j]);
					}


					//Check for issues, e.g. very large selection
					//Most of this should be caught at a different point in the code. Might be unnecssary here.
					bool issueFound = false;

					for(int r=0; r<data->getNumReplicates(); r++) {
		
						Model Mrep = currentModels[j].getModelRep(r);
						for(int g=0; g<Mrep.getNumGenes(); g++) {

							Model::modelSingleGene MSG = Mrep.getModelSingleGene(g);
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
					
						} else if(::isnan(currentModels[j].getBICbest())!=0)  {
							cout << "Current model # " << j << " has integral BIC which is NaN. No issues were found. Might need to lower tolerance to achieve convergence. Exiting\n";
							exit(1);
						}
					}


					//Print outcome
					#pragma omp critical (print1) //Ensure printing in correct order
					{
						cout << "Running (multi replicate) analysis for " << i << " fitted parameters. Currently at step " << j+1 << " of " << currentModels.size() << ". Results for combination " << j+1 << " below. Thread used: " << omp_get_thread_num() << "\n";

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

				}
			}

			cout << "-----------------------------------------------------\n";
			cout << "Best model for " << i << " parameters:\n";
			bestModels[i].print();
		        cout << "BIC best: " << bestModels[i].getBICbest() << "\n";
			cout << "-----------------------------------------------------\n";

			//Compare with best selection model for previous number of fitted parameters
			if(bestModels[i].getBICbest() > bestModels[i-1].getBICbest() - ap->getBICpenalty()) {

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

vector<vector<int> > AnalyserPHMR::findAllPos() {

	vector<vector<int> > allPos;
	
	vector<vector<vector<int> > > physicalPos = data->getPhysicalPos();

	for(unsigned int i=0; i<physicalPos.size(); i++) { //Loop over replicates


		for(unsigned int j=0; j<physicalPos[i].size(); j++) { //Loop over genes

			if(i==0) { //First time around, add empty vector for gene j
				vector<int> empty;
				allPos.push_back(empty);
			}

			for(unsigned int k=0; k<physicalPos[i][j].size(); k++) { //Loop over positions
			
				bool posAddedAlready = false;
				for(unsigned int l=0; l<allPos[j].size(); l++) {
					if(allPos[j][l] == physicalPos[i][j][k]) { posAddedAlready = true; break; }
				}

				if(posAddedAlready == false) { //Add new position to gene
					allPos[j].push_back(physicalPos[i][j][k]);
				}
			}
		}
	}


	//Sort allPos
	for(unsigned int i=0; i<allPos.size(); i++) { //Loop over genes
		sort(allPos[i].begin(), allPos[i].end());
	}

	return allPos;
}

vector<vector<vector<int> > > AnalyserPHMR::createMapFromAllPosToReps(vector<vector<int> > &allPos) {


	vector<vector<vector<int> > > map;
	vector<vector<vector<int> > > physicalPos = data->getPhysicalPos();

	//Loop over genes
	for(unsigned int i=0; i<allPos.size(); i++) {

		vector<vector<int> > empty1;
		map.push_back(empty1);

		//Loop over replicates
		for(unsigned int j=0; j<physicalPos.size(); j++) {

			vector<int> empty2;
			map[i].push_back(empty2);

			//Loop over positions in allPos
			for(unsigned int k=0; k<allPos[i].size(); k++) {

				int indexInRep = -1; //Minus 1 == not found

				//Loop over positions in replicate
				for(unsigned int l=0; l<physicalPos[j][i].size(); l++) {

					if(physicalPos[j][i][l] == allPos[i][k]) { 
						indexInRep = l;  break;
					}
				}
				
				map[i][j].push_back(indexInRep);
			}	
		}
	}

	return map;

}


vector<DiploidSequence> AnalyserPHMR::createCollapsedFullHapsReplicates(vector<vector<vector<Sequence> > > &fullHaps, vector<vector<vector<int> > > &map) {


//	cout << "In create collapsed full haps replicates.\n";
	vector<DiploidSequence> collapsedFullHaps;
	
	for(unsigned int i=0; i<map.size(); i++) { //Loop over genes

		DiploidSequence DSgene = DiploidSequence((int) map[i][0].size()); //Initialise empty diploid sequence
//		cout << "DSgene for gene " << i << "\n";
//		DSgene.print();
//		cout << "Map for gene " << i << " rep 0 \n";
//		printIntVector(map[i][0]);

		for(unsigned int j=0; j<map[i][0].size(); j++) { //Loop over positions

			bool posFound = false;
			char allele1 = '-';
			char  allele2 = '-';
			for(unsigned int k=0; k<map[i].size(); k++) { //Loop over replicates until replicate containing position j is found

				if(map[i][k][j] != -1) {
					
					posFound = true;
					allele1 = fullHaps[k][i][0].getBase(map[i][k][j]); //WLOG

					for(unsigned int l=0; l<fullHaps[k][i].size(); l++) { //Loop over full haplotypes in this rep

						char currentAllele = fullHaps[k][i][l].getBase(map[i][k][j]);
						if(currentAllele != allele1) {
							allele2 = currentAllele; break;
						}
					}
					break;
				}
			}

			if(posFound == true) {

				DSgene.setMajor(j,allele1); //Note that major and minor doesn't mean anything here
				DSgene.setMinor(j,allele2);

			} else { //THIS SHOULD NEVER HAPPEN 
					cout << "ERROR IN AnalyserPHMR::createCollapsedFullHapsReplicates: position not found amongst replicates.\n";
			}
		}

		collapsedFullHaps.push_back(DSgene);
	}

	return collapsedFullHaps;

}


vector<vector<int> > AnalyserPHMR::findSharedPos() {

	vector<vector<int> > sharedPos;
	
	vector<vector<vector<int> > > physicalPos = data->getPhysicalPos();
	
	vector<vector<int> > physicalPos0 = physicalPos[0]; //Pos for replicate 0
	
	for(unsigned int i=0; i<physicalPos0.size(); i++) { //Loop over genes

		vector<int> sharedPosGene;
		vector<int> physicalPos0i = physicalPos0[i]; //Pos for rep 0, gene i
		
		for(unsigned int j=0; j<physicalPos0i.size(); j++) { //Loop over positions

			bool shared = true;
			int currentPos = physicalPos0i[j];

			for(unsigned int k=1; k<physicalPos.size(); k++) { //Loop over remaining replicates


				vector<int> physicalPoski = physicalPos[k][i]; //Pos for rep k, gene i
				bool posFound = false;
				for(unsigned int l=0; l<physicalPoski.size(); l++) { //Loop over pos in rep k, gene i

					if(physicalPoski[l] == currentPos) {
						posFound = true;
						break;
					}					
				}

				if(posFound == false) {
					shared = false;
					break;
				}
			}

			if(shared == true) {
				sharedPosGene.push_back(currentPos);
			}

		}

		sharedPos.push_back(sharedPosGene);

	}
		
	return sharedPos;
}

//Create map from sharedPos to replicatePos
//Is this one deprecated? At least it is not in use any longer...
vector<vector<vector<int> > > AnalyserPHMR::createMap(vector<vector<int> > & sharedPos) {

	vector<vector<vector<int> > > map; //Levels: replicate, gene, map from shared to physical pos
	vector<vector<vector<int> > > physicalPos = data->getPhysicalPos();

	for(unsigned int i=0; i<physicalPos.size(); i++) { //Loop over replicates
	
		vector<vector<int> > mapReplicate;
		for(unsigned int j=0; j<sharedPos.size(); j++) { //Loop over genes within replicate

			vector<int> mapGene;
			for(unsigned int k=0; k<sharedPos[j].size(); k++) { //Loop over shared positions

				for(unsigned int l=0; l<physicalPos[i][j].size(); l++) { //Loop over physical pos
				
					if(sharedPos[j][k] == physicalPos[i][j][l]) {

						mapGene.push_back(l); //Pos l in physical pos corresponds to pos k in shared for current gene in current replicate
						break;
					}
				}
			}
			mapReplicate.push_back(mapGene);
		}
		map.push_back(mapReplicate);
	}

	return map;
}


void AnalyserPHMR::analyseSelectionCombination(AnaParam *ap, vector<vector<vector<Sequence> > > &fullHaps, gsl_rng *r, ModelMR &MMR) { 

	MMR.initialiseAndRandomiseCoefs(r, fullHaps);
	if(ap->getSelectionCapPresent() == true) {
		MMR.setSelectionCap(ap->getSelectionCap());
	}
	
	//Describe the update loop
	int nUpdate = 100;
	double delta = 1;
	int maxNt = ap->getMaxNt();
	int nRunsTotal = 0;

	int k = 0;
	bool discreteSolutionComputed = false;
	while(delta > 0.01 || discreteSolutionComputed == false) {

		if(delta > 0.01) {

			k++;
			if(k%nUpdate==0) {
				double acceptanceRate=MMR.updateAcceptanceRates();
				//cout << "Acceptance rate: " << acceptanceRate << "\tDelta: " << delta << "\tnUpdate: " << nUpdate << "\n";
				double tempDelta=delta*(0.90+acceptanceRate);
				if (tempDelta < 10) { //If delta > 10, computations get out of hand
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

				cout << "Best BIC: " << MMR.getBICbest() << "\tdelta: " << delta << "\tnUpdate: " << nUpdate <<  "\n";
				vector<int> currentNts = MMR.getNtBest();
				cout << "NtBest: "; printIntVector(currentNts);
//TBR?				cout << "NumParamsToBeFitted: " << MMR.getNumParamsToBeFitted() << "\n";
//				Model::modelSingleGene modelZero = MMR.getModelSingleGene(0);
///				vector<double> modelZeroSelCoefs = modelZero.getSelCoefsBest();
//				cout << "SelCoefs gene 0: "; printDoubleVector(modelZeroSelCoefs);
//				vector<double> modelZeroEpiCoefs = modelZero.getEpiCoefsBest();
//				if(modelZeroEpiCoefs.size() > 0) {
//					cout << "EpiCoefs gene 0: "; printDoubleVector(modelZeroEpiCoefs);
//				}
//			//	cout << "(Inferred) full haplotypes for gene 0:\n"; 
//			//	for(unsigned int i=0; i<fullHaps[0][0].size(); i++) {
//			//		fullHaps[0][0][i].print();
//			//	}
//				vector<double> hapFitTZero = modelZero.computeHapFitTransmission(fullHaps[0][0]);
//				cout << "HapFitTransmission for gene 0: "; printDoubleVector(hapFitTZero);
//				cout << "\n";
//				cout << "Current virtual memory usage: " << getValue() << "\n";

			}

			//Update values
			MMR.updateMR(r, delta, maxNt, ap->getUseSharedBottleneck());

		} else { //Final round, compute discrete solution

			discreteSolutionComputed = true;
			ap->setUseIntegralApproach(true);
			cout << "Best BIC before computing discrete solution: " << MMR.getBICbest() << "\tdelta: " << delta << "\tnUpdate: " << nUpdate << "\n";

		}


		double Ltot = 0;
		int Ntot = 0;
		for(unsigned int r=0; r<fullHaps.size(); r++) { //Loop over replicates

			Model Mrep = MMR.getModelRep(r); //Get the model for replicate i
			DataPHMG* dPHMG = data->getReplicate(r);

			for(int g=0; g<MMR.getNumGenes(); g++) { //Loop over genes

				if(dPHMG->getGenePresent(g) == true) { //Data available for current gene
					DataPH* dPH = dPHMG->getGene(g);
					AnalyserPH aPH;
					aPH.loadData(dPH);
			
					//Compute likelihood
					double Lfreq = 0;
					if(Mrep.getSelectionPresent(g) == true ) {

						vector<double> hapFitT = Mrep.computeHapFitTransmission(g, fullHaps[r][g]); //Compute hapfit for gene j only

						//There are two subcases: with and without within-host selection
//Old						if(ap->getWithinHostSelectionPresent(g) == true) {
						if(dPHMG->getWHSelPresent(g) == true) {

//Old							vector<double> hapFitG = ap->getHapFitG(r,g); //Probably needs fixing
							Lfreq = aPH.computeLSelTSelGVar(ap, Mrep.getNtNew(), Mrep.getqBmeanNew(g), Mrep.getqBvarNew(g), hapFitT, hapFitGMR[r][g]);
						} else{

							Lfreq = aPH.computeLSelTVar(ap, Mrep.getNtNew(), Mrep.getqBmeanNew(g), Mrep.getqBvarNew(g), hapFitT);
						}

					} else {
						
						//There are two subcases: with and without within-host selection
//						if(ap->getWithinHostSelectionPresent(g) == true) { 
						if(dPHMG->getWHSelPresent(g) == true) {

//Old							vector<double> hapFitG = ap->getHapFitG(r,g);
							Lfreq = aPH.computeLNeutralSelGVar(ap,Mrep.getNtNew(), Mrep.getqBmeanNew(g), Mrep.getqBvarNew(g), hapFitGMR[r][g]);

						} else {

							Lfreq = aPH.computeLNeutralVar(ap,Mrep.getNtNew(), Mrep.getqBmeanNew(g), Mrep.getqBvarNew(g));
						}

					}
            
					Ltot += Lfreq;
					Ntot += dPH->getNtot();
		
				}
			}

			Mrep.deallocateqB(); //Remove the matrices allocated in ModelMR::getModelRep()
		
		}


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
			MMR.setBICbest(numeric_limits<double>::max());
			ap->setUseIntegralApproach(false); //Shouldn't be necessary, but reset nonetheless

			return; //Exit optimising this selection model

		}

		if(discreteSolutionComputed == false) { //Update BIC as normal

			//Updates the best BIC value if new BIC is better. If not, it reverts changes. Also updates number of acceptances.
			MMR.updateBICbest(-2*Ltot+MMR.getNumParamsToBeFitted()*log(Ntot));

		} else { //When computing the discrete solution (final step before exiting), update only the BIC value, not the parameters

			MMR.setBICbest(-2*Ltot+MMR.getNumParamsToBeFitted()*log(Ntot));
			cout << "Best BIC after computing discrete solution: " << MMR.getBICbest() << "\tdelta: " << delta << "\tnUpdate: " << nUpdate << "\n";
			break;
		}

		nRunsTotal++;
		
	}

	cout << "Best BIC at end of run: " << MMR.getBICbest() << "\n";
	ap->setUseIntegralApproach(false); //Reset

	cout << "Total number of runs: " << nRunsTotal << "\n";
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
//Assumes dataPHMR has been loaded in to this AnalyserPHMR object
vector<vector<vector<Sequence> > > AnalyserPHMR::filterHaps3(gsl_rng* rng, double C, string filePathFolder /* = ""*/) {

	vector<vector<vector<Sequence> > > fullHaps;
	for(int r=0; r<data->getNumReplicates(); r++) {

		DataPHMG* dPHMG = data->getReplicate(r);
		AnalyserPHMG aPHMG;
		aPHMG.loadData(dPHMG);
		vector<vector<Sequence> > fullHapsRep = aPHMG.filterHaps3(rng,C,filePathFolder,r);

		fullHaps.push_back(fullHapsRep);
	}

	return fullHaps;

}

