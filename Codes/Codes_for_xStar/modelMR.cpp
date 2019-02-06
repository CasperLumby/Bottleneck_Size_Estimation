//
//  modelMR.cpp
//  FluTransmissionProject
//
//  Created by Casper Lumby on 27/01/2016.
//  Copyright © 2016 Casper Lumby. All rights reserved.
//

#include "modelMR.hpp"
#include <math.h>
#include <iostream>
#include "misc.h"
#include <limits>

using namespace std;

ModelMR::ModelMR() : Model() { }
ModelMR::~ModelMR() { }


void ModelMR::setMap(vector<vector<vector<int> > > &m) { map = m; }
void ModelMR::setAllPos(vector<vector<int> > &ap) { allPos = ap; }
void ModelMR::setCollapsedFullHaps(vector<DiploidSequence> &cfh) { collapsedFullHaps = cfh; }

vector<vector<vector<int> > > ModelMR::getMap() { return map; }
vector<vector<int> > ModelMR::getAllPos() { return allPos; }
vector<DiploidSequence> ModelMR::getCollapsedFullHaps() { return collapsedFullHaps; }

int ModelMR::getNtNew(int index) { return NtNewVec[index]; }
int ModelMR::getNtBest(int index) { return NtBestVec[index]; }
vector<int> ModelMR::getNtNew() { return NtNewVec; }
vector<int> ModelMR::getNtBest() { return NtBestVec; }

void ModelMR::print() {

	for(unsigned int i=0; i<modelAllGenes.size(); i++) {

		cout << "Gene " << i << ":\n";
		modelAllGenes[i].print();
	}
	cout << "Bottleneck for replicates:\n";
	for(unsigned int i=0; i<NtBestVec.size(); i++) {

		cout << NtBestVec[i] << " ";		
	}
	cout << "\n";
}

void ModelMR::printNew() {

	for(unsigned int i=0; i<modelAllGenes.size(); i++) {

		cout << "Gene " << i << ":\n";
		modelAllGenes[i].printNew();
	}
	cout << "Bottleneck for replicates:\n";
	for(unsigned int i=0; i<NtNewVec.size(); i++) {

		cout << NtNewVec[i] << " ";		
	}
	cout << "\n";
}

void ModelMR::initialiseAndRandomiseCoefs(gsl_rng *r, vector<vector<vector<Sequence> > > &haps) {

	//Number of bottleneck sizes to be fitted (one for now)
	//Consider using map[0].size() to get number of replicates - saves passing haps
	//Need to check that map[0] has the right number of replicates in all scenarios
	numNtToBeFitted = (int) haps.size(); //I.e. number of replicates

	//Initialise qB frequencies to the same value
	numqBtoBeFitted = 0;

//Temporarily commented out - assuming qB fixed by qBpre (with mean and variance)
//
// IF EVER UNCOMMENTED, IT NEEDS TO GET ADJUSTED FOR JOINT OPTIMISATION IN MULTIREPLICATE SCENARIO!!!
//
//	for(unsigned int i=0; i<modelAllGenes.size(); i++) { //Loop over genes
//		int dim = (int) haps[i].size();
//		vector<double> qBgene (dim,1); //vector of length # of haps in gene i with values 1
//		rescale(qBgene); //Rescale to make into frequency
//		modelAllGenes[i].setqBmeanNew(qBgene);
//		modelAllGenes[i].setqBmeanBest(qBgene);
//
//		//Add each qB element (i.e. haplotype) to numqBtoBeFitted
//		numqBtoBeFitted += dim; //Number of haplotypes for gene i
//
//		gsl_matrix* qBvarInitial = gsl_matrix_alloc(dim,dim);
//
//		//Initialise qBvar
//		for(int j=0; j<dim; j++) {
//
//			for(int k=0; k<dim; k++) {
//
//				//Matrix diagonals set to 0.1 - everything else set to zero
//				if(j==k) {
//					gsl_matrix_set(qBvarInitial, j, k, 0.000000001);
//				} else {
//					gsl_matrix_set(qBvarInitial, j, k, 0.0);
//				}
//			}
//			numqBtoBeFitted++;
//		}
//			
//		modelAllGenes[i].setqBvarNew(qBvarInitial); //Uses memcpy
//		modelAllGenes[i].setqBvarBest(qBvarInitial);
//		gsl_matrix_free(qBvarInitial);
//
//		//cout << "Test:\n";
//		//gsl_matrix* a = modelAllGenes[i].getqBvarNew();
//		//cout << "Dim a: " << a->size1 << "\n";
//		//gsl_matrix_set(a, 2,2, 500);
//		//gsl_matrix* b = modelAllGenes[i].getqBvarNew();
//		//cout << "B[2,2] = " << gsl_matrix_get(b,2,2) << "\n";	
//	}
//End of temporary comment


	//Initialise selection coefficients to random values in [0,5-0,5]
	selCoefsToBeFitted.clear(); //Clean up before counting again
	epiCoefsToBeFitted.clear();
	for(unsigned int i=0; i<modelAllGenes.size(); i++) { //Loop over genes
        
		Sequence currentGeneMajor = modelAllGenes[i].getModel().getMajor();
        
		vector<double> dummy(currentGeneMajor.getLength());
		modelAllGenes[i].setSelCoefsBest(dummy);
		modelAllGenes[i].setSelCoefsNew(dummy);
	
        
		//Randomise selection coefficients
		for(int j=0; j<currentGeneMajor.getLength();j++) { //Loop over positions
			if(currentGeneMajor.getBase(j) != '-') {

				//Actual code
				double selStrength = gsl_rng_uniform(r) -0.5;
				modelAllGenes[i].setSelCoefsBest(j, selStrength);
				modelAllGenes[i].setSelCoefsNew(j, selStrength);

				//Add coefficient to selCoefsToBeFitted
				vector<int> selCoefsPos = {i,j}; //{gene, pos}
				selCoefsToBeFitted.push_back(selCoefsPos);

				
			} else {
				modelAllGenes[i].setSelCoefsBest(j, 0);
				modelAllGenes[i].setSelCoefsNew(j, 0);
			}
		}
         
		//Randomise epistasis coefficients
		vector<epistasis> epiVec = modelAllGenes[i].getEpiModel();
		if(epiVec.size() > 0) {
			vector<double> dummyEpi(epiVec.size());
			modelAllGenes[i].setEpiCoefsBest(dummyEpi);
			modelAllGenes[i].setEpiCoefsNew(dummyEpi);
            
			for(unsigned int j=0; j<epiVec.size(); j++) {
				double epiStrength = gsl_rng_uniform(r) -0.5;
				modelAllGenes[i].setEpiCoefsBest(j, epiStrength);
				modelAllGenes[i].setEpiCoefsNew(j, epiStrength);


				//Add epi to epiCoefsToBeFitted
				vector<int> epiCoefsNum = {i, j}; //{gene,coefficient #}
				epiCoefsToBeFitted.push_back(epiCoefsNum);
			}
		}
	}
    
	//Initialise Nt
	NtNewVec.clear();
	NtBestVec.clear();
	for(int i=0; i<numNtToBeFitted; i++) {
		NtNewVec.push_back(100);
		NtBestVec.push_back(100);
	}

	//Initialse BIC
	BICbest=numeric_limits<double>::max();

	//Set initial acceptance rates proportional to the number of coefficients to be updated in each category
	//E.g. Nt only has one variable to be updated, but there may be multiple selection coefficients, so this
	//category should be updated more frequently
	double numNt = (double) numNtToBeFitted;
	double numqB = (double) numqBtoBeFitted;
	double numSelCoefs = (double) selCoefsToBeFitted.size();
	double numEpiCoefs = (double) epiCoefsToBeFitted.size();
	numParamsToBeFitted = numSelCoefs + numEpiCoefs; //Don't include numqB in this (for now at least)

	//Construct acceptanceRates vector which has dim for each set of parameters to be fitted (e.g. Nt only = 1, Nt and qB = 2 etc.)
	//Note that there is an assumed order, e.g. if there is selection there must also by qB and Nt or if there is epistasis, there
	//must also be selection, qB and Nt.
	acceptanceRates.clear();
	if(numNt>0) {
		acceptanceRates.push_back(numNt);
	}
	if(numqB > 0) {
		acceptanceRates.push_back(numqB);
	}
	if(numSelCoefs > 0) {
		acceptanceRates.push_back(numSelCoefs);
	}
	if(numEpiCoefs > 0) {
		acceptanceRates.push_back(numEpiCoefs);
	}
	rescaleMin(acceptanceRates,0.1); //Rescales and changes 0 to 0.01
//	cout << "Initial acceptanceRates: "; printDoubleVector(acceptanceRates);

	//Reset acceptance/tries counters
	accepts.clear();
	tries.clear();
	for(unsigned int i=0; i<acceptanceRates.size(); i++) { //Loop over 0=Nt, 1=qB, 2=sel, 3=epi
		accepts.push_back(0);
		tries.push_back(0);
	}
	acceptsOverall = 0;
	triesOverall = 0;
}

//Update method for multi replicates
void ModelMR::updateMR(gsl_rng *r, double delta, int maxNt, bool useSharedBottleneck) {

	currentlyBeingUpdated = -1; //-1=unassigned, 0=Nt, 1=qB, 2=sel, 3=epi

	double prob = gsl_rng_uniform(r);

	if(acceptanceRates.size() == 1) { //Single parameter to be updated

		if(numNtToBeFitted > 0) {
			currentlyBeingUpdated = 0; //Updating index 0 of acceptanceRates
			updateNt(r,delta,maxNt, useSharedBottleneck);		
		} else if(numqBtoBeFitted > 0) {
			currentlyBeingUpdated = 0;
			updateqB(r,delta);		
		} else if(selCoefsToBeFitted.size() > 0) {
			currentlyBeingUpdated = 0;
			updateSelCoefs(r,delta);		
		} else {
			cout << "Error in Model::update: acceptanceRates.size() == 1, but nothing to be updated. \n";
		}

		//Note: We assume here that we will never have to update epistasis without also updating selection
		//If that is the case, just add another clause
	}
	else if(acceptanceRates.size() == 2) { //Two types of parameters being updated (e.g. Nt and qB, Nt and sel, qB and sel or sel and epi)

		if(numNtToBeFitted>0 && numqBtoBeFitted>0) {

			if(prob < acceptanceRates[0]) {
				currentlyBeingUpdated = 0;
				updateNt(r, delta, maxNt, useSharedBottleneck);
			} else {
				currentlyBeingUpdated = 1;
				updateqB(r, delta);
			}
		} else if(numNtToBeFitted>0 && selCoefsToBeFitted.size()>0) {

			if(prob < acceptanceRates[0]) {
				currentlyBeingUpdated = 0;
				updateNt(r, delta, maxNt, useSharedBottleneck);
			} else {
				currentlyBeingUpdated = 1;
				updateSelCoefs(r, delta);
			}
		} else if(numqBtoBeFitted>0 && selCoefsToBeFitted.size()>0) {

			if(prob < acceptanceRates[0]) {
				currentlyBeingUpdated = 0;
				updateqB(r, delta);
			} else {
				currentlyBeingUpdated = 1;
				updateSelCoefs(r, delta);
			}
		} else if(selCoefsToBeFitted.size()>0 && epiCoefsToBeFitted.size()>0) {

			if(prob < acceptanceRates[0]) {
				currentlyBeingUpdated = 0;
				updateSelCoefs(r, delta);
			} else {
				currentlyBeingUpdated = 1;
				updateEpiCoefs(r, delta);
			}
		} else {
			cout << "Error in Model::update: acceptanceRates.size() == 2, but nothing to update. \n";
		}

	}
	else if(acceptanceRates.size() == 3) { //Three types of parameters being updated (e.g. Nt,qB and sel, Nt, sel and epi, qB,sel and epi)
		
		if(numNtToBeFitted>0 && numqBtoBeFitted>0 && selCoefsToBeFitted.size()>0) {
			if(prob < acceptanceRates[0]) {
				currentlyBeingUpdated = 0;
				updateNt(r, delta, maxNt, useSharedBottleneck);
			} else if (prob < (acceptanceRates[0]+acceptanceRates[1]))  {
				currentlyBeingUpdated = 1;
				updateqB(r, delta);
			} else {
				currentlyBeingUpdated = 2;
				updateSelCoefs(r, delta);
			}
		} else if(numNtToBeFitted>0 && selCoefsToBeFitted.size()>0 && epiCoefsToBeFitted.size()>0) {
			if(prob < acceptanceRates[0]) {
				currentlyBeingUpdated = 0;
				updateNt(r, delta, maxNt, useSharedBottleneck);
			} else if (prob < (acceptanceRates[0]+acceptanceRates[1]))  {
				currentlyBeingUpdated = 1;
				updateSelCoefs(r, delta);
			} else {
				currentlyBeingUpdated = 2;
				updateEpiCoefs(r, delta);
			}
		} else if(numqBtoBeFitted>0 && selCoefsToBeFitted.size()>0 && epiCoefsToBeFitted.size()>0) {
			if(prob < acceptanceRates[0]) {
				currentlyBeingUpdated = 0;
				updateqB(r, delta);
			} else if (prob < (acceptanceRates[0]+acceptanceRates[1]))  {
				currentlyBeingUpdated = 1;
				updateSelCoefs(r, delta);
			} else {
				currentlyBeingUpdated = 2;
				updateEpiCoefs(r, delta);
			}
		} else {
			cout << "Error in Model::update: acceptanceRates.siz() == 3, but nothing to update.\n";
		}
	}
	else if(acceptanceRates.size() == 4) { //All four: Nt, qB,sel and epi
	
		if(prob < acceptanceRates[0]) {
			currentlyBeingUpdated = 0;
			updateNt(r, delta, maxNt, useSharedBottleneck);
		} else if (prob < (acceptanceRates[0]+acceptanceRates[1]))  {
			currentlyBeingUpdated = 1;
			updateqB(r, delta);
		} else if (prob < (acceptanceRates[0]+acceptanceRates[1]+acceptanceRates[2])) {
			currentlyBeingUpdated = 2;
			updateSelCoefs(r, delta);
		} else {
			currentlyBeingUpdated = 3;
			updateEpiCoefs(r, delta);
		}
	}
	else { cout << "Error in Model::update(): acceptanceRates dimensions incorrect. Dimensions are: " << acceptanceRates.size() << ".\n"; }
}


//Override updateNt from model.cpp to reflect replicates
void ModelMR::updateNt(gsl_rng *r, double delta, int maxNt, bool useSharedBottleneck) {
    
	int repToBeUpdated = gsl_rng_uniform_int(r,(int)NtBestVec.size());
	int ceilDelta = (int) ceil(delta);
	int factor;
	if(NtBestVec[repToBeUpdated] < 100) {
		factor = 1;
	} else if(NtBestVec[repToBeUpdated] < 500) { //Between 100 and 500, jump in increments of maximum 5 (min 1)
		factor = gsl_rng_uniform_int(r,5)+1;
	} else if(NtBestVec[repToBeUpdated] < 1000) { //Between 500 and 1000, increment by maximum 10
		factor = gsl_rng_uniform_int(r,10)+1;
	} else {
		factor = gsl_rng_uniform_int(r,50)+1;
	}

	if(gsl_rng_uniform(r) > 0.5) {
		if(NtBestVec[repToBeUpdated] < (maxNt-ceilDelta*factor)) {NtNewVec[repToBeUpdated] += ceilDelta*factor; }
	} else {
		if(NtBestVec[repToBeUpdated] > (1+ceilDelta*factor)) { NtNewVec[repToBeUpdated] -= ceilDelta*factor; }
	}

	//Could be written more efficiently, but not a huge issues. Also reads easier.
	if(useSharedBottleneck == true) { //all bottlenecks must be the same, so update to same values as repToBeUpdated

		for(unsigned int i=0; i<NtBestVec.size(); i++) {
		
			NtNewVec[i] = NtNewVec[repToBeUpdated];
		}
	}
}


//Returns a Model object corresponding to the replicate identified by repIndex
//Output is essentially a projection from the MMR space to the M basis
Model ModelMR::getModelRep(int repIndex) {

//	cout << "In getModelRep.\n";
	Model M;
	M.setNtNew(NtNewVec[repIndex]);
	M.setNtBest(NtBestVec[repIndex]);
	vector<vector<int> > selCoefsToBeFittedRep;
//	cout << "M.getNtNew(): " << M.getNtNew() << "\n";

	//Loop over genes
	for(unsigned int i=0; i<map.size(); i++) {

		//Get selection model for single gene from modelMR, then change parameters accordingly
		Model::modelSingleGene MSG = modelAllGenes[i]; //MSG from MMR

		//Add qB for current replicate
		MSG.setqBmeanBest(qBmean[repIndex][i]);
		MSG.setqBmeanNew(qBmean[repIndex][i]);
		MSG.setqBvarBest(qBvar[repIndex][i]); //Uses memcpy
		MSG.setqBvarNew(qBvar[repIndex][i]); //Uses memcpy


		//Sort out selection coefficients
		DiploidSequence model = MSG.getModel();
//		cout << "model before correction:\n";
//		model.print();
		vector<double> selCoefsNew = MSG.getSelCoefsNew();
		vector<double> selCoefsBest = MSG.getSelCoefsBest();
//		cout << "selCoefsNew before correction:\n";
//		printDoubleVector(selCoefsNew);
		for(int j=(int)map[i][repIndex].size()-1; j>=0; j--) { //Go through positions in selModel in reverse
			if(map[i][repIndex][j] == -1) {
				model.removeBase(j);
				selCoefsNew.erase(selCoefsNew.begin()+j);
				selCoefsBest.erase(selCoefsBest.begin()+j);
			}
		}
//		cout << "model after correction:\n";
//		model.print();
//		cout << "selCoefsNew after correction:\n";
//		printDoubleVector(selCoefsNew);

		//Sort out selCoefsToBeFitted
		MSG.setSelectionPresent(false); //Assume false, then check below
		for(int j=0; j<model.getLength(); j++) { //Loop over pos in model
		
			if(model.getMajor().getBase(j) != '-') { //Selection present at locus j in gene i
				vector<int> currentSel = { i, j }; //Gene, position
				selCoefsToBeFittedRep.push_back(currentSel);
				MSG.setSelectionPresent(true);
//				cout << "Selection found. Printing model:\n";
//				model.print();
			}
		}
//		cout << "SelectionPresent in getModelRep for gene " << i << ": " << MSG.getSelectionPresent() << "\n";
		
		//Update MSG with the corrected values
		MSG.setModel(model);
		MSG.setSelCoefsNew(selCoefsNew);
		MSG.setSelCoefsBest(selCoefsBest);

		//Sort out epistatic effects - make sure the positions correspond to those in selModel
		vector<Model::epistasis> epiModel = MSG.getEpiModel();
		vector<double> epiCoefsNew = MSG.getEpiCoefsNew();
		vector<double> epiCoefsBest = MSG.getEpiCoefsBest();
		vector<int> episToBeRemoved;
		for(unsigned int j=0; j<epiModel.size(); j++) { //Loop over epistatic effects
		//	cout << "In epi correction loop.\n";
		//	cout << "EpiModel.size(): " << epiModel.size() << "\n";
			int numRemoved = 0;
			int epiCounter = 0;
			vector<int> epiPos = epiModel[j].getPositions();
		//	cout << "EpiPos: \n";
		//	printIntVector(epiPos);
		//	cout << "map[i][repIndex]:\n";
		//	printIntVector(map[i][repIndex]);
			for(unsigned int k=0; k<map[i][repIndex].size(); k++) { //Go through positions in map
			
		//		cout << "k: " << k << "\n";
				if(map[i][repIndex][k] == -1) {
					numRemoved++;
				}


				//If epistatic effect at this position
				if(epiPos[epiCounter] == (int)k) {
				
					//Check if conflict
					if(map[i][repIndex][k] == -1) {
						episToBeRemoved.push_back(j); break; //This site is NOT in this replicate model, so epistatic effect invalid
					}

					//If no conflict, update the epistatic position to reflect this replicate
					epiPos[epiCounter] -= numRemoved; //Adjust for removed loci
					epiCounter++;
				}
							
				if(epiCounter == (int) epiPos.size()) { break; } //No more positions to update	

			}
		}

		//Remove any invalid effects
		for(int j=(int)episToBeRemoved.size()-1; j>=0; j--) { //Go through values in reverse
		//	cout << "episToBeRemoved:\n";
		//	printIntVector(episToBeRemoved);
		//	cout << "j: " << j << "\n";
			epiModel.erase(epiModel.begin() + episToBeRemoved[j]);
			epiCoefsNew.erase(epiCoefsNew.begin() + episToBeRemoved[j]);
			epiCoefsBest.erase(epiCoefsBest.begin() + episToBeRemoved[j]);
		}

		//Update values in MSG with corrected epistatic effects
		MSG.setEpiModel(epiModel);
		MSG.setEpiCoefsNew(epiCoefsNew);
		MSG.setEpiCoefsBest(epiCoefsBest);
		
		//Add MSG to selection model
		M.addModelSingleGene(MSG);
	}	
	
	//Update selCoefsToBeFitted - is it essential to know gene and pos information to compute hapfits?
	//Seems silly to count selection twice: In the above and below when counting parameters.
	//Maybe we can simplify this
	M.setSelCoefsToBeFitted(selCoefsToBeFittedRep);

	M.countNumParamsToBeFitted();

	return M;
}

//Needs updated if ever start jointly updating qB, Nt and sel
void ModelMR::updateBICbest(double BICnew) {

	if(BICnew < BICbest) {

		NtBestVec = NtNewVec;
		for(unsigned int i=0; i<modelAllGenes.size(); i++) {

			modelAllGenes[i].setSelCoefsBestToNew();
			modelAllGenes[i].setEpiCoefsBestToNew();
		}

		BICbest = BICnew;

		//Update acceptance rates
		accepts[currentlyBeingUpdated] = accepts[currentlyBeingUpdated] + 1;
                acceptsOverall++;
	} else {

		NtNewVec = NtBestVec;
		for(unsigned int i=0; i<modelAllGenes.size(); i++) {

			modelAllGenes[i].setSelCoefsNewToBest();
			modelAllGenes[i].setEpiCoefsNewToBest();
		}
	}

	tries[currentlyBeingUpdated] = tries[currentlyBeingUpdated] + 1;
        triesOverall++;
}

void ModelMR::setqBmean(vector<vector<vector<double> > >& QBMEAN) { qBmean = QBMEAN; }
void ModelMR::setqBvar(vector<vector<gsl_matrix*> >& QBVAR) { qBvar = QBVAR; } //This does not use memcpy. Memcpy will be relevant if we optimise jointly over variables in the future (update method need change too in this case)
