//
//  model.cpp
//  FluTransmissionProject
//
//  Created by Casper Lumby on 27/01/2016.
//  Copyright © 2016 Casper Lumby. All rights reserved.
//

#include "model.hpp"
#include <math.h>
#include <iostream>
#include "misc.h"
#include <limits>
#include <gsl/gsl_linalg.h>

using namespace std;

Model::Model() : numParamsToBeFitted(-1), BICbest(numeric_limits<double>::max()), NtBest(-1), NtNew(-1), selectionCap(10) { }
Model::~Model() {

	//Deallocate any gsl_matrix objects
	for(unsigned int i=0; i<modelAllGenes.size(); i++) {

//		modelAllGenes[i].eallocateqB();
	}

}

void Model::deallocateqB() {

	for(unsigned int i=0; i<modelAllGenes.size(); i++) {
		modelAllGenes[i].deallocateqB();
	}
}

void Model::modelSingleGene::deallocateqB() {

	gsl_matrix_free(qBvarNew);
	gsl_matrix_free(qBvarBest);
}

Model::modelSingleGene::modelSingleGene() : qBvarBestAllocated(false), qBvarNewAllocated(false) { }
Model::modelSingleGene::~modelSingleGene() { }

void Model::modelSingleGene::countNumParamsToBeFitted() {
    
    int result = 0;
    
    //Selection parameters to be fitted
    for(int i=0; i<model.getMajor().getLength(); i++) {
        if(model.getMajor().getBase(i) != '-') {
            result++;
        }
    }
    
    //Epistatic parameters to be fitted
    for(unsigned int i=0; i<epiModel.size(); i++) {
        result++;
    }
    
    numParamsToBeFitted = result;
}

int Model::modelSingleGene::getNumParamsToBeFitted()  {
 
    if(numParamsToBeFitted != -1) {
        return numParamsToBeFitted;
    } else {
        return -1;
    }
    
}

void Model::countNumParamsToBeFitted() {
    int result = 0;
    for(unsigned int i=0; i<modelAllGenes.size(); i++) {
        modelAllGenes[i].countNumParamsToBeFitted();
        result += modelAllGenes[i].getNumParamsToBeFitted();
    }
    numParamsToBeFitted = result;
}

int Model::getNumParamsToBeFitted() {
   
    if(numParamsToBeFitted != -1) {
        return numParamsToBeFitted;
    } else {
        return -1;
    }
}

double Model::getBICbest() { return BICbest; }


void Model::modelSingleGene::setModel(DiploidSequence & ds) {
    model = ds;
}
void Model::modelSingleGene::setEpiModel(vector<Model::epistasis> &em) {
	epiModel = em;
}
void Model::modelSingleGene::addEpiModel(Model::epistasis & epi) {
    epiModel.push_back(epi);
}
void Model::modelSingleGene::setSelCoefsNew(vector<double> &SCN) {
    selCoefsNew = SCN;
}
void Model::modelSingleGene::setSelCoefsBest(vector<double> &SCB) {
    selCoefsBest = SCB;
}
void Model::modelSingleGene::setSelCoefsNew(int index, double value) {
    selCoefsNew[index] = value;
}
void Model::modelSingleGene::setSelCoefsBest(int index, double value) {
    selCoefsBest[index] = value;
}
DiploidSequence Model::modelSingleGene::getModel() { return model; }
vector<Model::epistasis> Model::modelSingleGene::getEpiModel() {
    return epiModel;
}
void Model::modelSingleGene::setEpiCoefsNew(int index, double value) {
    epiCoefsNew[index] = value;
}
void Model::modelSingleGene::setEpiCoefsBest(int index, double value) {
    epiCoefsBest[index] = value;
}
vector<double> Model::modelSingleGene::getSelCoefsNew() {
    return selCoefsNew;
}
double Model::modelSingleGene::getSelCoefsNew(int index) { return selCoefsNew[index]; }

int Model::getNumGenes() { return (int) modelAllGenes.size(); }
int Model::getNtNew() { return NtNew; }
int Model::getNtBest() { return NtBest; }
void Model::setNtNew(int n) { NtNew = n; }
void Model::setNtBest(int n) { NtBest = n; }
void Model::updateBICbest(double BICNew) {
    
	if(BICNew < (BICbest)) { 	

		//cout << "New best BIC found: " << BICNew << ".\tNew Nt: " << NtNew << ".\n"; //TEMP
		NtBest = NtNew;
		for(unsigned int i=0; i<modelAllGenes.size(); i++) {

			modelAllGenes[i].setSelCoefsBestToNew();
			modelAllGenes[i].setEpiCoefsBestToNew();

			vector<double> qBmeanNew = modelAllGenes[i].getqBmeanNew();
			modelAllGenes[i].setqBmeanBest(qBmeanNew);

			gsl_matrix* qBvarNewGene = modelAllGenes[i].getqBvarNew();
			modelAllGenes[i].setqBvarBest(qBvarNewGene); //Uses memcpy
		}

		BICbest = BICNew;

		//Update acceptance rates
		accepts[currentlyBeingUpdated] = accepts[currentlyBeingUpdated] + 1;
		acceptsOverall++;
	} else {
		NtNew = NtBest;
		for(unsigned int i=0; i<modelAllGenes.size(); i++) {
			modelAllGenes[i].setSelCoefsNewToBest();
			modelAllGenes[i].setEpiCoefsNewToBest();
			
			vector<double> qBmeanBest = modelAllGenes[i].getqBmeanBest();
			modelAllGenes[i].setqBmeanNew(qBmeanBest);

			gsl_matrix* qBvarBestGene = modelAllGenes[i].getqBvarBest();
			modelAllGenes[i].setqBvarNew(qBvarBestGene); //Uses memcpy
		}
	}

	tries[currentlyBeingUpdated] = tries[currentlyBeingUpdated] + 1;
	triesOverall++;
}
void Model::modelSingleGene::setSelCoefsBestToNew() {
    //cout << "Inside selCoefsBestToNew. Before: "; printDoubleVector(selCoefsBest);
    selCoefsBest = selCoefsNew;
    //cout << "Inside selCoefsBestToNew. After: "; printDoubleVector(selCoefsBest);
}
void Model::modelSingleGene::setSelCoefsNewToBest() {
    selCoefsNew = selCoefsBest;
}
void Model::modelSingleGene::setEpiCoefsBestToNew() {
    epiCoefsBest = epiCoefsNew;
}
void Model::modelSingleGene::setEpiCoefsNewToBest() {
    epiCoefsNew = epiCoefsBest;
}

bool Model::modelSingleGene::epistasisPresent() {
    if(epiModel.size() > 0) { return true; }
    else { return false; }
}
bool Model::epistasisPresent() {
    bool result = false;
    for(unsigned int i=0; i<modelAllGenes.size();i++) {
        if(modelAllGenes[i].epistasisPresent() == true) {
            result = true;
            break;
        }
    }
    return result;
}

int Model::modelSingleGene::getEpiModelLength() {
    return (int) epiModel.size();
}
double Model::modelSingleGene::getEpiCoefsBest(int index) {
    return epiCoefsBest[index];
}
double Model::modelSingleGene::getEpiCoefsNew(int index) {
    return epiCoefsNew[index];
}

void Model::modelSingleGene::setEpiCoefsBest(vector<double> &ECB) { epiCoefsBest = ECB; }
void Model::modelSingleGene::setEpiCoefsNew(vector<double> &ECN) { epiCoefsNew = ECN; }
vector<double> Model::modelSingleGene::getEpiCoefsNew() { return epiCoefsNew; }
vector<double> Model::modelSingleGene::getEpiCoefsBest() { return epiCoefsBest; }

void Model::modelSingleGene::setqBmeanBest(vector<double> & q) { qBmeanBest = q; }
void Model::modelSingleGene::setqBmeanNew(vector<double> & q) { qBmeanNew = q; }
vector<double>& Model::modelSingleGene::getqBmeanBest() { return qBmeanBest; }
vector<double>& Model::modelSingleGene::getqBmeanNew() { return qBmeanNew; }


//Epistasis related methods
void Model::epistasis::setDim(int d) { dim = d; }
int Model::epistasis::getDim() { return (int) positions.size(); }
void Model::epistasis::addPositions(int pos) { positions.push_back(pos); }
void Model::epistasis::setPositions(vector<int> &pos) { positions = pos; }
vector<int> Model::epistasis::getPositions() { return positions; }
int Model::epistasis::getPosition(int index) { return positions[index]; }


void Model::addModelSingleGene(Model::modelSingleGene &MSG) {
    modelAllGenes.push_back(MSG);
}

void Model::setModelSingleGene(int index, modelSingleGene &MSG) {

	modelAllGenes[index] = MSG;
}

void Model::initialiseAndRandomiseCoefs(gsl_rng *r, vector<vector<Sequence> > &haps) {

	//Number of bottleneck sizes to be fitted (one for now)
	numNtToBeFitted = 1;

	//Initialise qB frequencies to the same value
	numqBtoBeFitted = 0;

//Temporarily commented out - assuming qB fixed by qBpre
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
	for(int i=0; i<((int) modelAllGenes.size()); i++) { //Loop over genes
        
		Sequence currentGeneMajor = modelAllGenes[i].getModel().getMajor();
        
		vector<double> dummy(currentGeneMajor.getLength());
		modelAllGenes[i].setSelCoefsBest(dummy);
		modelAllGenes[i].setSelCoefsNew(dummy);

//Not necessary anymore		//Set selection present to false, then alter if selection found
//Remove when convenient		modelAllGenes[i].setSelectionPresent(false);
	
        
		//Randomise selection coefficients
		for(int j=0; j<currentGeneMajor.getLength();j++) { //Loop over positions
			if(currentGeneMajor.getBase(j) != '-') {

				//Actual code
				double selStrength = 0;
				modelAllGenes[i].setSelCoefsBest(j, selStrength);
				modelAllGenes[i].setSelCoefsNew(j, selStrength);

				//Add coefficient to selCoefsToBeFitted
				vector<int> selCoefsPos = {i,j}; //{gene, pos}
				selCoefsToBeFitted.push_back(selCoefsPos);

//This shouldn't be necessary any more				//Update selectionPresent
//Taken care of somewhere else. Remove when convenient				modelAllGenes[i].setSelectionPresent(true);
				
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
            
			for(int j=0; j<((int) epiVec.size()); j++) {
				double epiStrength = 0;
				modelAllGenes[i].setEpiCoefsBest(j, epiStrength);
				modelAllGenes[i].setEpiCoefsNew(j, epiStrength);


				//Add epi to epiCoefsToBeFitted
				vector<int> epiCoefsNum = {i, j}; //{gene,coefficient #}
				epiCoefsToBeFitted.push_back(epiCoefsNum);
			}
		}
	}
    
	//Initialise Nt
	NtBest = 100;
	NtNew = 100;


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

void Model::update(gsl_rng *r, double delta, int maxNt) {

	currentlyBeingUpdated = -1; //-1=unassigned, 0=Nt, 1=qB, 2=sel, 3=epi

	double prob = gsl_rng_uniform(r);

	if(acceptanceRates.size() == 1) { //Single parameter to be updated

		if(numNtToBeFitted > 0) {
			currentlyBeingUpdated = 0; //Updating index 0 of acceptanceRates
			updateNt(r,delta,maxNt);		
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
				updateNt(r, delta, maxNt);
			} else {
				currentlyBeingUpdated = 1;
				updateqB(r, delta);
			}
		} else if(numNtToBeFitted>0 && selCoefsToBeFitted.size()>0) {

			if(prob < acceptanceRates[0]) {
				currentlyBeingUpdated = 0;
				updateNt(r, delta, maxNt);
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
				updateNt(r, delta, maxNt);
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
				updateNt(r, delta, maxNt);
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
			updateNt(r, delta, maxNt);
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

void Model::updateNt(gsl_rng *r, double delta, int maxNt) {
    int ceilDelta = (int) ceil(delta);
	int factor;
	if(NtBest < 100) {
		factor = 1;
	} else if(NtBest < 500) { //Between 100 and 500, jump in increments of maximum 5 (min 1)
		factor = gsl_rng_uniform_int(r,5)+1;
	} else if(NtBest < 1000) { //Between 500 and 1000, increment by maximum 10
		factor = gsl_rng_uniform_int(r,10)+1;
	} else {
		factor = gsl_rng_uniform_int(r,50)+1;
	}
    if(gsl_rng_uniform(r) > 0.5) {
        if(NtBest <= (maxNt-ceilDelta*factor)) { NtNew+= ceilDelta*factor; }
    } else {
        if(NtBest >= (1+ceilDelta*factor)) { NtNew-= ceilDelta*factor; }
    }
}

void Model::updateCoefs(gsl_rng *r, double delta) {

    if(gsl_rng_uniform(r) < 0.50) {
        updateSelCoefs(r, delta);
	} else {
        updateEpiCoefs(r, delta);
    }
}

void Model::updateSelCoefs(gsl_rng *r, double delta) {

	if(selCoefsToBeFitted.size() > 0) { //Only try updating if there is something that needs updating
		int pickedCoef = gsl_rng_uniform_int(r,selCoefsToBeFitted.size());
		vector<double> selCoefsNew = modelAllGenes[selCoefsToBeFitted[pickedCoef][0]].getSelCoefsNew();

		//Keep finding new selCoefs until all coefs are within range (-10,10)
		bool withinRange = false;		
		while(withinRange == false) {
		
			double newValue = selCoefsNew[selCoefsToBeFitted[pickedCoef][1]] + (gsl_rng_uniform(r)-0.5)*delta;
			if(newValue > selectionCap || newValue < -selectionCap) {
				withinRange = false;
			} else{
				withinRange = true;
				selCoefsNew[selCoefsToBeFitted[pickedCoef][1]] = newValue;
			}
		}
		modelAllGenes[selCoefsToBeFitted[pickedCoef][0]].setSelCoefsNew(selCoefsNew);
	}
    
}

//This one might not work when including data for empty genes..
//If ever used again, update it
void Model::updateqB(gsl_rng *r, double delta) {
   
	//Pick random gene
	int numGenes = (int) modelAllGenes.size();
	int pickedGene = (int) gsl_rng_uniform_int(r, numGenes);
        
	if(gsl_rng_uniform(r) < 0.5) { //Update frequencies
		vector<double> qBmeanNew;
		vector<double> qBold = modelAllGenes[pickedGene].getqBmeanBest();
		addRandom(qBold,qBmeanNew,delta,r);

		modelAllGenes[pickedGene].setqBmeanNew(qBmeanNew);
	} else { //Update variance

		gsl_matrix * qBvarBestGene = modelAllGenes[pickedGene].getqBvarBest();
		gsl_matrix * qBvarNewGene = modelAllGenes[pickedGene].getqBvarNew();
		int dim = (int) qBvarBestGene->size1;

		bool positiveDefinite = false;

		while(positiveDefinite == false) {
			int index1 = gsl_rng_uniform_int(r,dim);
			int index2 = index1; // gsl_rng_uniform_int(r,dim);
			double oldValue = gsl_matrix_get(qBvarBestGene, index1, index2);
			double newValue;
			do {
				newValue = oldValue + delta*(gsl_rng_uniform_pos(r) - 0.5);
			} while(newValue < 10e-10);

			if(index1 == index2) {
				gsl_matrix_set(qBvarNewGene, index1, index1, newValue);
			} else {
				gsl_matrix_set(qBvarNewGene, index1, index2, newValue);
				gsl_matrix_set(qBvarNewGene, index2, index1, newValue);	
			}
		
			//Check if positive-definite
			gsl_matrix* qBvarCholesky = gsl_matrix_alloc(dim,dim);
			gsl_matrix_memcpy(qBvarCholesky, qBvarNewGene); //Copy into qBvarCholesky as matrix is modified by action below

			//The function below is deprecated by gsl. New function is gsl_linalg_cholesky_decomp1(), but this version of gsl doesn't support that.
			//Hence we use the deprecated version
			gsl_set_error_handler_off();
			int status = gsl_linalg_cholesky_decomp(qBvarCholesky);
			if(status == 0) { //Success indicates that matrix is symmetric positive definite
				positiveDefinite = true;
			} else {
					
				gsl_matrix_memcpy(qBvarNewGene, qBvarBestGene); //Not positive definite, so revert changes
			}
	
		}	

		//Shouldn't need to set qBvarNew as it is a pointer

	}
    
}
void Model::updateEpiCoefs(gsl_rng *r, double delta) {
    
    if(epistasisPresent() == true) { //Checks if there is any epistasis at all across all genes
        
        //int ceilDelta = (int) ceil(delta);
        
        bool changeOccurred = false;
        while(changeOccurred == false) {
            
            //Pick random gene
            int numGenes = (int) modelAllGenes.size();
            int pickedGene = (int) gsl_rng_uniform_int(r, numGenes);
            
          //Pick random epistatic coefficient
            int epiLength = modelAllGenes[pickedGene].getEpiModelLength();
            if (epiLength > 0) { //Checks if there is epistasis on the picked gene
                int pickedEpiCoef = (int) gsl_rng_uniform_int(r, epiLength);
                
                double oldValue = modelAllGenes[pickedGene].getEpiCoefsBest(pickedEpiCoef);
		double newValue = 0;
		
		//Keep finding new epistasis value until it is within range (-10,10)
		bool withinRange = false;
		while(withinRange == false) {
	
			newValue = oldValue + (gsl_rng_uniform(r)*delta)-delta*0.5;
		
			if(newValue > selectionCap || newValue < -selectionCap) {
				withinRange = false;
			} else {
				withinRange = true;
                    		modelAllGenes[pickedGene].setEpiCoefsNew(pickedEpiCoef, newValue);
				changeOccurred = true;
			}
		}
            }
        }
    }
}

vector<double> Model::computeHapFitTransmission(int gene, vector<Sequence> &fullHaps) {
	
	return modelAllGenes[gene].computeHapFitTransmission(fullHaps);
}

vector<double> Model::modelSingleGene::computeHapFitTransmission(vector<Sequence> &fullHaps) {

	vector<double> hapFit(fullHaps.size(),0);
	Sequence modelSeq = model.getMajor();
       

	//Check for epistasis first so we can enforce haplotype fitness max limit in selection below
	if(epiModel.size() > 0) { //Possibility of epistasis
        
		for(unsigned int j=0;j<fullHaps.size();j++) { //Loop over haplotypes
            
			for(unsigned int k=0; k<epiModel.size(); k++) { //Loop over epistatic models
                   
				vector<int> positions = epiModel[k].getPositions();
				bool epistasisPresent = true; //Assume it is true, then check if the haplotype fits epistatic model
				for(unsigned int l=0; l<positions.size();l++) { //Loop over positions in epi model
                    
					if(fullHaps[j].getBase(positions[l]) != modelSeq.getBase(positions[l])) {
						epistasisPresent = false;
						break;
					}
				}
                   
				if(epistasisPresent == true) {
					hapFit[j] += epiCoefsNew[k];
				}
			}
		}
	}


	//Check for selection
	for(unsigned int j=0;j<fullHaps.size();j++) { //Loop over haplotypes
		for(int k=0;k<fullHaps[j].getLength();k++) { //For each locus in haplotype
               
			if(modelSeq.getBase(k) != '-') { //Locus under selection
				if(fullHaps[j].getBase(k) == modelSeq.getBase(k)) { //Haplotype under selection
					hapFit[j] = hapFit[j] + selCoefsNew[k];
				}
			}
		}

//// This is no longer in use as it impacted inference dynamics in a bad way
//
//		//Ensure that no haplotype fitness are larger than +/- 10
//		if(hapFit[j] > 10) {
//			hapFit[j] = 10;
//		} else if(hapFit[j] < -10) {
//			hapFit[j] = -10;
//		}

	}
       
       
	return hapFit;
}

Model::modelSingleGene Model::getModelSingleGene(int index) {
    return modelAllGenes[index];
}

vector<double> Model::modelSingleGene::getSelCoefsBest() { return selCoefsBest; }


void Model::setqBmeanNew(vector<vector<double> > &qB) {

	for(unsigned int i=0; i<modelAllGenes.size(); i++) {
		modelAllGenes[i].setqBmeanNew(qB[i]);
	}			
}

void Model::setqBmeanBest(vector<vector<double> > &qB) {

	for(unsigned int i=0; i<modelAllGenes.size(); i++) {
		modelAllGenes[i].setqBmeanBest(qB[i]);
	}			
}

void Model::setqBmeanNew(int gene, vector<double> &qB) {
	modelAllGenes[gene].setqBmeanNew(qB);
}
void Model::setqBmeanBest(int gene, vector<double> &qB) {
	modelAllGenes[gene].setqBmeanBest(qB);
}

vector<vector<double> > Model::getqBmeanNew() {

	vector<vector<double> > result;
	for(unsigned int i=0; i<modelAllGenes.size(); i++) {
		result.push_back(modelAllGenes[i].getqBmeanNew());
	}
	
	return result;
}

vector<vector<double> > Model::getqBmeanBest() {

	vector<vector<double> > result;
	for(unsigned int i=0; i<modelAllGenes.size(); i++) {
		result.push_back(modelAllGenes[i].getqBmeanBest());
	}
	
	return result;
}


vector<double>& Model::getqBmeanNew(int gene) {

	return modelAllGenes[gene].getqBmeanNew();
}

vector<double>& Model::getqBmeanBest(int gene) {

	return modelAllGenes[gene].getqBmeanBest();
}

void Model::setBICbest(double BICnew) {
	BICbest = BICnew;
}

//Check if this is a subset of M
//Assumes that this and M come from the same overall set of selection models
//i.e. that this and M have the same number of genes, same number of loci
//with same defintion of major and minor at each locus
//If used in a different context, this method will induce errors
bool Model::isASubsetOf(Model &M) {

	//Check each gene in turn
	for(unsigned i=0; i<modelAllGenes.size(); i++) {

		modelSingleGene MSG1 = modelAllGenes[i];
		modelSingleGene MSG2 = M.getModelSingleGene(i);

		//Check selection first
		DiploidSequence s1 = MSG1.getModel();
		DiploidSequence s2 = MSG2.getModel();
		for(int j=0; j<s1.getLength(); j++) {
			
			if(s1.getMajor().getBase(j) != '-') { //Locus under selection in s1
			
				if(s2.getMajor().getBase(j) == '-') { //Locus in s2 NOT under selection
					return false; //this cannot be a subset of M
				}
			}
		}
		
		//Check for epistasis
		vector<epistasis> epi1 = MSG1.getEpiModel();
		vector<epistasis> epi2 = MSG2.getEpiModel();
		for(unsigned int j=0; j<epi1.size(); j++) { //Loop over epistasis in MSG1
			
			vector<int> epi1Pos = epi1[j].getPositions();
			bool found = false;
	
			//Check if this instance of epistasis exists in MSG2 - if not, then can't be a subset
			for(unsigned int k=0; k<epi2.size(); k++) {

				vector<int> epi2Pos = epi2[k].getPositions();

				if(epi1Pos.size() == epi2Pos.size()) { //Only compare epistasis of same length
					if(identicalIntVectors(epi1Pos,epi2Pos) == true) { //Check if epistasis positions are identical
						found = true; break;
					}
				}
			}

			if(found == false) { return false; } //Cannot be a subset if an instance of epistasis in MSG1 wasn't found in MSG2
		}
		
	}


	return true;
}


void Model::print() {

	for(unsigned int i=0; i<modelAllGenes.size(); i++) {

		cout << "Gene " << i << ":\n";
		modelAllGenes[i].print();
	//	cout << "qBmean: "; printDoubleVector(modelAllGenes[i].getqBmeanBest());
	//	cout << "qBvar: "; printMatrixMathematica(modelAllGenes[i].getqBvarBest());
	}
	cout << "Nt: " << NtBest << endl;
}

void Model::printNew() {

	for(unsigned int i=0; i<modelAllGenes.size(); i++) {

		cout << "Gene " << i << ":\n";
		modelAllGenes[i].printNew();
	}
	cout << "Nt: " << NtNew << "\n";
}

void Model::modelSingleGene::print() {

	cout << "SelModel:\n";
	model.print();
	if(selCoefsBest.size() > 0) { //Only print out values if model has been initialised
		printDoubleVector(selCoefsBest);
	}

	if(epiModel.size() > 0) {

		cout << "Epistasis:\n";

		for(unsigned int i=0; i<epiModel.size(); i++) {
			epiModel[i].print();
			if(epiCoefsBest.size() > 0) { //Only print values if model has been initialised
				cout << "Value: " << epiCoefsBest[i] << "\n";
			}
		}
	}
}

void Model::modelSingleGene::printNew() {

	cout << "SelModel:\n";
	model.print();
	if(selCoefsNew.size() > 0) { //Only print out values if model has been initialised
		printDoubleVector(selCoefsNew);
	}

	if(epiModel.size() > 0) {

		cout << "Epistasis:\n";

		for(unsigned int i=0; i<epiModel.size(); i++) {
			epiModel[i].print();
			if(epiCoefsNew.size() > 0) { //Only print values if model has been initialised
				cout << "Value: " << epiCoefsNew[i] << "\n";
			}
		}
	}
}


void Model::epistasis::print() {

	printIntVector(positions);
}


void Model::resetAccepts() {

	//Reset accepts for Nt, qB, sel and epi
	for(unsigned int i=0; i<accepts.size(); i++) { 

		accepts[i] = 0;
		tries[i] = 0;
	}
	acceptsOverall = 0;
	triesOverall = 0;
}

double Model::updateAcceptanceRates() {

	acceptanceRates.clear();
	for(unsigned int i=0; i<accepts.size(); i++) {

		if(accepts[i] == 0) {
			acceptanceRates.push_back(0);
		} else {
			acceptanceRates.push_back((accepts[i]+0.)/(tries[i]+0.));
		}
	}

	acceptanceRateOverall = (acceptsOverall+0.)/(triesOverall+0.);
	
//	cout << "New accepts: "; printIntVector(accepts);
//	cout << "New tries: "; printIntVector(tries);
//	cout << "New acceptance rates: "; printDoubleVector(acceptanceRates);

	resetAccepts(); //Reset the counts, not the rates
	rescaleMin(acceptanceRates, 0.1); //Makes the rates into probabilities, but with a minimum value of 0.01

//	cout << "New probs: "; printDoubleVector(acceptanceRates);

	return acceptanceRateOverall;
}


gsl_matrix* Model::getqBvarBest(int gene) {

	return modelAllGenes[gene].getqBvarBest();
}
gsl_matrix* Model::getqBvarNew(int gene) {

	return modelAllGenes[gene].getqBvarNew();
}
void Model::setqBvarBest(int gene, gsl_matrix* var) { modelAllGenes[gene].setqBvarBest(var); }
void Model::setqBvarNew(int gene, gsl_matrix* var) { modelAllGenes[gene].setqBvarNew(var); }


void Model::modelSingleGene::setqBvarBest(gsl_matrix* var) {
	if(qBvarBestAllocated == false) {
		int dim = var->size1;
		qBvarBest = gsl_matrix_alloc(dim,dim);
		qBvarBestAllocated = true;
	}
	gsl_matrix_memcpy(qBvarBest, var);
}
void Model::modelSingleGene::setqBvarNew(gsl_matrix* var) {
	if(qBvarNewAllocated == false) {
		int dim = var->size1;
		qBvarNew = gsl_matrix_alloc(dim,dim);
		qBvarNewAllocated = true;
	}
	gsl_matrix_memcpy(qBvarNew, var);
}
gsl_matrix* Model::modelSingleGene::getqBvarBest() { return qBvarBest; }
gsl_matrix* Model::modelSingleGene::getqBvarNew() { return qBvarNew; }

void Model::setSelCoefsToBeFitted(vector<vector<int> >& s) { selCoefsToBeFitted = s; }

bool Model::modelSingleGene::getSelectionPresent() { return selectionPresent; }
bool Model::getSelectionPresent(int gene) { return modelAllGenes[gene].getSelectionPresent(); }
void Model::modelSingleGene::setSelectionPresent(bool b) { selectionPresent = b; }
void Model::setSelectionCap(double value) { selectionCap = value; }
