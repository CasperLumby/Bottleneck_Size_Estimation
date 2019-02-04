//
//  analyserPHWH.cpp
//
//

#include "analyserPHWH.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include "misc.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <limits>
#include <math.h>
#include <gsl/gsl_linalg.h>	
#include <algorithm>

using namespace std;

AnalyserPHWH::AnalyserPHWH() { } //Constructor does nothing
AnalyserPHWH::~AnalyserPHWH() { } //Destructor does nothing

void AnalyserPHWH::loadData(Data *d) {
    
    //cout << "Loading data for partial haplotypes. \n";
    data = dynamic_cast<DataPHgen *>(d);
}

void AnalyserPHWH::runAnalysis(AnaParam *ap) {
    
	//This method needs to get updated (if ever in use again)
	cout << "Running analysis for partial haplotypes. - CURRENTLY NOT IMPLEMENTED!\n";
    
}

void AnalyserPHWH::optimiseqMean(int numFHs, gsl_rng *r, int C, vector<vector<double> >* qMean, vector<double>* pL) {

	int dim = numFHs;
	qMean->clear();
	pL->clear(); //pointer to Likelihoods


	for(int t=0; t<data->getNumOfTimePoints(); t++) { //Loop over time points, assuming that all SOPHs have the same time points

		cout << "\n\nOptimising qMean for time point " << t << ":\n";
		vector<double> qMeanTimePoint; //qMean for time point t

		//Initialise mean to random values, then renormalise
		for(int i=0; i<dim; i++) {

			qMeanTimePoint.push_back(gsl_rng_uniform_pos(r));
		}
		rescale(qMeanTimePoint);

  
		//Define "new" parameters
		vector<double> qMeanTimePointNew = qMeanTimePoint;

	
		double bestL = - std::numeric_limits<double>::max();
		double delta = 1;
		int accepts = 0;
		int trys = 0;
		int nUpdate = 1000; //Frequency for updating acceptance rate


		int k=0; //Counter
		while(delta > 0.000001) {
			k++;
			if(k%nUpdate==0) {

				double acceptanceRate=(accepts+0.)/(trys+0.);
				//cout << "Acceptance rate: " << acceptanceRate << "\tdelta: " << delta << "\n";


				double tempDelta=delta*(0.95+acceptanceRate);
				double tempnUpdate = ceil(nUpdate*(0.95+acceptanceRate));
				if(tempDelta < 10) { // If delta>10 computations get out of ahnd
					delta = tempDelta;
					if(tempnUpdate > 100 && tempnUpdate <= 1000) { nUpdate = tempnUpdate; }

					cout << "Delta: " << delta << "\tnUpdate: " << nUpdate << "\tL: " << bestL << "\n";
				}

				//Reset
				accepts=0;
				trys=0;
			}
	
			//Update mean
			addRandom(qMeanTimePoint,qMeanTimePointNew,delta,r); //Create a new vector newFreqs from adding a small delta*rndm(double) to a random index

//			double newL = computeLmean(&qMeanTimePointNew, r, C, t, false); 
			double newL = computeLmeanDirMult(&qMeanTimePointNew, C, t, false); 

			if(newL > bestL) {

				qMeanTimePoint = qMeanTimePointNew;

				bestL = newL;
				accepts++;
		
			} else {
				//Revert back any changes
				qMeanTimePointNew = qMeanTimePoint;
			}	
		
			trys++;
		}

		qMean->push_back(qMeanTimePoint);
		pL->push_back(bestL);
	}
} 



//Input mean of full dimensionality	
double AnalyserPHWH::computeLmean(vector<double>* qMean, gsl_rng* r, double C, int timeIndex, bool print) {

	double L = 0;
	int numSOPHs = data->getNumOfSOPHs();
	int dimFH = (int) qMean->size();


	for(int i=0; i<numSOPHs; i++) { //Loop over sets of partial haps

		DataPHgen::setOfPartialHaps* SOPHi = data->getSOPH(i);

		if(print == true) {
			cout << "SOPH " << i << "\n";
		}
		int dimPH = SOPHi->getNumOfPartialHaps();
        	vector<int> obs = SOPHi->getObs(timeIndex); //This might be a slow process. Consider making more efficient.
		int Ntot = sumOfVector(obs);
		if(print == true) {

			cout << "obs: "; printIntVector(obs);
		}

		if(Ntot == 0) { 
	
			//This SOPH has Ntot of zero for time point i, so not informative. As such we ignore it.
			continue;
		}

		vector<vector<int> > contribsPHset;
		for(int j=0; j<dimPH; j++) {

			contribsPHset.push_back(*(SOPHi->getContribs(j)));
		}


		//Construct matrix T
		gsl_matrix* T = gsl_matrix_alloc(dimPH,dimFH);
		for(int j=0; j<dimPH; j++) { //Loop over rows 

			if(print == true) {
					
				cout << "contribs PH set " << j << ": "; printIntVector(contribsPHset[j]);
			}

			int counter = 0;
			for(int k=0; k<dimFH; k++) { //Loop over columns

				if(contribsPHset[j][counter]==k && counter < (int) contribsPHset[j].size()) { //This assumes the contribs are ordered

					gsl_matrix_set(T,j,k,1);
					counter++;
				} else {
					gsl_matrix_set(T,j,k,0);
				}
			}
		}


		/*
		 * First, transform qMean into qMean_PH=TqMean
		 */
		//Turn qMean into gsl vector
		gsl_vector *qMeanGSL = gsl_vector_alloc(dimFH);
		for(int i=0;i<dimFH;i++) {
			gsl_vector_set(qMeanGSL,i,(*qMean)[i]);
		}

		//Create TqBmean through matrix vector multiplication
		gsl_vector *TqMeanGSL = gsl_vector_alloc(dimPH);
		gsl_blas_dgemv(CblasNoTrans, 1.0, T, qMeanGSL, 0.0, TqMeanGSL);

		//Convert back into C++ vector
		vector<double> TqMean;
		vector<double> NtotTqMean;
		for(int i=0; i<dimPH; i++) {
			TqMean.push_back(gsl_vector_get(TqMeanGSL, i));
			NtotTqMean.push_back(Ntot*TqMean[i]);
		}

	
		/*
		 * Next, create the variance expression 
		 */
		//Construct epsilonNtotM(TqMean) with epsilon = (Nb+C)/(1+C)
		gsl_matrix* epsilonNtotMTqMean = gsl_matrix_alloc(dimPH,dimPH);
		constructMatrixM(TqMean, epsilonNtotMTqMean);
		double epsilon = (Ntot+C)/((double)(1+C));
		gsl_matrix_scale(epsilonNtotMTqMean, epsilon*Ntot);

		
		double LphSet = - numeric_limits<double>::max();

		if(dimPH > 1) { //Most common case, i.e. multiple partial haps in a set
		
			//Reduce dimensionality by 1 to ensure non-degeneracy
			NtotTqMean.pop_back();
			gsl_matrix_view qVarPHreducedView = gsl_matrix_submatrix(epsilonNtotMTqMean, 0,0, dimPH-1, dimPH-1); //Remove row k and column k
			gsl_matrix * qVarPH = gsl_matrix_alloc(dimPH-1, dimPH-1);
			gsl_matrix_memcpy(qVarPH, &(qVarPHreducedView.matrix)); //Copy the reduced matrix view into a "normal" matrix varXaPH
			vector<double> obsReduced;
			for(int i=0; i<dimPH-1; i++) {
				obsReduced.push_back(obs[i]);
	        	}


			if(dimPH>2) { //Multivariate approach

				LphSet = logMultivariateNormalPDF(obsReduced,NtotTqMean, qVarPH);
		
			} else if(dimPH==2) {

				//Univariate case
				double qMeanUni = NtotTqMean[0];
				double qVarUni = gsl_matrix_get(qVarPH,0,0);
				double stDevqVarUni = sqrt(qVarUni);

				LphSet = -(obs[0]-qMeanUni)*(obs[0]-qMeanUni)/(2*qVarUni) - log(stDevqVarUni*sqrt(2*M_PI)); //log L
	
			}

			//Deallocate variables created in this scope			
			gsl_matrix_free(qVarPH);

		} else if (dimPH == 1) { 
			//Case where there is a single partial haplotype in a set.
			//This can be thought of as a case with a single observed partial haplotype
			//and one (or several) unobserved haplotypes.
			//If all loci in the partial haplotype set are variants (i.e. not monomorphic)
			//then there MUST exists at least one other haplotype containing the unobserved alleles,
			//however, this (or these) haplotype(s) are not observed.
			//As such, WLOG we can consider this a case of dimPH=2, which gets reduced to a one dimensional system
			//under reduction, i.e. similar to the dimPH==2 case of above.
				
			double qMeanUni = NtotTqMean[0];
			double qVarUni = gsl_matrix_get(epsilonNtotMTqMean,0,0);
			double stDevqVarUni = sqrt(qVarUni);

			LphSet = -(obs[0]-qMeanUni)*(obs[0]-qMeanUni)/(2*qVarUni) - log(stDevqVarUni*sqrt(2*M_PI)); //log L
		
	
		} else {

			cout << "ERROR in computeLmean: dimPH is <1\n";
		}


		//Check if L is still negative infinity (shouldn't happen)
		if(LphSet == -numeric_limits<double>::max() || ::isnan(LphSet) == true) { 
		
				cout << "ERROR in AnalyserPH::computeLmean: Likelihood is negative infinity.\n";
				exit(1);
		}
		else{

			L += LphSet;
		}

		//Clean up
		gsl_vector_free(qMeanGSL);
		gsl_vector_free(TqMeanGSL);
		gsl_matrix_free(T);
		gsl_matrix_free(epsilonNtotMTqMean);

	}

	return L;
} 



//Input mean of full dimensionality	
double AnalyserPHWH::computeLmeanDirMult(vector<double>* qMean, double C, int timeIndex, bool print) {

	double L = 0;
	int numSOPHs = data->getNumOfSOPHs();
	int dimFH = (int) qMean->size();


	for(int i=0; i<numSOPHs; i++) { //Loop over sets of partial haps

		DataPHgen::setOfPartialHaps* SOPHi = data->getSOPH(i);

		if(print == true) {
			cout << "SOPH " << i << "\n";
		}
		int dimPH = SOPHi->getNumOfPartialHaps();
        	vector<int> obs = SOPHi->getObs(timeIndex); //This might be a slow process. Consider making more efficient.
		int Ntot = sumOfVector(obs);
		if(print == true) {

			cout << "obs: "; printIntVector(obs);
		}

		if(Ntot == 0) { 
	
			//This SOPH has Ntot of zero for time point i, so not informative. As such we ignore it.
			continue;
		}

		vector<vector<int> > contribsPHset;
		for(int j=0; j<dimPH; j++) {

			contribsPHset.push_back(*(SOPHi->getContribs(j)));
		}


		//Construct matrix T
		gsl_matrix* T = gsl_matrix_alloc(dimPH,dimFH);
		for(int j=0; j<dimPH; j++) { //Loop over rows 

			if(print == true) {
					
				cout << "contribs PH set " << j << ": "; printIntVector(contribsPHset[j]);
			}

			int counter = 0;
			for(int k=0; k<dimFH; k++) { //Loop over columns

				if(contribsPHset[j][counter]==k && counter < (int) contribsPHset[j].size()) { //This assumes the contribs are ordered

					gsl_matrix_set(T,j,k,1);
					counter++;
				} else {
					gsl_matrix_set(T,j,k,0);
				}
			}
		}


		/*
		 * Construct alpha=C*qMean_PH = C*T*qMean
		 */
		//Turn qMean into gsl vector
		gsl_vector *qMeanGSL = gsl_vector_alloc(dimFH);
		for(int i=0;i<dimFH;i++) {
			gsl_vector_set(qMeanGSL,i,(*qMean)[i]);
		}

		//Create TqBmean through matrix vector multiplication
		gsl_vector *TqMeanGSL = gsl_vector_alloc(dimPH);
		gsl_blas_dgemv(CblasNoTrans, 1.0, T, qMeanGSL, 0.0, TqMeanGSL);

		//Convert back into C++ vector
		vector<double> alpha;
		for(int i=0; i<dimPH; i++) {
			alpha.push_back(C*gsl_vector_get(TqMeanGSL, i));
		}

		
		double LphSet = - numeric_limits<double>::max();

		if(dimPH > 1) { //Most common case, i.e. multiple partial haps in a set
			
			if(dimPH>2) { //Multivariate approach, i.e. Dirichlet Multinomial likelihood

				LphSet = logDirMultProb(alpha, obs, Ntot);

		
			} else if(dimPH==2) { //Univariate case, i.e. beta binomial

				//alpha, beta, x, n
				LphSet = logBetaBinProb(alpha[0], alpha[1], obs[0], Ntot);

	
			}


		} else if (dimPH == 1) { 
			//Case where there is a single partial haplotype in a set.
			//This can be thought of as a case with a single observed partial haplotype
			//and one (or several) unobserved haplotypes.
			//If all loci in the partial haplotype set are variants (i.e. not monomorphic)
			//then there MUST exists at least one other haplotype containing the unobserved alleles,
			//however, this (or these) haplotype(s) are not observed.
			//As such, WLOG we can consider this a case of dimPH=2, i.e. similar to the dimPH==2 case of above.
				

			//alpha=Cq
			//beta=C(1-q)=C-alpha
			double beta = C - alpha[0];

			//alpha, beta, x, n
			LphSet = logBetaBinProb(alpha[0], beta, obs[0], Ntot);
		
	
		} else {

			cout << "ERROR in computeLmean: dimPH is <1\n";
		}


		//Check if L is still negative infinity (shouldn't happen)
		if(LphSet == -numeric_limits<double>::max() || ::isnan(LphSet) == true) { 
		
				cout << "ERROR in AnalyserPH::computeLmean: Likelihood is negative infinity.\n";
				exit(1);
		}
		else{

			L += LphSet;
		}

		//Clean up
		gsl_vector_free(qMeanGSL);
		gsl_vector_free(TqMeanGSL);
		gsl_matrix_free(T);

	}

	return L;
} 

