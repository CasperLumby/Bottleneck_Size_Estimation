//
//  analyserPHMGWH.cpp
//

#include "analyserPHMGWH.hpp"
#include "analyserPHWH.hpp"
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

AnalyserPHMGWH::AnalyserPHMGWH() { } //Constructor does nothing
AnalyserPHMGWH::~AnalyserPHMGWH() { } //Destructor does nothing

//Data *d pointer passed by copy, not by reference. Can't change original data.
void AnalyserPHMGWH::loadData(Data *d) {
    
	cout << "Loading data for full haplotypes. \n";
	data = dynamic_cast<DataPHMGgen *>(d);
}
 
void AnalyserPHMGWH::runAnalysis(AnaParam *ap) {
    
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
						
	
	double C = ap->getC();
		
	//Obtain full haplotypes
	vector<vector<Sequence> > fullHaps;
	cout << "Number of genes: " << getNumGenes() << "\n";
	for(int g=0;g<getNumGenes();g++) {
		if(data->getGenePresent(g) == true) {
			DataPHgen* d = data->getGene(g);

			cout << "Printing data for gene " << g << ":\n";
			d->print();

			//Generate full haps from partial haps
			vector<Sequence> phs = d->getSequences(); //Get partial haplotypes
			vector<Sequence> fhs; //Full haplotypes for current gene
			

//			cout << "Printing partial haplotypes before constructing full haplotypes:\n";
//			for(unsigned int j=0; j<phs.size(); j++) {
//				phs[j].print();
//			}


			
			AnalyserPH aPH;
			aPH.constructFullHaps(&phs, &fhs);
			cout << "Minimum set of full haplotypes: \n";
			cout << "fhs.size(): " << fhs.size() << "\n";
			for(unsigned int j=0; j<fhs.size();j++) {
				fhs[j].print();
				cout << j+1 << " of " << fhs.size() << "\n";
			}
			d->computeContribs(fhs); //Updates contribs wrt full haplotypes
			fullHaps.push_back(fhs); //Add to full haps for all genes

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

//	//Infer frequencies q as a mean vector for each time point
	vector<vector<vector<double> > > qMean; //Levels: gene, time point, haplotype
//	vector<gsl_matrix*> qBvar;
	vector<vector<double> > LqB; //Likelihoods for qB optimisation. Levels: gene, time point

	for(int g=0;g<getNumGenes();g++) {
		if(data->getGenePresent(g) == true) {
			DataPHgen* d = data->getGene(g);

			cout << "Inferring full haplotype frequencies for gene " << g << ":\n";
			AnalyserPHWH aPHWH;
			aPHWH.loadData(d);
            
			//Compute qMean for all time points
			vector<vector<double> > qMeanGene; //Levels: time point, haplotype
//			gsl_matrix* qBvarGene = gsl_matrix_alloc(fullHaps[g].size(), fullHaps[g].size()); //Allocate here, otherwise can't be passed back
			vector<double> LqGene;
			cout << "List of full haps for gene " << g << ":\n";
			for(unsigned int i=0; i<fullHaps[g].size(); i++) {

				fullHaps[g][i].print();
			}
			aPHWH.optimiseqMean((int)fullHaps[g].size(), r, C, &qMeanGene, &LqGene); 


			cout << "Printing optimised qMean for gene " << g << ":\n";
			for(unsigned int i=0; i<qMeanGene.size(); i++) { //Loop over time points
				cout << "Time point " << i << ": "; printDoubleVector(qMeanGene[i]);

			}
//			cout << "Printing optimised qB variance for gene " << g << ":\n";
//			printMatrixMathematica(qBvarGene);

			qMean.push_back(qMeanGene);
//			qBvar.push_back(qBvarGene);
			LqB.push_back(LqGene);
	
		} else { //I.e. no data for current gene, so add empty variables

			vector<vector<double> > qMeanGene;
			qMean.push_back(qMeanGene);
//			gsl_matrix* qBvarGene = gsl_matrix_alloc(1,1);
//			gsl_matrix_set(qBvarGene, 0,0, 0); //Define a 1x1 matrix with element 0 as an "empty" matrix - no other way to do this, as can't have it unassigned or 0 size
//			qBvar.push_back(qBvarGene);
			vector<double> LqGene;
			LqB.push_back(LqGene);
		}
	}


	//Print qMean
	for(unsigned int i=0; i<qMean.size(); i++) { //Loop over genes

		ss.str(""); //Clear stringstream
		ss << filePathFolder << "qMean_gene" << i << ".dat";
		filePath = ss.str();
		outputFile.open(filePath.c_str());
	
		for(unsigned int j=0; j<qMean[i].size(); j++) { //Loop over time points

			outputFile << j << " ";
			for(unsigned int k=0; k<qMean[i][j].size(); k++) {

				outputFile << qMean[i][j][k] << " ";

			}
			outputFile << "\n";
		}
		outputFile.close();
	}
//	
//
//	//Print qBvar
//	ss.str(""); //Clear stringstream
//	ss << filePathFolder << "qBvarPre.dat";
//	filePath = ss.str();
//	outputFile.open(filePath.c_str());
//	for(unsigned int i=0; i<qBvar.size(); i++) { //Loop over genes
//		
//		if(data->getGenePresent(i) == true) {
//			string qBvarGeneString = printMatrixMathematicaToString(qBvar[i]);
//		
//			outputFile << qBvarGeneString << "\n";
//		
//		} else {
//			outputFile << "\n";
//		}
//	}
//	outputFile.close();
//
//	//Print LqB
//	ss.str(""); //Clear stringstream
//	ss << filePathFolder << "LqB.dat";
//	filePath = ss.str();
//	outputFile.open(filePath.c_str());
//	for(unsigned int i=0; i<LqB.size(); i++) { //Loop over genes
//		
//		if(data->getGenePresent(i) == true) {
//		
//			outputFile << LqB[i] << "\n";
//		
//		} else {
//			outputFile << "\n";
//		}
//	}
//	outputFile.close();
//
}

int AnalyserPHMGWH::getNumGenes() { return AnalyserPHMGWH::data->getNumGenes(); }




