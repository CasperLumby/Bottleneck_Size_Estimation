//
//  analyserPH.cpp
//  FluTransmissionProject
//
//  Created by Casper Lumby on 30/11/2015.
//  Copyright © 2015 Casper Lumby. All rights reserved.
//

#include "analyserPH.hpp"
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
#include <gsl/gsl_cdf.h>

using namespace std;

AnalyserPH::AnalyserPH() { } //Constructor does nothing
AnalyserPH::~AnalyserPH() { } //Destructor does nothing

void AnalyserPH::loadData(Data *d) {
    
    //cout << "Loading data for partial haplotypes. \n";
    data = dynamic_cast<DataPH *>(d);
}

void AnalyserPH::runAnalysis(AnaParam *ap) {
    
	//This method needs to get updated (if ever in use again)
	cout << "Running analysis for partial haplotypes. - CURRENTLY NOT IMPLEMENTED!\n";
    
}

void AnalyserPH::constructFullHaps(vector<Sequence> *phs, vector<Sequence> *fhs) {
    
    //Changing sequences in vectors of chars to work with Chris' code
 //   cout << "Partial haps:\n";
    vector<vector<char> > haps;
    for(unsigned int i=0; i<phs->size(); i++) {
        
    //    (*phs)[i].print();
        
        vector<char> h;
        for(int j=0; j< (*phs)[i].getLength(); j++) {
            h.push_back((*phs)[i].getBase(j));
        }
        
        haps.push_back(h);
    }

    
    //Cycle through overlap steps
    vector<vector<char> > haps2=haps;
//    for(unsigned int i=0; i<haps.size();i++) {
//        cout << "Char converted partial hap is: ";
//        for(unsigned int j=0; j<haps[i].size();j++) {
//            cout << haps[i][j];
//        }
//        cout << "\n";
//    }
    
    //Matching overlap joins
    for (int i=0;i<200;i++) { //Change to while loop and include a stopping criteria?
        haps2=haps;
        OverlapStepShared1(haps); //Remove haplotypes contained within each other
        OverlapStepShared2(haps);
        if (haps==haps2) {
            break;
        }
    }
    
    //Non-matching overlap joins
    for (int i=0;i<20;i++) {
        haps2=haps;
        OverlapStepNoShare(haps);
        OverlapStepShared1(haps);
        OverlapStepShared2(haps);
        if (haps==haps2) {
            break;
        }
    }
    
    
	cout << "C\n";
    
    //Non-overlap joins
    vector<AnalyserPH::par> sf;
    GetStartFinish(sf,haps);
    vector< vector<char> > new_haps;
    BuildPaths(sf,haps,new_haps);
    haps=new_haps;
    DuplicateStep(haps);
    
	cout << "D\n";
    
    //Converting from vector<char> to Sequence
  //  cout << "Full haps:\n";
    for(unsigned int i=0; i<haps.size(); i++) {
        Sequence s;
        for(unsigned int j=0; j<haps[i].size(); j++) {
            s.addBase(haps[i][j]);
        }
        //s.print();
        fhs->push_back(s);
    }
    
    
}

//Method for updating haps following a step (e.g. overlap step, combination step, etc.)
void AnalyserPH::ResetHaps (vector<int> incl, vector< vector<char> >& haps) {
    vector< vector<char> > new_haps;
    for (unsigned int i=0;i<haps.size();i++) {
        if (incl[i]==1) {
            new_haps.push_back(haps[i]);
        }
    }
    haps=new_haps;
}


//Remove duplicate haplotypes
void AnalyserPH::DuplicateStep (vector< vector<char> >& haps) {
    vector< vector<char> > new_haps;
    vector<int> uniq;
    for (unsigned int i=0;i<haps.size();i++) {
        uniq.push_back(1);
    }
    //cout << "Haplotypes " << haps.size() << "\n";
    for (unsigned int i=0;i<haps.size();i++) { //Simply loop over all combinations i,j>i and check if any of them are duplicates
        for (unsigned int j=i;j<haps.size();j++) {
            if (i!=j&&haps[i]==haps[j]) {
                uniq[j]=0;
            }
        }
    }
    for (unsigned int i=0;i<haps.size();i++) {
        if (uniq[i]==1) {
            new_haps.push_back(haps[i]);
        }
    }
    // cout << "Unique haplotypes " << new_haps.size() << "\n";
    haps=new_haps;
}

//Remove haplotypes contained within each other
void AnalyserPH::OverlapStepShared1 (vector< vector<char> >& haps) {
    //One contained within another - try doing these first before other overlaps...
    vector<int> incl (haps.size(),1); //Haplotypes that should be kept have a 1, those to be removed have a 0
    int size=haps.size();
    for (int i=0;i<size;i++) {
        for (int j=0;j<size;j++) {
            if (i!=j) {
                int match=1;   //Match at shared positions
                int n_match=0;  //Number of shared positions
                for	(unsigned int k=0;k<haps[i].size();k++) {
                    //Check whether j is contained in i.
                    if (haps[j][k]!='-') { //Here we are only considering the alleles k in hap j that isn't '-'
                        if (haps[i][k]!='-') {
                            n_match++;
                        }
                        if (haps[j][k]!=haps[i][k]) { //If just a single one of the non-'-' alleles in hap j is not identical to hap i, then no match
                            match=0;
                        }
                    }
                }
                if (match==1&&n_match>0) {
                    //Remove haps[j]
                    incl[j]=0;
                }
            }
        }
    }
    ResetHaps(incl,haps); //Reset haps, i.e. only keep haps with a 1 in incl
    DuplicateStep(haps); //Remote any possible duplicate haplotypes
}




//Remove shared overlap like AACC-- --CCTT
void AnalyserPH::OverlapStepShared2 (vector< vector<char> >& haps) {
    vector<int> incl (haps.size(),1);
    int size=haps.size();
    for (int i=0;i<size;i++) {
        for (int j=i+1;j<size;j++) {
            int match=1;   //Match at shared positions
            int n_match=0;  //Number of shared positions
            for	(unsigned int k=0;k<haps[i].size();k++) {
                //Check whether i and j have the same nucleotides wherever they overlap.
                if (haps[j][k]!='-'&&haps[i][k]!='-') { //Number of overlaps
                    n_match++;
                    if (haps[j][k]!=haps[i][k]) {
                        match=0;
                    }
                }
            }
            if (match==1&&n_match>0) {
//                for	(int k=0;k<haps[i].size();k++) {
//                    cout << haps[i][k];
//                }
//                cout << " overlaps ";
//                for	(int k=0;k<haps[j].size();k++) {
//                    cout << haps[j][k];
//                }
//                cout << "\n";
                vector<char> newhap;
                for	(unsigned int k=0;k<haps[i].size();k++) {
                    if (haps[j][k]!='-') {
                        newhap.push_back(haps[j][k]);
                    } else if (haps[i][k]!='-') {
                        newhap.push_back(haps[i][k]);
                    } else {
                        newhap.push_back('-');
                    }
                }
                //Add newhap, remove the other two
                haps.push_back(newhap);
                incl.push_back(1);
                incl[i]=0;
                incl[j]=0;
            }
        }
    }
    
    //Update haps and remove any possible duplicates
    ResetHaps(incl,haps);
    DuplicateStep(haps);
}


//Non-matching overlaps
void AnalyserPH::OverlapStepNoShare (vector< vector<char> >& haps) {
    
    vector<int> incl (haps.size(),1);
    int size=haps.size();
    for (int i=0;i<size;i++) {
        for (int j=i+1;j<size;j++) {
            int match=1;   //Match at shared positions
            int n_match=0;  //Number of shared positions
            for	(unsigned int k=0;k<haps[i].size();k++) {
                //Check whether i and j have the same nucleotides wherever they overlap.  N.B. Must have at least one locus not in common
                if (haps[j][k]!='-'&&haps[i][k]!='-') { //Number of overlaps
                    n_match++;
                    if (haps[j][k]!=haps[i][k]) {
                        match=0;
                    }
                }
            }
            int n_over=0;
            if (match==0&&n_match<(int) haps[i].size()&&n_match>0) { //some shared positions, but not a complete match
                for	(unsigned int k=0;k<haps[i].size();k++) {
                    if (haps[j][k]!='-'&&haps[i][k]=='-') {
                        n_over=1; //There is an overlap of some kind
                    }
                    if (haps[i][k]!='-'&&haps[j][k]=='-') {
                        n_over=1; //There is an overlap of some kind
                    }
                }
            }
            
            //Overlap present
            if (n_over==1) {
                
                vector<char> nh_start;
                //Make beginning
                int ov=0;
                for	(unsigned int k=0;k<haps[i].size();k++) {
                    if (haps[i][k]!='-'&&haps[j][k]!='-') {
                        ov=1;
                    }
                    if (ov==0) {
                        if (haps[j][k]!='-') {
                            nh_start.push_back(haps[j][k]);
                        } else if (haps[i][k]!='-') {
                            nh_start.push_back(haps[i][k]);
                        } else {
                            nh_start.push_back('-');
                        }
                    } else {
                        nh_start.push_back('-');
                    }
                }
                
                //Make end
                vector<char> nh_end;
                ov=0;
                for	(unsigned int k=0;k<haps[i].size();k++) {
                    if (haps[i][k]!='-'&&haps[j][k]!='-') {
                        ov=1;
                    }
                    if (ov==1) {
                        if (haps[j][k]!='-'&&haps[i][k]=='-') {
                            nh_end.push_back(haps[j][k]);
                        } else if (haps[i][k]!='-'&&haps[j][k]=='-') {
                            nh_end.push_back(haps[i][k]);
                        } else {
                            nh_end.push_back('-');
                        }
                    } else {
                        nh_end.push_back('-');
                    }
                }
                
                //Make c1 and c2
                vector<char> nh_c1;
                vector<char> nh_c2;
                for	(unsigned int k=0;k<haps[i].size();k++) {
                    if (nh_start[k]=='-'&&nh_end[k]=='-') {
                        nh_c1.push_back(haps[i][k]);
                        nh_c2.push_back(haps[j][k]);
                    } else {
                        nh_c1.push_back('-');
                        nh_c2.push_back('-');
                    }
                }
                
                //Construct start + c1 + end
                vector<char> newhap1;
                vector<char> newhap2;
                for (unsigned int k=0;k<haps[i].size();k++){
                    if (nh_start[k]!='-') {
                        newhap1.push_back(nh_start[k]);
                        newhap2.push_back(nh_start[k]);
                    } else if (nh_c1[k]!='-') {
                        newhap1.push_back(nh_c1[k]);
                        newhap2.push_back(nh_c2[k]);
                    } else if (nh_end[k]!='-') {
                        newhap1.push_back(nh_end[k]);
                        newhap2.push_back(nh_end[k]);
                    } else {
                        newhap1.push_back('-');
                        newhap2.push_back('-');
                    }
                }
                
                haps.push_back(newhap1);
                haps.push_back(newhap2);
                incl.push_back(1);
                incl.push_back(1);
                //cout << newhap1.size() << " " << newhap2.size() << "\n";
                incl[i]=0;
                incl[j]=0;
            }
        }
    }
    ResetHaps(incl,haps);
    DuplicateStep(haps);
}

void AnalyserPH::GetStartFinish(vector<AnalyserPH::par>& sf, vector< vector<char> >& haps) {
    for (unsigned int i=0;i<haps.size();i++) {
        AnalyserPH::par p;
        p.i1=-1;
        p.i2=-1;
        for (unsigned int j=0;j<haps[i].size();j++) {
            //cout << "i " << i << " j " << j << " " << haps[i][j] << "\n";
            if (haps[i][j]!='-'&&p.i1<0) {
                p.i1=j;
                //	cout << "Found start\n";
            }
            if (haps[i][j]=='-'&&p.i1>=0&&p.i2<0) {
                p.i2=j-1;
                //	cout << "Found end\n";
            }
        }
        if (p.i2<0) {
            p.i2=haps[i].size();
        }
        sf.push_back(p);
    }
}


void AnalyserPH::AddPath (int prev, int start, vector<char> v, vector<AnalyserPH::par>& sf, vector< vector<char> >& haps, vector< vector<char> >& new_haps) {
    //Add last haplotype to vector
    vector<char> v2=v;
    //cout << "Previous is " << prev << " " << sf[prev].i1 << " " << sf[prev].i2 << "\n";
    for (int j=sf[prev].i1;j<=sf[prev].i2;j++) {
        v2.push_back(haps[prev][j]);
    }
    //cout << "Vector v is: ";
    //for (int i=0;i<v2.size();i++) {
    //	cout << v2[i];
    //}
    //cout << "\n";
    
    //Check if this is the end of the haplotype and if so push to new_haps
    if (sf[prev].i2>=(int) haps[prev].size()) {
        //	cout << "Doing push_back\n";
        vector<char> v3;
        for (unsigned int i=0;i<haps[prev].size();i++) {//Trim length
            v3.push_back(v2[i]);
        }
        new_haps.push_back(v3);
    } else {
        //If not search for new haplotype and recall
        for (unsigned int i=0;i<sf.size();i++) {
            if (sf[i].i1==start) {
                AddPath(i,sf[i].i2+1,v2,sf,haps,new_haps);
            }
        }
    }
}

void AnalyserPH::BuildPaths(vector<AnalyserPH::par>& sf, vector< vector<char> > haps, vector< vector<char> >& new_haps) {
    vector<char> v;
    for (unsigned int i=0;i<sf.size();i++) {
        if (sf[i].i1==0) {
            //cout << "Detected " << i << "\n";
            AddPath(i,sf[i].i2+1,v,sf,haps,new_haps);
        }
    }
}




double AnalyserPH::computeLSelection(AnaParam *ap, int Nt, vector<double> &qBfhFreqs, vector<double> &fhHapFit) {

	double Ltot = 0;
   
	//Compute likelihood
	for(int i=0; i<data->numOfPHSets(); i++) { //Loop over (independent) partial haplotype sets within a simulation
        
		vector<int> obsA = data->getNa(i); //Bad choice of function name. Should be called getObsA(i);
		vector<int> obsB = data->getNb(i);	

		vector<vector<int> > contribsPHset; //This one could be created earlier to save reproducing every time
		for(int j=0; j<data->getNumOfPartialHaps(i); j++) {

			contribsPHset.push_back(*data->getContribs(i,j));
		}
        
		Ltot += ::computeLSelection(obsB, obsA, Nt, ap, qBfhFreqs, fhHapFit, &contribsPHset, false); //Defined in misc
	}

    
	return Ltot;

}

double AnalyserPH::computeLSelTVar(AnaParam *ap, int Nt, vector<double> &qBfhFreqs, gsl_matrix* qBvar, vector<double> &fhHapFitT) {

	double Ltot = 0;
   
	//Compute likelihood
	for(int i=0; i<data->numOfPHSets(); i++) { //Loop over (independent) partial haplotype sets within a simulation
        
		vector<int> obsA = data->getNa(i); //Bad choice of function name. Should be called getObsA(i);

		vector<vector<int> > contribsPHset;
		for(int j=0; j<data->getNumOfPartialHaps(i); j++) {

			contribsPHset.push_back(*data->getContribs(i,j));
		}

		Ltot += ::computeLSelTVar(obsA, Nt, ap, qBfhFreqs, qBvar, fhHapFitT, &contribsPHset, false); //Defined in misc
	}

	return Ltot;
}

double AnalyserPH::computeLSelTSelGVar(AnaParam *ap, int Nt, vector<double> &qBfhFreqs, gsl_matrix* qBvar, vector<double> &fhHapFitT, vector<double> &fhHapFitG) {

	double Ltot = 0;
   
	//Compute likelihood
	for(int i=0; i<data->numOfPHSets(); i++) { //Loop over (independent) partial haplotype sets within a simulation
        
		vector<int> obsA = data->getNa(i); //Bad choice of function name. Should be called getObsA(i);

		vector<vector<int> > contribsPHset;
		for(int j=0; j<data->getNumOfPartialHaps(i); j++) {

			contribsPHset.push_back(*data->getContribs(i,j));
		}

		Ltot += ::computeLSelTSelGVar(obsA, Nt, ap, qBfhFreqs, qBvar, fhHapFitT, fhHapFitG, &contribsPHset, false); //Defined in misc

	}

	return Ltot;
}

double AnalyserPH::computeLNeutral(AnaParam* ap, int Nt, vector<double> &qBfhFreqs) {

	double Ltot = 0;
   
	//Compute likelihood
	for(int i=0; i<data->numOfPHSets(); i++) { //Loop over (independent) partial haplotype sets within a simulation
        
		vector<int> obsA = data->getNa(i); //Bad choice of function name. Should be called getObsA(i);
		vector<int> obsB = data->getNb(i);	

		vector<vector<int> > contribsPHset;
		for(int j=0; j<data->getNumOfPartialHaps(i); j++) {

			contribsPHset.push_back(*data->getContribs(i,j));
		}

		Ltot += ::computeLNeutral(obsB, obsA, Nt, ap, qBfhFreqs, &contribsPHset, false); //Defined in misc
	}
    
	return Ltot;

}


double AnalyserPH::computeLNeutralVar(AnaParam* ap, int Nt, vector<double> &qBfhFreqs, gsl_matrix* qBvar) {

	double Ltot = 0;
   
	//Compute likelihood
	for(int i=0; i<data->numOfPHSets(); i++) { //Loop over (independent) partial haplotype sets within a simulation
        
		vector<int> obsA = data->getNa(i); //Bad choice of function name. Should be called getObsA(i);

		vector<vector<int> > contribsPHset;
		for(int j=0; j<data->getNumOfPartialHaps(i); j++) {

			contribsPHset.push_back(*data->getContribs(i,j));
		}

		Ltot += ::computeLNeutralVar(obsA, Nt, ap, qBfhFreqs, qBvar, &contribsPHset, false); //Defined in misc
//		double Ltemp;
//		if(i==0) {
//			cout << "obsA: ";  printIntVector(obsA);
//			Ltemp = ::computeLNeutralVar(obsA, Nt, ap, qBfhFreqs, qBvar, &contribsPHset, true); //Defined in misc
//		} else {
//
//			Ltemp = ::computeLNeutralVar(obsA, Nt, ap, qBfhFreqs, qBvar, &contribsPHset, false); //Defined in misc
//
//		}
//		cout << "Likelihood for SOPH " << i << ": " << Ltemp << "\n";
//		Ltot += Ltemp;
	
	}


	return Ltot;

}


double AnalyserPH::computeLNeutralVarAdvanced(AnaParam* ap, int Nt, vector<double> &qBfhFreqs, gsl_matrix* qBvar) {

	double Ltot = 0;
	int deltaDays = data->getDeltaDays();
   
	//Compute likelihood
	for(int i=0; i<data->numOfPHSets(); i++) { //Loop over (independent) partial haplotype sets within a simulation
        
		vector<int> obsA = data->getNa(i); //Bad choice of function name. Should be called getObsA(i);

		vector<vector<int> > contribsPHset;
		for(int j=0; j<data->getNumOfPartialHaps(i); j++) {

			contribsPHset.push_back(*data->getContribs(i,j));
		}

		Ltot += ::computeLNeutralVarAdvanced(obsA, deltaDays, Nt, ap, qBfhFreqs, qBvar, &contribsPHset, false); //Defined in misc
	
	}


	return Ltot;

}

double AnalyserPH::computeLNeutralSelGVar(AnaParam* ap, int Nt, vector<double> &qBfhFreqs, gsl_matrix* qBvar, vector<double> &hapFitG) {

	double Ltot = 0;
   
	//Compute likelihood
	for(int i=0; i<data->numOfPHSets(); i++) { //Loop over (independent) partial haplotype sets within a simulation
        
		vector<int> obsA = data->getNa(i); //Bad choice of function name. Should be called getObsA(i);

		vector<vector<int> > contribsPHset;
		for(int j=0; j<data->getNumOfPartialHaps(i); j++) {

			contribsPHset.push_back(*data->getContribs(i,j));
		}

		
		
		Ltot += ::computeLNeutralSelGVar(obsA, Nt, ap, qBfhFreqs, qBvar, hapFitG, &contribsPHset, false); //Defined in misc
//		double Ltemp;
//		if(i==21) {
//			cout << "obsA: ";  printIntVector(obsA);
//			Ltemp = ::computeLNeutralSelGVar(obsA, Nt, ap, qBfhFreqs, qBvar, hapFitG, &contribsPHset, true); //Defined in misc
//		} else {
//
//			Ltemp = ::computeLNeutralSelGVar(obsA, Nt, ap, qBfhFreqs, qBvar, hapFitG, &contribsPHset, false); //Defined in misc
//		}
//		cout << "Likelihood for soph " << i << ": " << Ltemp << "\n";
//		Ltot += Ltemp;
	}

	return Ltot;
}

vector<double> AnalyserPH::findOptimalFreqFHBefore(int numFHs, gsl_rng *r, int C, vector<double>* logFactStore) {
    
    vector<double> oldFreqs;
    do { //Ensure that the starting point doesn't have zero values
        oldFreqs.clear();
        for(int i=0; i<numFHs; i++) { oldFreqs.push_back(gsl_rng_uniform_pos(r)); }
        rescale(oldFreqs);
    } while(nonZeroDoubleVector(oldFreqs) == false);
   
	vector<double> newFreqs;
	double oldL = - std::numeric_limits<double>::max();
	double newL = oldL;
	double delta = 1;
	int accepts = 0;
	int trys = 0;
	int nUpdate = 50; //Frequency for updating acceptance rate
	
	int k=1; //Counter
	while(delta > 0.0001) {
		k++;
		if(k%nUpdate==0) {

			double acceptanceRate=(accepts+0.)/(trys+0.);
			//cout << "Acceptance rate: " << acceptanceRate << "\tdelta: " << delta << "\n";


			double tempDelta=delta*(0.95+acceptanceRate);
			if(tempDelta < 10) { // If delta>10 computations get out of ahnd
				delta = tempDelta;
			}

			//Update value for nUpdate, i.e. the number of rounds before the acceptance value is computed
			//For high acceptance values we want to spend more time with the current delta
			//For low acceptance value, we wish to progress towards smaller deltas faster
			if(acceptanceRate > 0.05) {
				nUpdate += 10;  //Increase nUpdate
			} else {
				if(nUpdate >= 20) {
					nUpdate -= 10; //Decrease nUpdate only if resulting value is >= 10
				}
			}

			//Reset
			accepts=0;
			trys=0;
		}
	
		do{
			addRandom(oldFreqs,newFreqs,delta,r); //Create a new vector newFreqs from adding a small delta*rndm(double) to a random index
		} while(nonZeroDoubleVector(newFreqs)==false); //Ensure that the new frequencies do not have zero elements

		newL = computeLb(newFreqs, r, C, logFactStore);

		if(newL > oldL) {
			oldFreqs = newFreqs;
			oldL = newL;
			accepts++;
        }	
		
		trys++;
	}
  
    return oldFreqs;
}

vector<double> AnalyserPH::findOptimalFreqFHAfter(int numFHs, gsl_rng *r, int C, vector<double>* logFactStore) {
    
    vector<double> oldFreqs;
    do { //Ensure that the starting point doesn't have zero values
        oldFreqs.clear();
        for(int i=0; i<numFHs; i++) { oldFreqs.push_back(gsl_rng_uniform_pos(r)); }
        rescale(oldFreqs);
    } while(nonZeroDoubleVector(oldFreqs) == false);
   
	vector<double> newFreqs;
	double oldL = - std::numeric_limits<double>::max();
	double newL = oldL;
	double delta = 1;
	int accepts = 0;
	int trys = 0;
	int nUpdate = 50; //Frequency for updating acceptance rate
	
	int k=1; //Counter
	while(delta > 0.0001) {
		k++;
		if(k%nUpdate==0) {

			double acceptanceRate=(accepts+0.)/(trys+0.);
			//cout << "Acceptance rate: " << acceptanceRate << "\tdelta: " << delta << "\n";


			double tempDelta=delta*(0.95+acceptanceRate);
			if(tempDelta < 10) { // If delta>10 computations get out of ahnd
				delta = tempDelta;
			}

			//Update value for nUpdate, i.e. the number of rounds before the acceptance value is computed
			//For high acceptance values we want to spend more time with the current delta
			//For low acceptance value, we wish to progress towards smaller deltas faster
			if(acceptanceRate > 0.05) {
				nUpdate += 10;  //Increase nUpdate
			} else {
				if(nUpdate >= 20) {
					nUpdate -= 10; //Decrease nUpdate only if resulting value is >= 10
				}
			}

			//Reset
			accepts=0;
			trys=0;
		}
	
		do{
			addRandom(oldFreqs,newFreqs,delta,r); //Create a new vector newFreqs from adding a small delta*rndm(double) to a random index
		} while(nonZeroDoubleVector(newFreqs)==false); //Ensure that the new frequencies do not have zero elements

		//newL = computeLb(newFreqs, r, C, logFactStore);
		newL = computeLa(newFreqs, r, C, logFactStore);

		if(newL > oldL) {
			oldFreqs = newFreqs;
			oldL = newL;
			accepts++;
        }	
		
		trys++;
	}
  
    return oldFreqs;
}

vector<double> AnalyserPH::findOptimalFreqFHBefore(int numFHs, double *L, gsl_rng *r, int C, vector<double>* logFactStore) {
    
    vector<double> oldFreqs;
    do { //Ensure that the starting point doesn't have zero values
        oldFreqs.clear();
        for(int i=0; i<numFHs; i++) { oldFreqs.push_back(gsl_rng_uniform_pos(r)); }
        rescale(oldFreqs);
    } while(nonZeroDoubleVector(oldFreqs) == false);
   
	vector<double> newFreqs;
	double oldL = - std::numeric_limits<double>::max();
	double newL = oldL;
	double delta = 1;
	int accepts = 0;
	int trys = 0;
	int nUpdate = 50; //Frequency for updating acceptance rate
	
	int k=1; //Counter
	while(delta > 0.001) {
		k++;
		if(k%nUpdate==0) {

			double acceptanceRate=(accepts+0.)/(trys+0.);
			//cout << "Acceptance rate: " << acceptanceRate << "\tdelta: " << delta << "\n";


			double tempDelta=delta*(0.95+acceptanceRate);
			if(tempDelta < 10) { // If delta>10 computations get out of ahnd
				delta = tempDelta;
			}

			//Update value for nUpdate, i.e. the number of rounds before the acceptance value is computed
			//For high acceptance values we want to spend more time with the current delta
			//For low acceptance value, we wish to progress towards smaller deltas faster
			if(acceptanceRate > 0.05) {
				nUpdate += 10;  //Increase nUpdate
			} else {
				if(nUpdate >= 20) {
					nUpdate -= 10; //Decrease nUpdate only if resulting value is >= 10
				}
			}

			//Reset
			accepts=0;
			trys=0;
		}
	
		do{
			addRandom(oldFreqs,newFreqs,delta,r); //Create a new vector newFreqs from adding a small delta*rndm(double) to a random index
		} while(nonZeroDoubleVector(newFreqs)==false); //Ensure that the new frequencies do not have zero elements

		newL = computeLb(newFreqs, r, C, logFactStore);

		if(newL > oldL) {
			oldFreqs = newFreqs;
			oldL = newL;
			*L = newL; //This is a way of exporting the likelihood
			accepts++;
        }	
		
		trys++;
	}
  
    return oldFreqs;
}

vector<double> AnalyserPH::findOptimalFreqFHBeforeFromMMS(vector<vector<Sequence> > &MMSfhs, vector<Sequence> *fhs, gsl_rng *r, int C, vector<double> *logFactStore) {


	vector<vector<double> > bestFreqsMMS;
	vector<vector<Sequence> > bestSeqsMMS;
	vector<double> bestBICsMMS;

	for(unsigned int i = 0; i<MMSfhs.size(); i++) { //Loop over MMS haplotype sets

		vector<Sequence> oldSeqs = MMSfhs[i];
		vector<Sequence> newSeqs;
		double L;
		data->computeContribs(oldSeqs); //Update the contributions to reflect the the original MMS
		vector<double> oldFreqs = findOptimalFreqFHBefore((int) oldSeqs.size(), &L, r, C, logFactStore); //Also updates L
		vector<double> newFreqs;
		double oldBIC = -2*L; //Starting BIC based on original MSS
		double newBIC;
		vector<Sequence> bestSeqs;

		//Find the set of seqs that is not int oldSeqs
		vector<Sequence> remainingSeqs = findRemainingSeqs(oldSeqs, *fhs);
		bool improvementFound = true;

		//Try adding haplotypes to best set
		//int counter = 0;
		while(improvementFound == true && remainingSeqs.size() > 0) {
			//if(counter == 1) { break; }

			improvementFound = false;
			cout << "oldBIC at start of round: " << oldBIC << "\n";
			for(unsigned int j=0; j<remainingSeqs.size(); j++) { //Loop over remaining haplotypes and try and include them one at a time
				
				cout << "Length of remainingSeqs: " << remainingSeqs.size() << "\n";
				newSeqs = oldSeqs;
				newSeqs.push_back(remainingSeqs[j]);
				data->computeContribs(newSeqs); //Update the contributions to reflect the new combination of full haplotypes
			
				newFreqs = findOptimalFreqFHBefore((int)newSeqs.size(), &L, r, C, logFactStore);

				int numOfParams = (int) (newSeqs.size() - MMSfhs[i].size());
				newBIC = -2*L + numOfParams*log(data->getNbTot());

				cout << "oldBic: " << oldBIC << "\n";
				cout << "newBic: " << newBIC << "\n";
				if(newBIC < oldBIC) {
					improvementFound = true;
					oldBIC = newBIC;
					bestSeqs = newSeqs;
					oldFreqs = newFreqs;
					
					cout << "Improvement found: \n";
					cout << "Length: " << newSeqs.size() << "\n";
					cout << "BIC: " << newBIC << "\n";
					for(unsigned int k=0; k<newSeqs.size(); k++) {
						newSeqs[k].print();
					}
				}
			}

			if(improvementFound == true) {
				oldSeqs = bestSeqs;
				remainingSeqs = findRemainingSeqs(oldSeqs, *fhs);
			}
			//counter ++;
		}

		//Try removing haplotypes from best set
		cout << "Trying to remove haplotypes from MMS:\n";
		improvementFound = true;
		while(improvementFound == true && oldSeqs.size() > 1) {

			improvementFound = false;
			//cout << "oldBIC at start of round: " << oldBIC << "\n";
			for(unsigned int j=0; j<oldSeqs.size(); j++) { //Loop over haplotypes and try and remove one at a time
				
				newSeqs = oldSeqs;
				newSeqs.erase(newSeqs.begin()+j,newSeqs.begin()+j+1); //Erase j'th haplotype

				//Need to check that newSeqs includes the correct partial haplotypes
				//If not, can't remove the current haplotype
				if(data->decomposeAndCheck(newSeqs) == true) {

					data->computeContribs(newSeqs); //Update the contributions to reflect the new combination of full haplotypes
				
					newFreqs = findOptimalFreqFHBefore((int)newSeqs.size(), &L, r, C, logFactStore);

					int numOfParams = (int) (newSeqs.size() - MMSfhs[i].size()); //NB: Can be negative! It is all relative
					newBIC = -2*L + numOfParams*log(data->getNbTot());
	
					//cout << "oldBic: " << oldBIC << "\n";
					//cout << "newBic: " << newBIC << "\n";
					if(newBIC < oldBIC) {
						improvementFound = true;
						oldBIC = newBIC;
						bestSeqs = newSeqs;
						oldFreqs = newFreqs;
						
						//cout << "Improvement found: \n";
						//cout << "Length: " << newSeqs.size() << "\n";
						//cout << "BIC: " << newBIC << "\n";
						//for(unsigned int k=0; k<newSeqs.size(); k++) {
						//	newSeqs[k].print();
						//}
					}
				}
			}

			if(improvementFound == true) {
				oldSeqs = bestSeqs;
			}
		}

		//Try to swap haplotypes
		//cout << "Trying to swap haplotypes from MMS and remaindingr:\n";
		//improvementFound = true;
		//remainingSeqs = findRemainingSeqs(oldSeqs, *fhs);
		//while(improvementFound == true && oldSeqs.size() > 0 && remainingSeqs.size()>0) {

		//	improvementFound = false;
		//	cout << "oldBIC at start of round: " << oldBIC << "\n";
		//	for(unsigned int j=0; j<oldSeqs.size(); j++) { //Loop over haplotypes and try and swap with a haplotype from remainingSeqs
	
		//		for(unsigned int k=0; k<remainingSeqs.size(); k++) {
		//			newSeqs = oldSeqs;
		//			newSeqs[j]=remainingSeqs[k]; //Swap jth haplotype of oldSeq with kth haplotype from remainingSeqs

		//			//Need to check that newSeqs includes the correct partial haplotypes
		//			//If not, can't remove the current haplotype
		//			if(data->decomposeAndCheck(newSeqs, simIndex) == true) {
	
		//				data->computeContribs(newSeqs, simIndex); //Update the contributions to reflect the new combination of full haplotypes
					
		//				newFreqs = findOptimalFreqFHBefore((int)newSeqs.size(), &L, r, simIndex, C, logFactStore);

		//				int numOfParams = (int) (newSeqs.size() - MMSfhs[i].size()); //NB: Can be negative! It is all relative
		//				newBIC = -2*L + numOfParams*log(data->getNbTot(simIndex));
	
		//				//cout << "oldBic: " << oldBIC << "\n";
		//				//cout << "newBic: " << newBIC << "\n";
		//				if(newBIC < oldBIC -0.1) { //Swapping can lead to changes with almost indistinguishable BICs, so be more conservative
		//					improvementFound = true;
		//					oldBIC = newBIC;
		//					bestSeqs = newSeqs;
		//					oldFreqs = newFreqs;
						
		//					cout << "Improvement found: \n";
		//					cout << "Length: " << newSeqs.size() << "\n";
		//					cout << "BIC: " << newBIC << "\n";
		//					for(unsigned int k=0; k<newSeqs.size(); k++) {
		//						newSeqs[k].print();
		//					}
		//				}
		//			}
		//		}
		//	}

		//	if(improvementFound == true) {
		//		oldSeqs = bestSeqs;
		//		remainingSeqs = findRemainingSeqs(oldSeqs, *fhs);
		//	}
		//}
		bestFreqsMMS.push_back(oldFreqs);
		bestSeqsMMS.push_back(oldSeqs);
		bestBICsMMS.push_back(oldBIC);
	}


	vector<double> bestFreqs;
	double bestBIC = std::numeric_limits<double>::max();
	for(unsigned int i=0; i<bestBICsMMS.size(); i++) {
		if(bestBICsMMS[i] < bestBIC) {
			bestBIC = bestBICsMMS[i];
			bestFreqs = bestFreqsMMS[i];
			*fhs = bestSeqsMMS[i]; //Return best seqs through inputted full haplotypes
		}
	}
	return bestFreqs;

}

void AnalyserPH::optimiseqBmeanFromPre(int numFHs, gsl_rng *r, int C, vector<double>* qBmean, double* pL) {

	int dim = numFHs;
	qBmean->clear();

	//Initialise mean to random values, then renormalise
	for(int i=0; i<dim; i++) {

		qBmean->push_back(gsl_rng_uniform_pos(r));
	}
	rescale(*qBmean);

  
	//Define "new" parameters
	vector<double> qBmeanNew = *qBmean;

	
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
		addRandom(*qBmean,qBmeanNew,delta,r); //Create a new vector newFreqs from adding a small delta*rndm(double) to a random index

		double newL = computeLmeanPreOnly(&qBmeanNew, r, C, false); 

		if(newL > bestL) {

			(*qBmean) = qBmeanNew;

			bestL = newL;
			accepts++;
		
		} else {
			//Revert back any changes
			qBmeanNew = *qBmean;
		}	
		
		trys++;
	}
 
	//Set pointer to L to bestL, so can be printed outside this method
	*pL = bestL;
} 





void AnalyserPH::optimiseqBmeanAndVarFromPre(int numFHs, gsl_rng *r, int C, vector<double>* qBmean, gsl_matrix* qBvar, double* pL) {

	int dim = numFHs;
	qBmean->clear();

	//Initialise mean and var
	double value = 0.1; //1e-10; //0.1; //For simplicity, initialise all values to 1/k
	for(int i=0; i<dim; i++) { 

		qBmean->push_back(gsl_rng_uniform_pos(r)); //Set mean to random values, then renormalise
		for(int j=0; j<dim; j++) {

			//Matrix diagonals set to 1/k - everything else set to zero
			if(j==i) {
				gsl_matrix_set(qBvar, i, j, value);
				//gsl_matrix_set(qBvar, i, j, gsl_rng_uniform_pos(r)*100);
			} else {
				gsl_matrix_set(qBvar, i, j, 0.0);
				//gsl_matrix_set(qBvar, i, j, gsl_rng_uniform_pos(r)*100);
			}
		}
	}
	rescale(*qBmean);

  
	//Define "new" parameters
	vector<double> qBmeanNew = *qBmean;
	gsl_matrix* qBvarNew = gsl_matrix_alloc(dim,dim);
	gsl_matrix_memcpy(qBvarNew, qBvar);

	
	double bestL = - std::numeric_limits<double>::max();
	double delta = 1;
	int accepts = 0;
	int trys = 0;
	int nUpdate = 1000; //Frequency for updating acceptance rate

	int k=1; //Counter
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

			//if(gsl_rng_uniform(r) < 0.005 && delta < 0.1) { delta *= 100; } //Random jump to higher delta

			//Reset
			accepts=0;
			trys=0;
		}
	
		if(gsl_rng_uniform(r) < 0.5) { //Update mean
	
			//addRandomMin(*qBmean,qBmeanNew,delta,r,1e-7); //Create a new vector newFreqs from adding a small delta*rndm(double) to a random index
			//addRandomMin(*qBmean,qBmeanNew,delta,r,1e-2); //Create a new vector newFreqs from adding a small delta*rndm(double) to a random index
			addRandom(*qBmean,qBmeanNew,delta,r); //Create a new vector newFreqs from adding a small delta*rndm(double) to a random index


		} else { //Update variance


			bool positiveDefinite = false;

			while(positiveDefinite == false) {
				int index1 = gsl_rng_uniform_int(r,dim);
				int index2 = index1; // = gsl_rng_uniform_int(r,dim);
				double oldValue = gsl_matrix_get(qBvar, index1, index2);
				double newValue;
				do {
					newValue = oldValue + delta*(gsl_rng_uniform_pos(r) - 0.5);
				} while(newValue < 10e-11); 

				if(index1 == index2) {
						gsl_matrix_set(qBvarNew, index1, index1, newValue);
				} else {
					gsl_matrix_set(qBvarNew, index1, index2, newValue);
					gsl_matrix_set(qBvarNew, index2, index1, newValue);
				}
		
				//Check if positive-definite
				gsl_matrix* qBvarCholesky = gsl_matrix_alloc(dim,dim);
				gsl_matrix_memcpy(qBvarCholesky, qBvarNew); //Copy into qBvarCholesky as matrix is modified by action below

				//The function below is deprecated by gsl. New function is gsl_linalg_cholesky_decomp1(), but this version of gsl doesn't support that.
				//Hence we use the deprecated version
				gsl_set_error_handler_off();
				int status = gsl_linalg_cholesky_decomp(qBvarCholesky);
				if(status == 0) { //Success indicates that matrix is symmetric positive definite

					positiveDefinite = true;
				} else {
					
					gsl_matrix_memcpy(qBvarNew, qBvar); //Not positive definite, so revert changes
				}
				
				gsl_matrix_free(qBvarCholesky);
			}
		}

		double newL = computeLmeanAndVarPreOnly(&qBmeanNew, qBvarNew, r, C); 
		
		if(newL > bestL) {

			//cout << "Better mean and var found.\n";

			(*qBmean) = qBmeanNew;
			gsl_matrix_memcpy(qBvar, qBvarNew);

	//		cout << "Priting qBmean: \n";
	//		printDoubleVector(*qBmean);
	//		cout << "Printing qBvar: \n";
	//		printMatrixMathematica(qBvar);
	
//			gsl_permutation *p = gsl_permutation_calloc(dim);		
//			int signum;
//			gsl_matrix * tmp_ptr = gsl_matrix_calloc(dim,dim);
//			gsl_matrix_memcpy(tmp_ptr,qBvar); //Copy elements of sigma into tmp_ptr
//  
//			//Get LU decomposition
//			gsl_linalg_LU_decomp(tmp_ptr,p,&signum);
//  
//			//Get determinant
//			double det = gsl_linalg_LU_det(tmp_ptr, signum);
//			double lndet = gsl_linalg_LU_lndet(tmp_ptr);
//	
//			cout << "Better mean and var found. Printing determinant: " << det << "\n";
//			cout << "Better mean and var found. Printing ln(determinant): " << lndet << "\n";
//			gsl_matrix_free(tmp_ptr);
//			gsl_permutation_free(p);
 
			bestL = newL;
			accepts++;
		
		} else {
			//Revert back any changes
			qBmeanNew = *qBmean;
			gsl_matrix_memcpy(qBvarNew, qBvar);
		}	
		
		trys++;
	}
 
	gsl_matrix_free(qBvarNew);

	//Set pointer to L to bestL, so can be printed outside this method
	*pL = bestL;
} 


void AnalyserPH::optimiseqBmeanAndVarSA(int numFHs, gsl_rng *r, int C, vector<double>* qBmean, gsl_matrix* qBvar, double* pL) {

	/*
	 * New simulated annealing  method
	 */

	//Set up various parameters
	int Nt = 100; //Steps required before temperature change. NOT A BOTTLENECK!
	int Ns = 20; //Steps required before step vector change.
	int dim = numFHs; //Dimensions
	int n = 2*dim; //Dim mean + dim var = 2*dim
	vector<int> nvec(n,0); //Number of accepts in the n directions
	vector<double> v(n,1); //Step vector, defines the size of steps in all n directions. Start with value of 1.
	double L = -std::numeric_limits<double>::max(); //Likelihood for current step
	double Lopt = -std::numeric_limits<double>::max(); //Optimal likelihood
	double T = 10; //Temperature
	double rT = 0.85; //Temperature decrease coefficient
	vector<double> Lvec; //Best Likelihoods for each temperature value (index k)
	double epsilon = 0.01; //Maximum difference in likelihoods required before stopping
	int Nepsilon = 4; //Number of likelihood values to be compared for stopping scenario, all of which must be within epsilon of each other for stopping
	vector<double> c(n,2); //Coefficient c = {2,2,...,2}. See paper for details.


	//Set up mean and variance
	qBmean->clear();

	//Initialise mean and var
	double value = 0.1; //1e-10; //0.1; //For simplicity, initialise all values to 1/k
	for(int i=0; i<dim; i++) { 

		qBmean->push_back(1.0/dim);
		for(int j=0; j<dim; j++) {

			//Matrix diagonals set to 1/k - everything else set to zero
			if(j==i) {
				gsl_matrix_set(qBvar, i, j, value);
				//gsl_matrix_set(qBvar, i, j, gsl_rng_uniform_pos(r)*100);
			} else {
				gsl_matrix_set(qBvar, i, j, 0.0);
				//gsl_matrix_set(qBvar, i, j, gsl_rng_uniform_pos(r)*100);
			}
		}
	}

  
	//Define "new" parameters
	vector<double> qBmeanNew = *qBmean;
	gsl_matrix* qBvarNew = gsl_matrix_alloc(dim,dim);
	gsl_matrix_memcpy(qBvarNew, qBvar);

	//Define optimal parameters
	vector<double> qBmeanOpt;
	gsl_matrix* qBvarOpt = gsl_matrix_alloc(dim,dim);




	//Loop over temperature changes
	for(int k=0; k<10; k++) { //In theory could be infinite loop, stopping only after a criteria is met
		cout << "k: " << k << "\n";


		for(int m=0; m<Nt; m++) { //Loop over step adjustments. NB, Nt doesn't mean bottlenck here!!
	//		cout << "m: " << m << "\n";

			for(int j=0; j<Ns; j++) { //Steps required before step vector update
				//cout << "j: " << j << "\n";

				for(int h=0; h<n; h++) { //Update each dimension in turn
				//	cout << "h: " << h << "\n";
				//	if(h>dim) { exit(1); } //TEMPORARY!

					/*
					 * Generate new change along dimension h ensuring matrix is positive definite
					 */
					qBmeanNew = *qBmean;
					//cout << "qBmean: "; printDoubleVector(*qBmean);
					gsl_matrix_memcpy(qBvarNew, qBvar);

					if(h < dim) { //Update mean

						do {
							qBmeanNew[h] = (*qBmean)[h] + (gsl_rng_uniform_pos(r)*2 - 1.0);
						} while(qBmeanNew[h] < 0 || qBmeanNew[h] > 1);

						rescale(qBmeanNew);
					//	if(k==0 && m==0 && j==0) { 
					//		cout << "qBmeanNew: "; printDoubleVector(qBmeanNew); 
					//		cout << "qBmean after qBmeanNew change: "; printDoubleVector(*qBmean); 
					//	}

					} else { //Update variance

						bool positiveDefinite = false;
						while(positiveDefinite == false) {

							double oldValue = gsl_matrix_get(qBvar, h-dim, h-dim);
							double newValue;
							do {
								newValue = oldValue + (gsl_rng_uniform_pos(r)*2 - 1.0) * v[h];
							} while(newValue < 10e-20); 
	
							gsl_matrix_set(qBvarNew, h-dim, h-dim, newValue);
				
							//Check if positive-definite
							gsl_matrix* qBvarCholesky = gsl_matrix_alloc(dim,dim);
							gsl_matrix_memcpy(qBvarCholesky, qBvarNew); //Copy into qBvarCholesky as matrix is modified by action below

							//The function below is deprecated by gsl. New function is gsl_linalg_cholesky_decomp1(), but this version of gsl doesn't support that.
							//Hence we use the deprecated version
							gsl_set_error_handler_off();
							int status = gsl_linalg_cholesky_decomp(qBvarCholesky);
							if(status == 0) { //Success indicates that matrix is symmetric positive definite

								positiveDefinite = true;
							} else {
					
								gsl_matrix_memcpy(qBvarNew, qBvar); //Not positive definite, so revert changes
							}
				
							gsl_matrix_free(qBvarCholesky);
						}
					}



					//Compute Likelihood L'(x')
					double Lprime = computeLmeanAndVarPreOnly(&qBmeanNew, qBvarNew, r, C); 
				//	cout << "Lprime: " << Lprime << "\n";

					//Update L if better
					if(Lprime > L) {
		
				//		cout << "Accept!\n";
						//Update values
						*qBmean = qBmeanNew;
						gsl_matrix_memcpy(qBvar, qBvarNew);
						L = Lprime;
						nvec[h]++;

						//If better than optimal values, update these
						if(L > Lopt) {
					
							qBmeanOpt = qBmeanNew;
							gsl_matrix_memcpy(qBvarOpt,qBvarNew);
							Lopt = L;
							cout << "Lopt: " << Lopt << "\n";
							cout << "qBmeanOpt: "; printDoubleVector(qBmeanOpt);
							cout << "qBvarOpt: "; printMatrixMathematica(qBvarOpt);
						}
	
					} else { //Update according to probabilistic framework
						
						double p = exp((Lprime - L)/T);
						double pprime = gsl_rng_uniform(r);
						
						if(pprime < p) { //Accept
		
				//			cout << "Accept by chance!\n";
							//Update values
							*qBmean = qBmeanNew;
							gsl_matrix_memcpy(qBvar, qBvarNew);
							L = Lprime;
							nvec[h]++;
					
						} else { //Reject

				//			cout << "Reject!\n";
							
						}
	
					}
				}
	
			}

			//Update step vector based on accepts - keep around 50%
			vector<double> vprime(n,0);
			for(int h=0; h<n; h++) {

				//Compute new value
				double vprimeh;
				if(nvec[h] > 0.6*Ns) { //Decrease step vector in this direction

					vprimeh = v[h]*(1+c[h]*((nvec[h]/((double)Ns)) - 0.6)/0.4);

				} else if(nvec[h] < 0.4*Ns) { //Increase step vector in this direction
	
					vprimeh = v[h]/(1+c[h]*(0.4-(nvec[h]/((double)Ns)))/0.4);
				} else { //Don't alter
			
					vprimeh = v[h];
				}

				//Update
				v[h] = vprimeh;
				nvec[h] = 0;
			}

		}

		//Reduce the temperature
		T = rT*T;
		
		//Add the best likelihood for temperature index k
		Lvec.push_back(L);

		//Check stopping criteria
		bool stopping = true;
		if(k>=Nepsilon) { //If there are enough data points for checking stop criteria, go ahead and check
		
			for(int u=1; u<=Nepsilon; u++) {

				if(fabs(Lvec[k]-Lvec[k-u]) > epsilon) { stopping = false; } //Check distance from Nepsilon other likelihoods
				if((Lopt - Lvec[k]) > epsilon) { stopping = false; } //Check distance from optimal likelihood
			}
		
			
		} else {
			stopping = false;
		}
	
		if(stopping == true) { break; }	

		//Get ready for next round
		qBmean = &qBmeanOpt;
		gsl_matrix_memcpy(qBvar, qBvarOpt);
		L = Lopt;
	}


	//Save optimal to qBmean and qBvar
	qBmean = &qBmeanOpt;
	gsl_matrix_memcpy(qBvar, qBvarOpt);
	*pL = Lopt; //Set pointer to L to bestL, so can be printed outside this method



	//Clean up
	gsl_matrix_free(qBvarNew);
	gsl_matrix_free(qBvarOpt);

} 



/*
 / Get 100 before frequencies and associated probabilities to be used as an approximation to a sum over all possible combinations.
 / The frequencies are chosen such that the associated likelihood is L > Lmax - 2 where Lmax is the likelihood of the most probable 
 / set of frequencies
 */
void AnalyserPH::getBeforeFrequenciesAndProbabilities(int numFHS, gsl_rng*r, int C, std::vector<double>* logFactStore,std::vector<std::vector<double> > *freqs, std::vector<double> *probs, int numPoints) {
    
    //First find most likely frequencies and associated likelihood
    vector<double> qBopt = findOptimalFreqFHBefore(numFHS,r,C,logFactStore);
    double Lmax = computeLb(qBopt, r, C, logFactStore);
    cout << "Lmax: " << Lmax << " qBopt: "; printDoubleVector(qBopt);
    
	 
    vector<vector<double> > qBs;
    vector<double> Ls;
    vector<double> oldFreqs = qBopt;
    vector<double> newFreqs;
    double L;
    double delta = 0.1;
    
	if(numPoints < 1000) {
    	//Sample 10000 frequency-likelihood sets by a random walk around a space restricted to be within L > L-2
    	for(int i=0;i<10000;i++) {
        
        	do{
            	    do{
                	addRandom(oldFreqs,newFreqs,delta,r); //Create a new vector newFreqs from adding a small delta*rndm(double) to a random index
            	} while(nonZeroDoubleVector(newFreqs)==false); //Ensure that the new frequencies do not have zero elements
            
            	L = computeLb(newFreqs,r,C,logFactStore); //This is a log likelihood

        	}while(L < Lmax-2); //Restrict explored area to be within L > Lmax-2
        
        	qBs.push_back(newFreqs);
        	Ls.push_back(L);
		oldFreqs=newFreqs;
    	}
    
    	//Pick numPoints of these randomly for further computation
    	vector<int> picked;
    	for(int i=0;i<numPoints;i++) {
        	int pick;
        	bool inPicked;
        	do{
            	    pick= gsl_rng_uniform_int(r,qBs.size());
            
            	    inPicked = false;
            	    for(unsigned int i=0; i<picked.size();i++) {
                	if(pick == picked[i]) {
                      	    inPicked = true;
                	}
            	    }
            
        	} while(inPicked == true); //Keep picking a new random integer until a unique one has been found
        
        	//Add picked frequency and likelihood to freqs and probs
       		picked.push_back(pick);
        	freqs->push_back(qBs[pick]);
      	  	probs->push_back(exp(Ls[pick])); //Now save the log likelihood as a likelihood
    	}
    
		rescale(*probs); //Convert likelihoods to probabilities by rescaling

	} else { //Sample the number of points directly

		//Add best value
		//freqs->push_back(qBopt);
		//probs->push_back(Lmax);
		for(int i=0;i<numPoints;i++) {
			
			do{
    			do{
       				addRandom(oldFreqs,newFreqs,delta,r); //Create a new vector newFreqs from adding a small delta*rndm(double) to a random index
 	       		} while(nonZeroDoubleVector(newFreqs)==false); //Ensure that the new frequencies do not have zero elements
        	
				L = computeLb(newFreqs,r,C,logFactStore); //This is a log likelihood

				//	cout << "L: " << L << " freqs: << "; printDoubleVector(newFreqs);
        	} while(L < Lmax-2); //Restrict explored area to be within L > Lmax-2
			
			//cout << "VALUE ADDED AT " << i << "!\n";
        
        	freqs->push_back(newFreqs);
        	probs->push_back(exp(L));
			oldFreqs=newFreqs;
    	}
		
		rescale(*probs); //Convert likelihoods to probabilities by rescaling
	}
}

void AnalyserPH::getBeforeLandscape(int numFHS, gsl_rng*r, int C, std::vector<double>* logFactStore,std::vector<std::vector<double> > *freqs, std::vector<double> *likelihoods, int numPoints, double diffFromMax) {
    
    //First find most likely frequencies and associated likelihood
    vector<double> qBopt = findOptimalFreqFHBefore(numFHS,r,C,logFactStore);
    double Lmax = computeLb(qBopt, r, C, logFactStore);
    cout << "Lmax: " << Lmax << " qBopt: "; printDoubleVector(qBopt);
    
	 
    vector<vector<double> > qBs;
    vector<double> Ls;
    vector<double> oldFreqs = qBopt;
    vector<double> newFreqs;
    double L;
    double delta = 0.1;
    
	if(numPoints < 1000) {
    	//Sample 10000 frequency-likelihood sets by a random walk around a space restricted to be within L > L-2
    	for(int i=0;i<10000;i++) {
        
        	do{
            	    do{
                	addRandom(oldFreqs,newFreqs,delta,r); //Create a new vector newFreqs from adding a small delta*rndm(double) to a random index
            	} while(nonZeroDoubleVector(newFreqs)==false); //Ensure that the new frequencies do not have zero elements
            
            	L = computeLb(newFreqs,r,C,logFactStore); //This is a log likelihood

        	}while(L < Lmax-diffFromMax); //Restrict explored area to be within L > Lmax-2
        
        	qBs.push_back(newFreqs);
        	Ls.push_back(L);
		oldFreqs=newFreqs;
    	}
    
    	//Pick numPoints of these randomly for further computation
    	vector<int> picked;
    	for(int i=0;i<numPoints;i++) {
        	int pick;
        	bool inPicked;
        	do{
            	    pick= gsl_rng_uniform_int(r,qBs.size());
            
            	    inPicked = false;
            	    for(unsigned int i=0; i<picked.size();i++) {
                	if(pick == picked[i]) {
                      	    inPicked = true;
                	}
            	    }
            
        	} while(inPicked == true); //Keep picking a new random integer until a unique one has been found
        
        	//Add picked frequency and likelihood to freqs and probs
       		picked.push_back(pick);
        	freqs->push_back(qBs[pick]);
      	  	likelihoods->push_back(Ls[pick]);
    	}
    

	} else { //Sample the number of points directly

		//Add best value
		//freqs->push_back(qBopt);
		//probs->push_back(Lmax);
		for(int i=0;i<numPoints;i++) {
			
			do{
    			do{
       				addRandom(oldFreqs,newFreqs,delta,r); //Create a new vector newFreqs from adding a small delta*rndm(double) to a random index
 	       		} while(nonZeroDoubleVector(newFreqs)==false); //Ensure that the new frequencies do not have zero elements
        	
				L = computeLb(newFreqs,r,C,logFactStore); //This is a log likelihood

				//	cout << "L: " << L << " freqs: << "; printDoubleVector(newFreqs);
        	} while(L < Lmax-diffFromMax); //Restrict explored area to be within L > Lmax-2
			
			//cout << "VALUE ADDED AT " << i << "!\n";
        
        	freqs->push_back(newFreqs);
       		likelihoods->push_back(L);
			oldFreqs=newFreqs;
    	}
	}
}

//Likelihood for before only
double AnalyserPH::computeLb(vector<double> &freqs, gsl_rng *r, int C, vector<double>* logFactStore) {
    
    double Lb = 0;
    int numDatasets = data->getNumOfDatasets();
	vector<double>& logFactStoreRef = *logFactStore; //Should copy by reference, should save time - consider simplifying in logDirMultProbC instead
    
    for(int i=0; i<numDatasets;i++) { //Loop over datasets (sets of partial haps) in a single simulation
        
        vector<double> inferredFrequencies;
        int numPHs = data->getNumOfPartialHaps(i);
        for(int j=0; j<numPHs; j++) { //Loop over partial haps in a set
            
            double underlyingFreq = 0;
            vector<int> * contribs = data->getContribs(i, j); //Contribs for the jth haplotype
            for(unsigned int k=0; k<contribs->size(); k++) {
                underlyingFreq += freqs[(*contribs)[k]];
            }
            
            if(underlyingFreq == 0) { cout << "ERROR IN computeLb: Underlying freq is zero!\n"; }
            inferredFrequencies.push_back(underlyingFreq);
        }
        rescale(inferredFrequencies); //Think this may be necessary if some haplotypes have been removed during import
        vector<int> obsB = data->getNb(i);
        //Lb += logMultinomialProb(inferredFrequencies, obsB);
        Lb += logDirMultProbC(C, obsB, inferredFrequencies, logFactStoreRef);
    }
    
    return Lb;
}

//Likelihood for after only
double AnalyserPH::computeLa(vector<double> &freqs, gsl_rng *r, int C, vector<double>* logFactStore) {
    
    double La = 0;
    int numDatasets = data->getNumOfDatasets();
	vector<double>& logFactStoreRef = *logFactStore; //Should copy by reference, should save time - consider simplifying in logDirMultProbC instead
    
    for(int i=0; i<numDatasets;i++) { //Loop over datasets (sets of partial haps) in a single simulation
        
        vector<double> inferredFrequencies;
        int numPHs = data->getNumOfPartialHaps(i);
        for(int j=0; j<numPHs; j++) { //Loop over partial haps in a set
            
            double underlyingFreq = 0;
            vector<int> * contribs = data->getContribs(i, j); //Contribs for the jth haplotype
            for(unsigned int k=0; k<contribs->size(); k++) {
                underlyingFreq += freqs[(*contribs)[k]];
            }
            
            if(underlyingFreq == 0) { cout << "ERROR IN computeLb: Underlying freq is zero!\n"; }
            inferredFrequencies.push_back(underlyingFreq);
        }
        rescale(inferredFrequencies); //Think this may be necessary if some haplotypes have been removed during import
        vector<int> obsA = data->getNa(i);
        La += logDirMultProbC(C, obsA, inferredFrequencies, logFactStoreRef);
    }
    
    return La;
}


//Input mean of full dimensionality	
double AnalyserPH::computeLmeanPreOnly(vector<double>* qBmean, gsl_rng* r, double C, bool print) {

	double L = 0;
	int numDatasets = data->getNumOfDatasets();
	int dimFH = (int) qBmean->size();


	for(int i=0; i<numDatasets; i++) { //Loop over sets of partial haps

		if(print == true) {
			cout << "SOPH " << i << "\n";
		}
		int dimPH = data->getNumOfPartialHaps(i);
        	vector<int> obsB = data->getNb(i);
		int Nb = sumOfVector(obsB);

		vector<vector<int> > contribsPHset;
		for(int j=0; j<dimPH; j++) {

			contribsPHset.push_back(*data->getContribs(i,j));
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
		 * First, transform qBmean into qBmean_PH=TqBmean
		 */
		//Turn qBmean into gsl vector
		gsl_vector *qBmeanGSL = gsl_vector_alloc(dimFH);
		for(int i=0;i<dimFH;i++) {
			gsl_vector_set(qBmeanGSL,i,(*qBmean)[i]);
		}

		//Create TqBmean through matrix vector multiplication
		gsl_vector *TqBmeanGSL = gsl_vector_alloc(dimPH);
		gsl_blas_dgemv(CblasNoTrans, 1.0, T, qBmeanGSL, 0.0, TqBmeanGSL);

		//Convert back into C++ vector
		vector<double> TqBmean;
		vector<double> NbTqBmean;
		for(int i=0; i<dimPH; i++) {
			TqBmean.push_back(gsl_vector_get(TqBmeanGSL, i));
			NbTqBmean.push_back(Nb*TqBmean[i]);
		}

	
		/*
		 * Next, create the variance expression 
		 */
		//Construct epsilonNbM(TqBmean) with epsilon = (Nb+C)/(1+C)
		gsl_matrix* epsilonNbMTqBmean = gsl_matrix_alloc(dimPH,dimPH);
		constructMatrixM(TqBmean, epsilonNbMTqBmean);
		double epsilon = (Nb+C)/((double)(1+C));
		gsl_matrix_scale(epsilonNbMTqBmean, epsilon*Nb);

		
		double LphSet = - numeric_limits<double>::max();

		if(dimPH > 1) { //Most common case, i.e. multiple partial haps in a set
		
			//Reduce dimensionality by 1 to ensure non-degeneracy
			NbTqBmean.pop_back();
			gsl_matrix_view qBvarPHreducedView = gsl_matrix_submatrix(epsilonNbMTqBmean, 0,0, dimPH-1, dimPH-1); //Remove row k and column k
			gsl_matrix * qBvarPH = gsl_matrix_alloc(dimPH-1, dimPH-1);
			gsl_matrix_memcpy(qBvarPH, &(qBvarPHreducedView.matrix)); //Copy the reduced matrix view into a "normal" matrix varXaPH
			vector<double> obsBreduced;
			for(int i=0; i<dimPH-1; i++) {
				obsBreduced.push_back(obsB[i]);
	        	}


			if(dimPH>2) { //Multivariate approach

				LphSet = logMultivariateNormalPDF(obsBreduced,NbTqBmean, qBvarPH);
		
			} else if(dimPH==2) {

				//Univariate case
				double qBmeanUni = NbTqBmean[0];
				double qBvarUni = gsl_matrix_get(qBvarPH,0,0);
				double stDevqBvarUni = sqrt(qBvarUni);

				LphSet = -(obsB[0]-qBmeanUni)*(obsB[0]-qBmeanUni)/(2*qBvarUni) - log(stDevqBvarUni*sqrt(2*M_PI)); //log L
	
			}

			//Deallocate variables created in this scope			
			gsl_matrix_free(qBvarPH);

		} else if (dimPH == 1) { 
			//Case where there is a single partial haplotype in a set.
			//This can be thought of as a case with a single observed partial haplotype
			//and one (or several) unobserved haplotypes.
			//If all loci in the partial haplotype set are variants (i.e. not monomorphic)
			//then there MUST exists at least one other haplotype containing the unobserved alleles,
			//however, this (or these) haplotype(s) are not observed.
			//As such, WLOG we can consider this a case of dimPH=2, which gets reduced to a one dimensional system
			//under reduction, i.e. similar to the dimPH==2 case of above.
				
			double qBmeanUni = NbTqBmean[0];
			double qBvarUni = gsl_matrix_get(epsilonNbMTqBmean,0,0);
			double stDevqBvarUni = sqrt(qBvarUni);

			LphSet = -(obsB[0]-qBmeanUni)*(obsB[0]-qBmeanUni)/(2*qBvarUni) - log(stDevqBvarUni*sqrt(2*M_PI)); //log L
		
	
		} else {

			cout << "ERROR in computeLmeanAndVarPreOnly: dimPH is <1\n";
		}


		//Check if L is still negative infinity (shouldn't happen)
		if(LphSet == -numeric_limits<double>::max() || ::isnan(LphSet) == true) { 
		
				cout << "ERROR in AnalyserPH::computeLmeanAndVarPreOnly: Likelihood is negative infinity.\n";
		}
		else{

			L += LphSet;
		}

		//Clean up
		gsl_vector_free(qBmeanGSL);
		gsl_vector_free(TqBmeanGSL);
		gsl_matrix_free(T);
		gsl_matrix_free(epsilonNbMTqBmean);

	}

	return L;
} 




//Input mean and variance are of full dimensionality	
double AnalyserPH::computeLmeanAndVarPreOnly(vector<double>* qBmean, gsl_matrix* qBvar, gsl_rng* r, double C) {

	double L = 0;
	int numDatasets = data->getNumOfDatasets();
	int dimFH = (int) qBmean->size();


	for(int i=0; i<numDatasets; i++) { //Loop over sets of partial haps

		int dimPH = data->getNumOfPartialHaps(i);
        	vector<int> obsB = data->getNb(i);
		int Nb = sumOfVector(obsB);

		vector<vector<int> > contribsPHset;
		for(int j=0; j<dimPH; j++) {

			contribsPHset.push_back(*data->getContribs(i,j));
		}


		//Construct matrix T
		gsl_matrix* T = gsl_matrix_alloc(dimPH,dimFH);
		for(int j=0; j<dimPH; j++) { //Loop over rows 

		
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
		 * First, transform qBmean into qBmean_PH=TqBmean
		 */
		//Turn qBmean into gsl vector
		gsl_vector *qBmeanGSL = gsl_vector_alloc(dimFH);
		for(int i=0;i<dimFH;i++) {
			gsl_vector_set(qBmeanGSL,i,(*qBmean)[i]);
		}

		//Create TqBmean through matrix vector multiplication
		gsl_vector *TqBmeanGSL = gsl_vector_alloc(dimPH);
		gsl_blas_dgemv(CblasNoTrans, 1.0, T, qBmeanGSL, 0.0, TqBmeanGSL);

		//Convert back into C++ vector
		vector<double> TqBmean;
		vector<double> NbTqBmean;
		for(int i=0; i<dimPH; i++) {
			TqBmean.push_back(gsl_vector_get(TqBmeanGSL, i));
			NbTqBmean.push_back(Nb*TqBmean[i]);
		}

		/*
		 * Next, create the variance expression 
		 */
		//Construct epsilonNbM(TqBmean) with epsilon = (Nb+C)/(1+C)
		gsl_matrix* epsilonNbMTqBmean = gsl_matrix_alloc(dimPH,dimPH);
		constructMatrixM(TqBmean, epsilonNbMTqBmean);
		double epsilon = (Nb+C)/((double)(1+C));
		gsl_matrix_scale(epsilonNbMTqBmean, epsilon*Nb);

		//Construct Nb(Nb-epsilon)TqBvarT^-1
		gsl_matrix* TqBvar = gsl_matrix_alloc(dimPH,dimFH);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, qBvar, 0.0, TqBvar);
		gsl_matrix* NbNbMinusEpsilonTqBvarTTrans = gsl_matrix_alloc(dimPH, dimPH);
		gsl_blas_dgemm(CblasNoTrans, CblasTrans,1.0, TqBvar, T, 0.0, NbNbMinusEpsilonTqBvarTTrans);
		gsl_matrix_scale(NbNbMinusEpsilonTqBvarTTrans, Nb*(Nb-epsilon));

		//Add the two terms together. Final matrix stored in first term
		gsl_matrix_add(epsilonNbMTqBmean, NbNbMinusEpsilonTqBvarTTrans);


		double LphSet = - numeric_limits<double>::max();

		if(dimPH > 1) { //Most common case, i.e. multiple partial haps in a set
		
			//Reduce dimensionality by 1 to ensure non-degeneracy
			NbTqBmean.pop_back();
			gsl_matrix_view qBvarPHreducedView = gsl_matrix_submatrix(epsilonNbMTqBmean, 0,0, dimPH-1, dimPH-1); //Remove row k and column k
			gsl_matrix * qBvarPH = gsl_matrix_alloc(dimPH-1, dimPH-1);
			gsl_matrix_memcpy(qBvarPH, &(qBvarPHreducedView.matrix)); //Copy the reduced matrix view into a "normal" matrix varXaPH
			vector<double> obsBreduced;
			for(int i=0; i<dimPH-1; i++) {
				obsBreduced.push_back(obsB[i]);
	        	}


			if(dimPH>2) { //Multivariate approach

				LphSet = logMultivariateNormalPDF(obsBreduced,NbTqBmean, qBvarPH);
		
			} else if(dimPH==2) {

				//Univariate case
				double qBmeanUni = NbTqBmean[0];
				double qBvarUni = gsl_matrix_get(qBvarPH,0,0);
				double stDevqBvarUni = sqrt(qBvarUni);

				LphSet = -(obsB[0]-qBmeanUni)*(obsB[0]-qBmeanUni)/(2*qBvarUni) - log(stDevqBvarUni*sqrt(2*M_PI)); //log L
	
			}

			//Deallocate variables created in this scope			
			gsl_matrix_free(qBvarPH);

		} else if (dimPH == 1) { 
			//Case where there is a single partial haplotype in a set.
			//This can be thought of as a case with a single observed partial haplotype
			//and one (or several) unobserved haplotypes.
			//If all loci in the partial haplotype set are variants (i.e. not monomorphic)
			//then there MUST exists at least one other haplotype containing the unobserved alleles,
			//however, this (or these) haplotype(s) are not observed.
			//As such, WLOG we can consider this a case of dimPH=2, which gets reduced to a one dimensional system
			//under reduction, i.e. similar to the dimPH==2 case of above.
				
			double qBmeanUni = NbTqBmean[0];
			double qBvarUni = gsl_matrix_get(epsilonNbMTqBmean,0,0);
			double stDevqBvarUni = sqrt(qBvarUni);

			LphSet = -(obsB[0]-qBmeanUni)*(obsB[0]-qBmeanUni)/(2*qBvarUni) - log(stDevqBvarUni*sqrt(2*M_PI)); //log L
		
	
		} else {

			cout << "ERROR in computeLmeanAndVarPreOnly: dimPH is <1\n";
		}


		//Check if L is still negative infinity (shouldn't happen)
		if(LphSet == -numeric_limits<double>::max() || ::isnan(LphSet) == true) { 
		
				cout << "ERROR in AnalyserPH::computeLmeanAndVarPreOnly: Likelihood is negative infinity.\n";
		}
		else{

			L += LphSet;
		}

		//Clean up
		gsl_vector_free(qBmeanGSL);
		gsl_vector_free(TqBmeanGSL);
		gsl_matrix_free(T);
		gsl_matrix_free(epsilonNbMTqBmean);
		gsl_matrix_free(TqBvar);
		gsl_matrix_free(NbNbMinusEpsilonTqBvarTTrans);

	}

	return L;
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
//RepIndex and geneIndex optional. If only one int is supplied, then geneIndex = int, repIndex = -1, hence ordering.
vector<Sequence > AnalyserPH::filterHaps3(gsl_rng* r, double C, string filePathFolder /* = ""*/, int geneIndex /* = -1 */, int repIndex /* = -1 */) {

	double cutOff = 2e-10; //Minimum cut-off that should be enforced at all times

	vector<Sequence> fullHaps;
	//Assumes dataPH have been loaded into AnalyserPH object
	if(data->getNumOfDatasets() > 0) { //Data present for this gene. Probably don't need to check here. Should be checked upstream.

		//Generate full haps from partial haps
		vector<Sequence> phs = data->getSequences(); //Get partial haplotypes
		constructFullHaps(&phs, &fullHaps);
		data->computeContribs(fullHaps); //Updates contribs wrt full haplotypes

		//If an output folder is provided, print the number of initial (i.e. potential) haplotypes
		if(filePathFolder.compare("") != 0) { //Only print if an output folder is supplied
			stringstream ss;
			if(repIndex != -1) {
				ss << filePathFolder << "numPotentialHaps_rep" << repIndex << ".dat";
			} else {
				ss << filePathFolder << "numPotentialHaps.dat";
			}
			string filePath = ss.str();
			ofstream outputFile;
			
			//If output file doesn't exist, create it
			if(fileExists(filePath)!= true) {
				outputFile.open(filePath.c_str());
			} else {
				outputFile.open(filePath.c_str(),ios_base::app); //Otherwise open in "append" mode
			}

			outputFile << fullHaps.size() << "\n";
			outputFile.close();
		}

		//Infer frequencies
		vector<double> qBmean;
		double LqBgene;
		optimiseqBmeanFromPre((int)fullHaps.size(), r, C, &qBmean, &LqBgene); 

		cout << "Number of haplotypes before filtering for replicate " << repIndex << " gene " << geneIndex << " : " << fullHaps.size() << "\n";

		//Sort haplotypes by size of frequencies
		typedef pair<double, int> hapFreq; //Temporary container for frequency and haplotype index
		vector<hapFreq> hapFreqs;
		for(unsigned int i=0; i<fullHaps.size(); i++) {

			hapFreqs.push_back(hapFreq(qBmean[i],i));
		}
		sort(hapFreqs.begin(), hapFreqs.end()); //Sorts by frequency, small to large


		vector<int> toBeRemoved;
		int hapsLimit = 100; //Less than 100 haplotypes needed
		bool hapsLimitReached = false;
		int kept = 0;
		double currentCutOff = cutOff; //Real currentCutOff may be lower in reality, but unimportant to us
		vector<Sequence> fullHapsOld = fullHaps;
		for(unsigned int i=0; i<hapFreqs.size(); i++) {
	
			cout << "Index " << i << " freq: " << hapFreqs[i].first << "\n";
			cout << "currentCutOff: " << currentCutOff << "\n";

			//Check that sufficiently many haplotypes have been removed such that #intial - #removed < limit
			if((int) (fullHaps.size() - toBeRemoved.size()) <= hapsLimit && hapsLimitReached == false) {

				cout << "Haps limit was reached:\n";
				cout << "Num removed: " << toBeRemoved.size() << "\n";
				cout << "Num kept: " << kept << "\n";
				cout << "Might continue removing haplotypes if cutOff limit hasn't been reached.\n";

				hapsLimitReached = true;

				//Print the required frequency for reaching the haps limit
				if(filePathFolder.compare("") != 0) { //Only print if an output folder is supplied
					stringstream ss;
					if(repIndex != -1 && geneIndex != -1) {
						ss << filePathFolder << "requiredFrequencyHapsLimit_" << hapsLimit << "_rep" << repIndex << "_gene" << geneIndex << ".dat";
					} else if(geneIndex != -1) {
						ss << filePathFolder << "requiredFrequencyHapsLimit_" << hapsLimit << "_gene" << geneIndex << ".dat";
					} else {
						ss << filePathFolder << "requiredFrequencyHapsLimit_" << hapsLimit << ".dat";
					}
					string filePath = ss.str();
					ofstream outputFile;
					outputFile.open(filePath.c_str());
					outputFile << currentCutOff;
					outputFile.close();
				}

			}

			//If haps limit has been reached, also check that we have removed haplotypes up to the cut-off limit, 
			//e.g. we may have met the hapsLimit but still have frequencies of < cutOff.
			//In this case we want to keep removing haplotypes until we have removed all haplotypes below the cut-off.
			currentCutOff = hapFreqs[i].first; //Get the current frequency cutoff for haplotype i
			if(hapsLimitReached == true && currentCutOff > cutOff) {

				//Save new haplotypes
				fullHaps = fullHapsOld;
				break;
			}


			/*
			 * If limits haven't been reached, continue/start removing haplotypes
			 */ 

			//Add the next haplotype to new list of removed haplotypes
			vector<int> toBeRemovedNew = toBeRemoved;
			toBeRemovedNew.push_back(hapFreqs[i].second);

			//Sort in ascending order and remove from fullHapsGeneNew						
			sort(toBeRemovedNew.begin(), toBeRemovedNew.end()); //Ascending order
			vector<Sequence> fullHapsNew = fullHaps;
			for(int j=((int) toBeRemovedNew.size()) -1; j>=0; j--) {

				fullHapsNew.erase(fullHapsNew.begin() + toBeRemovedNew[j]);	
			}

			//Check if fullHapsNew can de be decomposed to represent all partial haplotypes
			if(data->decomposeAndCheck(fullHapsNew) == true) {

				//Keep change
				toBeRemoved = toBeRemovedNew;
				cout << "Haplotype " << i << " removed.\n";
			} else {

				//Don't keep change
				cout << "Haplotype " << i << " kept.\n";
				kept++;
			}

			//Next round new = old
			fullHapsOld = fullHapsNew;

		}

		if(hapsLimitReached == true) {

			cout << "Haps and cut-off limit were reached. Continuing.\n";		
			
			//I don't think this is necessary as recomputed downstream.
			data->computeContribs(fullHaps); //Updates contribs wrt full haplotypes

		} else {

			cout << "Haps limit wasn't reached, so exiting.\n";
			exit(1);
		}
	}

	return fullHaps;
}

void AnalyserPH::analyseLeonard(vector<double>& logLNts) {

	cout << "Before computing freqs\n";
	vector<double> donorFreqs = data->getDonorSLFreqs();
	vector<double> recipientFreqs = data->getRecipientSLFreqs();	

	cout << "DonorFreqs: "; printDoubleVector(donorFreqs);
	cout << "RecipientFreqs: "; printDoubleVector(recipientFreqs);

	double qCut = 0; //0% frequency cut-off, i.e. use only the standard Leonard method, not the cumulative one

	for(unsigned int i=0; i<donorFreqs.size(); i++) { //Loop over variant loci

		if(recipientFreqs[i] >= qCut) { //Standard scenario

			for(int Nt=1; Nt<=((int) logLNts.size()); Nt++) {

				double LNt = 0;
				for(int j=0; j<=Nt; j++) {

					double temp;
					if(j==0) {
						temp = gsl_ran_beta_pdf(recipientFreqs[i],1e-10,Nt-j)*gsl_ran_binomial_pdf(j,donorFreqs[i],Nt);
				

					} else if(j==Nt) {

						temp = gsl_ran_beta_pdf(recipientFreqs[i],j,1e-10)*gsl_ran_binomial_pdf(j,donorFreqs[i],Nt);

						//The gsl_ran_beta_pdf(recipientFreqs[i],j,1e-10) can produce infinite or nan likelihood when
						//recipientFreqs[i] == 1. This is because the distribution is defined for [0,1] only
						//when a,b>=1. For a,b<1 then the distribution is only defined on (0,1). Regardless, even if
						//recipientFreqs[i] was set to 0.999999 this would still result in a very large value when b<1
						//as this distribution has an infinity tending tail as x -> 1. To avoid issues, we set the
						//likelihood to 1000 here. This is still a very large value, but it is finite. We obtain
						//the required dynamics.
						if(::isinf(temp)==true || ::isnan(temp)==true) {

							if(recipientFreqs[i]==1) {
		
								temp = 1000; //Some large value
							}			
						}
					} else {

						//LNt += gsl_ran_beta_pdf(recipientFreqs[i],j,Nt-j)*gsl_ran_binomial_pdf(j,donorFreqs[i],Nt);
						temp = gsl_ran_beta_pdf(recipientFreqs[i],j,Nt-j)*gsl_ran_binomial_pdf(j,donorFreqs[i],Nt);
					}
					//cout << "temp: " << temp << "\n";
					LNt += temp;
				}
				logLNts[Nt-1]+=log(LNt);
			}

		} else { //Cumulative scenario, i.e. we can't be sure if the observation in the recipient is really there

			for(int Nt=1; Nt<=((int) logLNts.size()); Nt++) {

				double LNt = 0;
				for(int j=0; j<=Nt; j++) {

					double temp;
					if(j==0) {
						temp = gsl_cdf_beta_P(qCut,1e-10,Nt-j)*gsl_ran_binomial_pdf(j,donorFreqs[i],Nt);
				
					} else if(j==Nt) {

						temp = gsl_cdf_beta_P(qCut,j,1e-10)*gsl_ran_binomial_pdf(j,donorFreqs[i],Nt);
					} else {

						//LNt += gsl_cdf_beta_P(qCut,j,Nt-j)*gsl_ran_binomial_pdf(j,donorFreqs[i],Nt);
						temp = gsl_cdf_beta_P(qCut,j,Nt-j)*gsl_ran_binomial_pdf(j,donorFreqs[i],Nt);
					}
					//cout << "temp: " << temp << "\n";
					LNt += temp;
				}
				logLNts[Nt-1]+=log(LNt);
			}
		}
	}	

}


void AnalyserPH::analyseLeonardExact(vector<double>& logLNts) {

	vector<double> donorFreqs = data->getDonorSLFreqs();
	vector<int> recipientObs = data->getRecipientSLObs();
	vector<int> recipientObsTot = data->getRecipientSLObsTot();
	cout << "DonorFreqs: "; printDoubleVector(donorFreqs);
	cout << "RecipientObs: "; printIntVector(recipientObs);
	cout << "RecipientObsTot: "; printIntVector(recipientObsTot);

	double qCut = 0; //0% frequency cut-off, i.e. use only the standard Leonard method, not the cumulative one

	for(unsigned int i=0; i<donorFreqs.size(); i++) { //Loop over variant loci

		if(recipientObs[i] >= qCut*recipientObsTot[i]) { //Standard scenario

			for(int Nt=1; Nt<=((int) logLNts.size()); Nt++) {

				double LNt = 0;
				for(int j=0; j<=Nt; j++) {

					double L;
					if(j==0) {
						L = exp(logBetaBinProb(1e-10,Nt-j,recipientObs[i],recipientObsTot[i]))*gsl_ran_binomial_pdf(j,donorFreqs[i],Nt);
				
					} else if(j==Nt) {

						L = exp(logBetaBinProb(j,1e-10,recipientObs[i],recipientObsTot[i]))*gsl_ran_binomial_pdf(j,donorFreqs[i],Nt);
					} else {

						L = exp(logBetaBinProb(j,Nt-j,recipientObs[i],recipientObsTot[i]))*gsl_ran_binomial_pdf(j,donorFreqs[i],Nt);
					}
				//	cout << "L: " << L << "\n";
				
					//Sometimes logBetaBinProb gets so small (e.g. if logBetaBin=-500 then exp(logBetaBin)~10^-218) that taking the exponential yields an outcome of zero.
					//These are very unlikely scenarios so we don't really care too much about precision here. As such, we truncate L at 10^-20:
					if(L < 1e-20) { L=1e-20; }
					LNt += L;

				}
				logLNts[Nt-1]+=log(LNt);
			}

		} else { //Cumulative scenario, i.e. we can't be sure if the observation in the recipient is really there

			for(int Nt=1; Nt<=((int) logLNts.size()); Nt++) {

				double LNt = 0;
				for(int j=0; j<=Nt; j++) {

					double L;
					if(j==0) {
						L = BetaBinCDF(1e-10,Nt-j,qCut*recipientObsTot[i],recipientObsTot[i])*gsl_ran_binomial_pdf(j,donorFreqs[i],Nt);
				
					} else if(j==Nt) {

						L = BetaBinCDF(j,1e-10,qCut*recipientObsTot[i],recipientObsTot[i])*gsl_ran_binomial_pdf(j,donorFreqs[i],Nt);
					} else {

						L = BetaBinCDF(j,Nt-j,qCut*recipientObsTot[i],recipientObsTot[i])*gsl_ran_binomial_pdf(j,donorFreqs[i],Nt);
					}
				//	cout << "L: " << L << "\n";
					LNt += L;
				}
				logLNts[Nt-1]+=log(LNt);
			}
		}
	}	

}










