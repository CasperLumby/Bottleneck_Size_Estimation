//
//  dataPHgen.cpp
//
//

#include "dataPHgen.hpp"
#include <stdio.h>
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "misc.h"
#include <math.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <algorithm> //Apparently necessary for sort command
#include <sys/stat.h> //For folder management

using namespace std;

//Class related methods + constructors
DataPHgen::DataPHgen() { } //Constructor does nothing
DataPHgen::~DataPHgen() {} //Deconstructor

void DataPHgen::readData(Path * p) {

	//Load paths from path object
	path = static_cast<PathPH *>(p);
	string multi_traj_path = path->getMultiLocusTrajFilePath(); //TODO: Implement this
	
	//Only try and add data if present. If not data, DataPHMGgen takes care of this.
	ifstream multi_trajs_file;
	vector<partialHap> phs;
	if(multi_traj_path.empty() != true) {
	
		multi_trajs_file.open(multi_traj_path.c_str());
		string line;
		while(getline(multi_trajs_file, line)) {
		
			//Read in all the partial haps, then sort into sophs later
			partialHap ph;
			ph.readTraj(line); //See method for assumptions on input format
			phs.push_back(ph);

			//Check if ph covers new position, if so, add to physicalPos
			vector<int> posCurrent = ph.getPos();
			for(unsigned int i=0; i<posCurrent.size(); i++) {
		
				bool inPhysPos = false;
				for(unsigned int j=0; j<physicalPos.size(); j++) {
					
					if(physicalPos[j] == posCurrent[i]) {
						inPhysPos = true; break;
					}
				}
		
				if(inPhysPos == false) {
					
					physicalPos.push_back(posCurrent[i]);
					sort(physicalPos.begin(), physicalPos.end()); //Always sort when adding an element
				}		
			}
		}

		//Sort partial haps into sophs based on the loci they cover
		SOPHs.clear();
		for(unsigned int i=0; i<phs.size(); i++) {

			vector<int> posCurrent = phs[i].getPos();
			Sequence seqCurrent = phs[i].getSeq();

			//First we update the sequence of the haplotype so it matches the length of physicalPos
			//E.g. sequence might be ACT covering pos 300, 400, 500, but needs to be --ACT-
			//if physPos are 100, 200, 300, 400, 500, 600
			int minPosPH = posCurrent[0];
			int maxPosPH = posCurrent.back();
			for(unsigned int k=0; k<physicalPos.size(); k++) {
				if(physicalPos[k] < minPosPH) {
					seqCurrent.addBaseFront('-');
							
				} else if(physicalPos[k] > maxPosPH) {
					seqCurrent.addBase('-');
				}
							
			}
			phs[i].setSeq(seqCurrent); //Set the partial haplotype seq to the updated seq


			//Then check which soph it belongs to
			bool phInSOPH = false;
			for(unsigned int j=0; j<SOPHs.size(); j++) {

				vector<int> sophCurrentPos = SOPHs[j].getPos();
				if(sophCurrentPos.size() == posCurrent.size()) { //Cover the same number of positions

					if(identicalIntVectors(sophCurrentPos,posCurrent) == true) { //Cover the same positions

						//Partial hap matches with soph, so add it.
						phInSOPH = true;
						SOPHs[j].addPartialHap(phs[i]);
						SOPHs[j].countObsTot();
						break;
					}
				}
			}

			if(phInSOPH == false) { //Partial hap not found, so create new soph

				setOfPartialHaps newSOPH;
				newSOPH.addPartialHap(phs[i]);
				newSOPH.setPos(posCurrent);
				newSOPH.countObsTot();

				//Find the loci covered for this soph
				vector<int> lociCovered;
				for(int j=0; j<seqCurrent.getLength(); j++) {

					if(seqCurrent.getBase(j) != '-') {

						lociCovered.push_back(j);
					}
				}
				newSOPH.setLociCovered(lociCovered);
				SOPHs.push_back(newSOPH);
			}

		}
	}


}
void DataPHgen::simulateData(SimParam *sp) { } //Implement later


//Assumes a string of the form: Number of Loci [Set of loci] Variant #times [Time #observations]
//Assumes everything is ordered in increasing size, e.g. Loci = 100, 200, 300, 400 ...
//and [Time #observations] = [1 23] [2 44] [4 37] ... (for time points 1,2,4)
void DataPHgen::partialHap::readTraj(string &trajString) {

	//Read traj string into vector of "words"
	stringstream ss(trajString);
	istream_iterator<std::string> begin(ss);
	istream_iterator<std::string> end;
	vector<string> words(begin, end);

	//Get the loci positions
	int numLoci = stoi(words[0]);
	pos.clear();
	for(int i=1; i<(numLoci+1); i++) {

		pos.push_back(stoi(words[i]));
	}

	seq = Sequence(words[numLoci+1]);
	int numTimePoints = stoi(words[numLoci+2]);

	times.clear();
	obs.clear();
	for(int i=numLoci+3; i<(numLoci+3+2*numTimePoints); i++) {

		times.push_back(stoi(words[i]));
		i++;
		obs.push_back(stoi(words[i]));
		
	}

}

Sequence DataPHgen::partialHap::getSeq() { return seq; }
vector<int> DataPHgen::partialHap::getPos() { return pos; }
int DataPHgen::partialHap::getObs(int obsIndex) { return obs[obsIndex]; }
void DataPHgen::partialHap::setSeq(Sequence &s) { seq = s; }
int DataPHgen::partialHap::getNumOfTimePoints() { return (int) times.size(); }
void DataPHgen::partialHap::clearContributions() { contribs.clear(); }
void DataPHgen::partialHap::addContrib(int c) { contribs.push_back(c); }
vector<int>* DataPHgen::partialHap::getContribs() { return &contribs; }
void DataPHgen::partialHap::print() {

	cout << "Sequence: ";
	seq.print();

	cout << "Observations: "; printIntVector(obs);

	//cout << "Contributions ";
	//printIntVector(contribs);

}

void DataPHgen::setOfPartialHaps::addPartialHap(partialHap &p) { phVec.push_back(p); }
vector<int> DataPHgen::setOfPartialHaps::getPos() { return pos; }
void DataPHgen::setOfPartialHaps::setPos(vector<int> &p) { pos = p; }
int DataPHgen::setOfPartialHaps::getNumOfPartialHaps() { return (int) phVec.size(); }
vector<int> DataPHgen::setOfPartialHaps::getObs(int timeIndex) { //Get vector of observations from all PHs in SOPH at time point timeIndex

	vector<int> obsT;
	for(unsigned int i=0; i<phVec.size(); i++) {

		obsT.push_back(phVec[i].getObs(timeIndex));
	}
	return obsT;
}
void DataPHgen::setOfPartialHaps::addContrib(int index, int contrib) { phVec[index].addContrib(contrib); }
vector<int>* DataPHgen::setOfPartialHaps::getContribs(int index) { return phVec[index].getContribs(); }
void DataPHgen::setOfPartialHaps::clearContributions() {
	for(unsigned int i=0;i<phVec.size();i++) {
		phVec[i].clearContributions();
	}
}
void DataPHgen::setOfPartialHaps::countObsTot() {

	obsTot.clear();
	if(phVec.size() > 0) {
		for(int i=0; i<phVec[0].getNumOfTimePoints(); i++) { //Loop over time points

			obsTot.push_back(0);
	
			for(unsigned int j=0; j<phVec.size(); j++) {
				obsTot[i] += phVec[j].getObs(i);
			}
		}

	} else { //No data!

		cout << "No data for this soph, so can't compute obsTot. Investigate. Exiting.\n";
		exit(1);
	}
}

void DataPHgen::setOfPartialHaps::print() {

	for(unsigned int i=0; i<phVec.size(); i++) {
		phVec[i].print();
	}

	cout << "Total observations: "; printIntVector(obsTot);
}
void DataPHgen::setOfPartialHaps::setLociCovered(vector<int> &lc) { lociCovered = lc; }
int DataPHgen::getNumOfSOPHs() { return (int) SOPHs.size(); }
int DataPHgen::getNumOfTimePoints() { return SOPHs[0].getNumOfTimePoints(); }
DataPHgen::setOfPartialHaps* DataPHgen::getSOPH(int SOPHindex) { return &(SOPHs[SOPHindex]); }



vector<Sequence> DataPHgen::getSequences() {

	vector<Sequence> allSeqs;
	for(unsigned int i=0;i<SOPHs.size();i++) {
		vector<Sequence> phsSeqs = SOPHs[i].getSequences();
		for(unsigned int j=0; j<phsSeqs.size();j++) {
			allSeqs.push_back(phsSeqs[j]);
		}
	}

	return allSeqs;
}

vector<Sequence> DataPHgen::setOfPartialHaps::getSequences() {

	vector<Sequence> allSeqs;
	for(unsigned int i=0;i<phVec.size();i++) {
		allSeqs.push_back(phVec[i].getSeq());
	}

	return allSeqs;
}

int DataPHgen::setOfPartialHaps::getNumOfTimePoints() { return phVec[0].getNumOfTimePoints(); }


void DataPHgen::computeContribs(vector<Sequence> &fhs) {

	for(unsigned int i=0; i<SOPHs.size(); i++) {

		setOfPartialHaps* soph = &(SOPHs[i]);
		soph->clearContributions();
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


void DataPHgen::print() {
	cout << "Printing " << SOPHs.size() << "sophs:\n";
	for(unsigned int i=0;i<SOPHs.size();i++) {
		SOPHs[i].print();
	}
}
