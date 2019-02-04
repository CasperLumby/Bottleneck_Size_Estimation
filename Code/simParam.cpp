
//Forward declared dependencies

//Included dependencies
#include "simParam.h"
#include <stdio.h>
#include <iostream>

using namespace std;

//Constructors
SimParam::SimParam() : Nb(1000), Na(1000), numLoci(3), Nt(100), numGenes(1), numGenesWithData(1), geneLength(-1), seed(1), withinHostSelectionPresent(false), useExternalReadStatistic(false), removeMonomorphicSim(true), filterHaplotypesMethod(0) {
    cout << "Param instantiated with default values: Na=1000, Nb=1000, numLoci=3, Nt=100, numGenes=1., seed=1\n";
}
SimParam::~SimParam() {} //Deconstructor


/*
 * Setters
 */
void SimParam::setNa(int NA) { Na = NA; }
void SimParam::setNb(int NB) { Nb = NB; }
void SimParam::setNumLoci(int num) { numLoci = num; }
void SimParam::setNt(int NT) { Nt = NT; }
void SimParam::setNumGenes(int num) { numGenes = num; }
void SimParam::setNumGenesWithData(int num) { numGenesWithData = num; }
void SimParam::setSeed(unsigned long int s) { seed = s; }
void SimParam::setOutputFile(ofstream* o) { outputFile = o; }
void SimParam::setNumReps(int num) { numReps = num; }
void SimParam::setNumLociPerGene(int num) { numLociPerGene = num; }
void SimParam::setWithinHostSelectionPresent(bool b) { withinHostSelectionPresent = b; }
void SimParam::setUseExternalReadStatistic(bool b) { useExternalReadStatistic = b; }
void SimParam::setGeneLength(int gl) { geneLength = gl; }
void SimParam::setFilterHaplotypesMethod(int fhm) { filterHaplotypesMethod = fhm; }
void SimParam::setRemoveMonomorphicSim(bool b) { removeMonomorphicSim = b; }

//Transmission
void SimParam::setSel(Sequence s) { sel = s; }
void SimParam::setSelMagVec(vector<double> & SMV) { selMagVec = SMV; }
void SimParam::setSeqVec(vector<Sequence> & SV) { seqVec = SV; }
void SimParam::setSelMulti(vector<Sequence> & SM) { selMulti = SM; }
void SimParam::setSelMagVecMulti(vector<vector<double> > & SMVM) { selMagVecMulti = SMVM; }
void SimParam::setSeqVecMulti(vector<vector<Sequence> > & SVM) { seqVecMulti = SVM; }
void SimParam::setEpiPosVec(vector<vector<int> > &EPV) { epiPosVec = EPV; }
void SimParam::setEpiMagVec(vector<double> &EMV) { epiMagVec = EMV; }
void SimParam::setEpiPosVecMulti(vector<vector<vector<int> > > &EPVM) { epiPosVecMulti = EPVM; }
void SimParam::setEpiMagVecMulti(vector<vector<double> > &EMVM) { epiMagVecMulti = EMVM; }

//Growth
void SimParam::setSelG(Sequence s) { selG = s; }
void SimParam::setSelMagVecG(vector<double> & SMV) { selMagVecG = SMV; }
void SimParam::setSelMultiG(vector<Sequence> & SM) { selMultiG = SM; }
void SimParam::setSelMagVecMultiG(vector<vector<double> > & SMVM) { selMagVecMultiG = SMVM; }
void SimParam::setEpiPosVecG(vector<vector<int> > &EPV) { epiPosVecG = EPV; }
void SimParam::setEpiMagVecG(vector<double> &EMV) { epiMagVecG = EMV; }
void SimParam::setEpiPosVecMultiG(vector<vector<vector<int> > > &EPVM) { epiPosVecMultiG = EPVM; }
void SimParam::setEpiMagVecMultiG(vector<vector<double> > &EMVM) { epiMagVecMultiG = EMVM; }



/*
 * Getters
 */
int SimParam::getNa() { return Na; }
int SimParam::getNb() { return Nb; }
int SimParam::getNumLoci() { return numLoci;}
int SimParam::getNt() { return Nt; }
int SimParam::getNumGenes() { return numGenes; }
int SimParam::getNumGenesWithData() { return numGenesWithData; }
unsigned long int SimParam::getSeed() { return seed; }
ofstream* SimParam::getOutputFile() { return outputFile; }
int SimParam::getNumReps() { return numReps; }
int SimParam::getNumLociPerGene() { return numLociPerGene; }
bool SimParam::getWithinHostSelectionPresent() { return withinHostSelectionPresent; }
bool SimParam::getUseExternalReadStatistic() { return useExternalReadStatistic; }
int SimParam::getGeneLength() { return geneLength; }
int SimParam::getFilterHaplotypesMethod() { return filterHaplotypesMethod; }
bool SimParam::getRemoveMonomorphicSim() { return removeMonomorphicSim; }

//Transmission
char SimParam::getSel(int index) { return sel.getBase(index); }
Sequence* SimParam::getSel() { return &sel; }
vector<double>* SimParam::getSelMagVec() { return &selMagVec; }
vector<Sequence>* SimParam::getSeqVec() { return &seqVec; }
vector<vector<int> >* SimParam::getEpiPosVec() { return &epiPosVec; }
vector<double>* SimParam::getEpiMagVec() { return &epiMagVec; }

vector<Sequence>* SimParam::getSelMulti() {
    if(selMulti.size()>0) {
        return &selMulti;
    } else {
        cout << "ERROR: selMulti is not initialised. Are you sure you are using multiple genes?\n";
        return &selMulti;
    }
}


vector<vector<Sequence> >* SimParam::getSeqVecMulti() {
    if(seqVecMulti.size()>0) {
        return &seqVecMulti;
    } else {
        cout << "ERROR: selSeqVecMulti is not initialised.\n";
        return &seqVecMulti;
    }
    
}

vector<vector<double> >* SimParam::getSelMagVecMulti() {
    
    if(selMagVecMulti.size()>0) {
        return &selMagVecMulti;
    } else {
        cout << "ERROR: selMagVecMulti is not initialised.\n";
        return &selMagVecMulti;
    }
    
}

vector<vector<vector<int> > >* SimParam::getEpiPosVecMulti() {

	if(epiPosVecMulti.size() == 0) {
		cout << "ERROR: epiPosVecMulti is not initialised.\n";
	}
	return &epiPosVecMulti;
}

vector<vector<double> >* SimParam::getEpiMagVecMulti() {

	if(epiMagVecMulti.size() == 0) {
		cout << "ERROR: epiMagVecMulti is not initialised.\n";
	}
	return &epiMagVecMulti;
}

//Growth
Sequence* SimParam::getSelG() { return &selG; }
vector<double>* SimParam::getSelMagVecG() { return &selMagVecG; }
vector<vector<int> >* SimParam::getEpiPosVecG() { return &epiPosVecG; }
vector<double>* SimParam::getEpiMagVecG() { return &epiMagVecG; }


vector<Sequence>* SimParam::getSelMultiG() {
    if(selMultiG.size()>0) {
        return &selMultiG;
    } else {
        cout << "ERROR: selMultiG is not initialised. Are you sure you are using multiple genes?\n";
        return &selMultiG;
    }
}

vector<vector<double> >* SimParam::getSelMagVecMultiG() {
    
    if(selMagVecMultiG.size()>0) {
        return &selMagVecMultiG;
    } else {
        cout << "ERROR: selMagVecMultiG is not initialised.\n";
        return &selMagVecMultiG;
    }
    
}

vector<vector<vector<int> > >* SimParam::getEpiPosVecMultiG() {

	if(epiPosVecMultiG.size() == 0) {
		cout << "ERROR: epiPosVecMultiG is not initialised.\n";
	}
	return &epiPosVecMultiG;
}

vector<vector<double> >* SimParam::getEpiMagVecMultiG() {

	if(epiMagVecMultiG.size() == 0) {
		cout << "ERROR: epiMagVecMultiG is not initialised.\n";
	}
	return &epiMagVecMultiG;
}
