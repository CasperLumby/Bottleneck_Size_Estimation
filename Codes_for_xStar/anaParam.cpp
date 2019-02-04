
//Forward declared dependencies

//Included dependencies
#include "anaParam.h"
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <sys/stat.h>
#include "misc.h"


using namespace std;

//Constructors
AnaParam::AnaParam() : c(200), maxNt(10000), seed(1), useGreedyMMS(false),  withinHostSelectionFolder(""), BICpenalty(10), selectionCapPresent(false), useIntegralApproach(false), maxNumFittedParams(10), terminateEarly(0), filterHaplotypesMethod(0), noVar(false), meanOnly(false), growthFactor(22), analyseSL(false) {
    //Do nothing
}



AnaParam::~AnaParam() { }; //Deconstructor


//Getters
int AnaParam::getMaxNt() { return maxNt; }
unsigned long int AnaParam::getSeed() { return seed; }
void AnaParam::setSeed(int s) { seed = s; }
vector<Sequence>* AnaParam::getSelVec() { return &selVec; }
vector<vector<double> >* AnaParam::getSelMagVec() { return &selMagVec; }
double AnaParam::getC() { return c; }
void AnaParam::setC(double C) { c = C; }
vector<double>* AnaParam::getLogFactStore() { return logFactStore; }
int AnaParam::getNumPrePoints() { return numPrePoints; }


void AnaParam::loadOutputFolder(string & fullPathToFolder) { pathToFolder=fullPathToFolder; }
void AnaParam::loadUseGreedyMMS(bool value) { useGreedyMMS = value; }

string AnaParam::getOutputFolder() {
    if(pathToFolder.length() > 0) { return pathToFolder; }
    else {
        cout << "ERROR: Path to output folder not loaded. Expect segmentation fault";
        return pathToFolder;
    }
}



bool AnaParam::getUseGreedyMMS() { return useGreedyMMS; }

void AnaParam::setUseSharedBottleneck(bool value) { useSharedBottleneck = value; }
bool AnaParam::getUseSharedBottleneck() { return useSharedBottleneck; }


void AnaParam::loadWithinHostSelectionFolder(string& fullPathToFolder) { withinHostSelectionFolder = fullPathToFolder; }
string AnaParam::getWithinHostSelectionFolder() { return withinHostSelectionFolder; }

//Don't think these are in use anymore. Within-host selection is read in by data file (dPHMG or dPHMR).
void AnaParam::loadWithinHostSelection(string & fullPathToFolder, vector<vector<int> >& physicalPos) { 

	selGvec.clear();
	selGmagVec.clear();

	vector<string> genes {"HA", "MP", "NA", "NP", "NS", "PA", "PB1", "PB2"};

	//Add within host selection coefficients for each gene in turn
	for(unsigned int i=0; i<genes.size(); i++) {

		bool selFoundGene = false;
		loadWithinHostSelectionGene(fullPathToFolder, genes[i], physicalPos[i], &selFoundGene);


		if(selFoundGene == true) { //Selection present for current gene

			withinHostSelectionPresent.push_back(true);
		} else {
			withinHostSelectionPresent.push_back(false);
		}
	} 

}

//Get information about within host selection from file
//As some data may have been ignored due to low counts or monorphic variants, the within host selection information
//may vary from the data in dataPHMG. As such, we use physical position information to ensure we only keep the selection
//coefficients relevant for the data we have.
void AnaParam::loadWithinHostSelectionGeneOld(string & fullPathToFolder, string & gene, vector<int>& physicalPos) { 

	//File path
	stringstream ssSelFile;
	ssSelFile << fullPathToFolder << "/" << gene << "/Selection.out";	
	string selFileString = ssSelFile.str();

	//Selection information containers
	Sequence selG;
	vector<double> selGmag;
	vector<Model::epistasis> epiG;
	vector<double> epiGmag;
	
	if(fileExists(selFileString)) {

		ifstream selFile;
		selFile.open(selFileString.c_str());


		for(string line; getline(selFile, line); ) {

			vector<string> words = split(line, ' '); //Defined in misc
			
			if(words[0].compare("Epi") != 0) { //Not an epistatic coefficient

				int physPosCurrent = atoi(words[0].c_str());
				bool physPosFound = false;

				for(unsigned int i=0; i<physicalPos.size(); i++) {
				
					if(physPosCurrent == physicalPos[i]) {
						
						physPosFound = true; break;
					}
				}

				if(physPosFound == true) { //Current loci is in dataset
		
					if(words.size() > 1) { //This loci has selection

						//The following assumes that the physical pos in data is a subset of the physical pos in the within-host file, i.e. there are
						//no physical pos in the dataset which isn't also in the within-host file.
						//Furthermore, it assumes that physical positions in the within-host file are in an increasing order and that epi cofficients
						//follow after all selection coefficients have been established.
						//In short, this allows us to use push_back methods.
						selG.addBase(words[1][0]); //This is a string of length 1 (e.g. "A" or "C"), hence the [0] to get the char value ('A' or 'C')
						selGmag.push_back(atof(words[2].c_str()));
					
					} else { //No selection at this loci

						selG.addBase('-');
						selGmag.push_back(0);
					}

				}
				
			} else { //Epistatic coefficient

				Model::epistasis epi;
				vector<int> positions;
				bool epistasisValid = true;

				//Loop over epistasis positions
				for(unsigned int i=1; i<(words.size()-1); i++) { //E.g. words = {"Epi", "741" "1036" "2.59084"} so only {1,2} are positions


					int physPosCurrent = atoi(words[i].c_str());

					//Find the locus numder corresponding to this physical position
					int locusCurrent = -1;
					for(unsigned int j=0; j<physicalPos.size(); j++) {

						if(physPosCurrent == physicalPos[j]) { locusCurrent = j; break; }
					}

					//Add locus to positions. Again, assumes loci are ordered by size
					if(locusCurrent != -1) {
						positions.push_back(locusCurrent);
					} else{
						epistasisValid = false; break;
					}
				}

				//If all epistasis positions are in the set of physical pos in the dataset, then add it to the epi model. Else ignore
				if(epistasisValid == true) {
					epi.setPositions(positions);
					epi.setDim(positions.size());

	
					epiG.push_back(epi);
					epiGmag.push_back(atof(words[words.size()-1].c_str()));
				}
				
			}
		}
	
	}

	selGvec.push_back(selG);
	selGmagVec.push_back(selGmag);
	epiGvec.push_back(epiG);
	epiGmagVec.push_back(epiGmag);
}




//Get information about within host selection from file
//As some data may have been ignored due to low counts or monorphic variants, the within host selection information
//may vary from the data in dataPHMG. As such, we use physical position information to ensure we only keep the selection
//coefficients relevant for the data we have.
void AnaParam::loadWithinHostSelectionGene(string & fullPathToFolder, string & gene, vector<int>& physicalPos, bool* selFound) { 

	//File path
	stringstream ssSelFile;
	ssSelFile << fullPathToFolder << "/" << gene << "/Sel_params.out";	
	string selFileString = ssSelFile.str();

	//Selection information containers
	Sequence selG;
	vector<double> selGmag;
	vector<Model::epistasis> epiG;
	vector<double> epiGmag;
	
	if(fileExists(selFileString)) {

		cout << "Within host selection file exists for gene " << gene << ".\n";
		cout << "Printing true physical positions: "; printIntVector(physicalPos);

		ifstream selFile;
		selFile.open(selFileString.c_str());
		int physPosCounter = 0; //Counts the position in physicalPos were are at


		for(string line; getline(selFile, line); ) {

			cout << "Printing current line: " << line << "\n";
			vector<string> words = split(line, ' '); //Defined in misc
			
			if(words[0][0] != 'E') { //Not an epistatic coefficient


				//Need to get the physical position
				//Either words is of the form number OR charNumberChar coefficient, e.g. 223 or G400A -0.2374
				int physPosCurrent = -1;
				if(words.size() == 1) { //words[0] is a position
					physPosCurrent = atoi(words[0].c_str());
				
				} else if(words.size() == 2) { //Words[0] of the form charNumberChar

					string physPosCurrentString = words[0].substr(1,words[0].length()-2);
					physPosCurrent = atoi(physPosCurrentString.c_str());	

				} else {

					cout << "Error in reading in within host selection information from file (anaPara.cpp). Exiting.\n";
					exit(1);
				}

				bool physPosFound = false;

				for(unsigned int i=0; i<physicalPos.size(); i++) {
				
					if(physPosCurrent == physicalPos[i]) {
						
						physPosFound = true; break;
						cout << "physical position found!\n";
					}
				}

				if(physPosFound == true) { //Current loci is in dataset


					//Some physical pos is not found in the within host selection file
					//i.e. sometimes physPosCurrent will be ahead of the postion in 
					//physicalPos. As such, advance selG and selGmag with neutral data.
					while(physPosCurrent > physicalPos[physPosCounter]) {

						selG.addBase('-');
						selGmag.push_back(0);
						physPosCounter++;
					}					

		
					if(words.size() > 1) { //This loci has selection

						//The following assumes that the physical pos in data is a subset of the physical pos in the within-host file, i.e. there are
						//no physical pos in the dataset which isn't also in the within-host file.
						//Furthermore, it assumes that physical positions in the within-host file are in an increasing order and that epi cofficients
						//follow after all selection coefficients have been established.
						//In short, this allows us to use push_back methods.
						selG.addBase(words[0][words[0].length()-1]); //Getting the last char in word
						selGmag.push_back(atof(words[1].c_str()));
						*selFound = true;
					
					} else { //No selection at this loci

						selG.addBase('-');
						selGmag.push_back(0);
					}
					
					physPosCounter++;

				}

//				cout << "Priting selG for gene " << gene << " at intermediate stage: ";
//				selG.print();
//				printDoubleVector(selGmag);
				
			} else { //Epistatic coefficient


				//Assumes that epistasis acts on same alleles as selection does
				//E.g. if selection acts on alleles T and A at positions 2 and 4,
				//then epi most also act on these alleles (and not e.g. alleles C and A
				//at position 2 and 4). This seems to be the case with the output from
				//Chris' code.

				Model::epistasis epi;
				vector<int> positions;
				bool epistasisValid = true;

				//Loop over epistasis positions
				for(unsigned int i=1; i<(words.size()-1); i++) { //E.g. words = {"Epi", "868T", "1263A", "-0.376633"} so only {1,2} are positions


					int physPosCurrent = atoi(words[i].c_str()); //Ignores the allele at the end (see atoi documentation)
//					cout << "physPosCurrent: " << physPosCurrent << "\n";

					//Find the locus number corresponding to this physical position
					int locusCurrent = -1;
					for(unsigned int j=0; j<physicalPos.size(); j++) {

						if(physPosCurrent == physicalPos[j]) { locusCurrent = j; break; }
					}

					//Add locus to positions. Again, assumes loci are ordered by size
					if(locusCurrent != -1) {
						positions.push_back(locusCurrent);
					} else{
						epistasisValid = false; break;
					}
				}

				//If all epistasis positions are in the set of physical pos in the dataset, then add it to the epi model. Else ignore
				if(epistasisValid == true) {

					cout << "Epistasis valid.\n";
					epi.setPositions(positions);
					epi.setDim(positions.size());

	
					epiG.push_back(epi);
					epiGmag.push_back(atof(words[words.size()-1].c_str()));
				}
				
			}
		}
	
		//Some physical pos is not found in the within host selection file
		//i.e. sometimes selG and selGmag needs to have neutral data added
		//to the end after going through the within host selection file.
		while(physPosCounter < (int) physicalPos.size()) {

			selG.addBase('-');
			selGmag.push_back(0);
			physPosCounter++;
		}					

	} else { //No within host selection file for current gene

		cout << "Within host selection file doesn't exists for gene " << gene << ".\n";
		cout << "Printing true physical positions: "; printIntVector(physicalPos);
		//It is possible that whilst there is no within host selection present for the current gene,
		//there can still be data for this gene. As such, we need to create a neutral within host
		//selection model for such genes.
		for(unsigned int i=0; i<physicalPos.size(); i++) {

			selG.addBase('-');
			selGmag.push_back(0);
		}
		
	
	}

	
	cout << "Priting final selG for gene " << gene << ": ";
	selG.print();
	printDoubleVector(selGmag);
				

	selGvec.push_back(selG);
	selGmagVec.push_back(selGmag);
	epiGvec.push_back(epiG);
	epiGmagVec.push_back(epiGmag);
}

bool AnaParam::getWithinHostSelectionPresent(int gene) { 

	if(withinHostSelectionPresent.size() > 0) { //Within host selection defined
		return withinHostSelectionPresent[gene]; 
	} else {
		return false; //Within host selection not defined, so return default of no selG for this gene
	}
}
bool AnaParam::getWithinHostSelectionPresent() {
	
	for(unsigned int i=0; i<withinHostSelectionPresent.size(); i++) {

		if(withinHostSelectionPresent[i] == true) { return true; }
	}

	return false;
 }
void AnaParam::setWithinHostSelectionPresent(int gene, bool b) { withinHostSelectionPresent[gene] = b; }
void AnaParam::addWithinHostSelectionPresent(bool b) { withinHostSelectionPresent.push_back(b); }

vector<Sequence> AnaParam::getSelGvec() { return selGvec; }
void AnaParam::setSelGvec(vector<Sequence>& SGV) { selGvec = SGV; }
vector<vector<double> > AnaParam::getSelGmagVec() { return selGmagVec; }
void AnaParam::setSelGmagVec(vector<vector<double> >& SGMV) { selGmagVec = SGMV; }
vector<vector<Model::epistasis> > AnaParam::getEpiGvec() { return epiGvec; }
void AnaParam::setEpiGvec(std::vector<std::vector<Model::epistasis> >& EGV) { epiGvec = EGV; }
vector<vector<double> > AnaParam::getEpiGmagVec() { return epiGmagVec; }
void AnaParam::setEpiGmagVec(std::vector<std::vector<double> >& EGMV) { epiGmagVec = EGMV; }



void AnaParam::setHapFitG(vector<vector<double> > & HFG) { hapFitG = HFG; }
void AnaParam::setHapFitG(vector<vector<vector<double> > > & HFG) { hapFitGmultiRep = HFG; }
vector<double> AnaParam::getHapFitG(int gene) { return hapFitG[gene]; }
vector<double> AnaParam::getHapFitG(int rep, int gene) { return hapFitGmultiRep[rep][gene]; }

void AnaParam::setBICpenalty(double p) { BICpenalty = p; }
double AnaParam::getBICpenalty() { return BICpenalty; }
void AnaParam::setSelectionCap(double value) { selectionCap = value; selectionCapPresent = true; }
double AnaParam::getSelectionCap() { return selectionCap; }
bool AnaParam::getSelectionCapPresent() { return selectionCapPresent; }
bool AnaParam::getUseIntegralApproach() { return useIntegralApproach; }
void AnaParam::setUseIntegralApproach(bool b) { useIntegralApproach = b; }
void AnaParam::setMaxNumFittedParams(int num) { maxNumFittedParams = num; }
int AnaParam::getMaxNumFittedParams() { return maxNumFittedParams; }
void AnaParam::setMaxNt(int mNt) { maxNt = mNt; }
void AnaParam::setTerminateEarly(int t) { terminateEarly = t; } //E.g. 1==After hap inference, 2==After freq inference, 3==...
int AnaParam::getTerminateEarly() { return terminateEarly; }
void AnaParam::setFilterHaplotypesMethod(int fhm) { filterHaplotypesMethod = fhm; }
int AnaParam::getFilterHaplotypesMethod() { return filterHaplotypesMethod; }
void AnaParam::setNoVar(bool b) { noVar = b; }
bool AnaParam::getNoVar() { return noVar; }
void AnaParam::setMeanOnly(bool b) { meanOnly = b; }
bool AnaParam::getMeanOnly() { return meanOnly; }
void AnaParam::setGrowthFactor(int gf) { growthFactor = gf; }
int AnaParam::getGrowthFactor() { return growthFactor; }
void AnaParam::setAnalyseSL(bool b) { analyseSL = b; }
bool AnaParam::getAnalyseSL() { return analyseSL; }
