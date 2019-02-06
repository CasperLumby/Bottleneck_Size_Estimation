#include "dataPHMG.hpp"
#include <stdio.h>
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "misc.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace std;

//Class related methods + constructors
DataPHMG::DataPHMG() { } //Constructor does nothing
DataPHMG::~DataPHMG() {} //Deconstructor

void DataPHMG::readData(Path * p) {
    
    path = * static_cast<PathPHMG *>(p);
    
    //Add data gene by gene
    for(int i=0; i<path.getNumOfPHPaths(); i++) {

	DataPH dph;
        PathPH pph = path.getPHPath(i);
        dph.readData(&pph);

	if(dph.getNumOfDatasets() > 0) {
		cout << "Data added for gene " << i << "\n";
	}
        

	if(path.getPhysicalPosMatter()==true) { //store physical positions for all genes

		vector<int> physicalPosCurrentGene = dph.getPhysicalPos();
		physicalPos.push_back(physicalPosCurrentGene);
		
	}
        
	//Check if gene contains any data
	if(dph.getNumOfDatasets() > 0) {
		genePresent.push_back(true);
	} else {

		cout << "Gene data not present for gene: " << i << "\n";
		genePresent.push_back(false);
	}
        	
	//Add the gene    
	data.addGene(dph);
	//cout << i << " added: ";
 
    }
    
}

void DataPHMG::simulateData(SimParam *sp) {
    
	cout << "Simulating full haplotype data for multiple genes...\n";
	SimParamPH* spPH = dynamic_cast<SimParamPH*>(sp);
	int dim = spPH->getDim();
	if(dim<3) {
		cout << "Error, the dimension must be above 2. Aborting simulation.";
		return;
	}
	int Na = spPH->getNa();
	int Nb = spPH->getNb();
	double C = spPH->getC();
	int Nt = spPH->getNt();
	int numGenes= spPH->getNumGenes(); //Default 8
	int numGenesWithData = spPH->getNumGenesWithData(); //Number of genes with actual data, e.g.3
	int seed = spPH->getSeed();

	//Transmission
	vector<Sequence> selMulti = *(spPH->getSelMulti());
	vector<vector<double> > selMagVecMulti = *(spPH->getSelMagVecMulti());
	vector<vector<vector<int> > > epiPosVecMulti = *(spPH->getEpiPosVecMulti());
	vector<vector<double> > epiMagVecMulti = *(spPH->getEpiMagVecMulti());



	//Set up RNG
	gsl_rng * rng;
	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	rng = gsl_rng_alloc (T);
	gsl_rng_set(rng, seed);

	//Reset data container if no external read information inputted
	if(spPH->getUseExternalReadStatistic() == false) {
		data.clearData();
	}
    
	for(int i=0; i<numGenes; i++) { //Loop over numGenes genes of which only numGenesWithData genes contain data
        
		bool geneContainsData = false;
		if(spPH->getUseExternalReadStatistic() == false) {
		
			if(i<numGenesWithData) { //Gene contains data
				geneContainsData = true;
			}
		} else {

			if(genePresent[i] == true) {
				geneContainsData = true;
			}
		}
		
		if(geneContainsData == true) { //If gene i contains actual data, simulate it
			cout << "Simulating full haplotype data for gene #" << i << ".\n";
			DataPH* g;
			DataPH dPH; //Empty DataPH object
			
			//If using external read information, read data has already been added to dataPHMG
			//object data. Hence g must be associated with this read information.
			if(spPH->getUseExternalReadStatistic() == true) {
				
				g = data.getGene(i);
			
			} else {

				//Need to instantiate an empty DataPH object and point g to it
//				DataPH dPH; //Does this go out of scope outside the else { }? It seems to work, but maybe definte it outside to be sure. Look into this!
				g = &dPH;
			}

			int numLoci = selMulti[i].getLength();
			cout << "Num loci for gene " << i << ": " << numLoci << "\n";
	
			//Set up output folder for simulated data (for possible subsequent manual analysis)
			ofstream outputFile;
			string filePathFolder = spPH->getPathToFolder();	
			if(spPH->getRepScenario() == true) {
				string replicateID = spPH->getRepID();
				stringstream ss;
				ss << filePathFolder << "SimulatedDataRep_" << replicateID << "_Gene_" << i << ".dat";
				outputFile.open((ss.str()).c_str()); 
			
			} else { //Non-replicate scenario

				stringstream ss;
				ss << filePathFolder << "SimulatedData_Gene_" << i << ".dat";
				outputFile.open((ss.str()).c_str());
			}        


			SimParamPH spSingleGene; //Could in theory copy spPH
			spSingleGene.setGeneIndex(i);
			spSingleGene.setPathToFolder(filePathFolder);
			spSingleGene.setFormat(spPH->getFormat()); //E.g. "Mahan"
			spSingleGene.setNb(Nb);
			spSingleGene.setNa(Na);
			spSingleGene.setC(C);
			spSingleGene.setNumLoci(numLoci);
			spSingleGene.setNt(Nt);
			spSingleGene.setSeed(gsl_rng_uniform_int(rng,1000000)); //In theory, should use full integer range, but this should suffice. See https://www.gnu.org/software/gsl/manual/html_node/Sampling-from-a-random-number-generator.html for more info.

			//Transmission
			spSingleGene.setSel(selMulti[i]);
			spSingleGene.setSelMagVec(selMagVecMulti[i]);

			//Not a perfect solution, but fine for now
			if(epiPosVecMulti.size() > 0 ) { //If epistasis has been initiated
				spSingleGene.setEpiPosVec(epiPosVecMulti[i]);
				spSingleGene.setEpiMagVec(epiMagVecMulti[i]);
			}

			//Growth
			if(spPH->getWithinHostSelectionPresent() == true) {
			
				//Multi gene parameters
				vector<Sequence> selMultiG = *(spPH->getSelMultiG());
				vector<vector<double> > selMagVecMultiG = *(spPH->getSelMagVecMultiG());
				vector<vector<vector<int> > > epiPosVecMultiG = *(spPH->getEpiPosVecMultiG());
				vector<vector<double> > epiMagVecMultiG = *(spPH->getEpiMagVecMultiG());
	
				//Single gene parameters
				spSingleGene.setWithinHostSelectionPresent(true);
				spSingleGene.setSelG(selMultiG[i]);
				spSingleGene.setSelMagVecG(selMagVecMultiG[i]);

				//Not a perfect solution, but fine for now
				if(epiPosVecMultiG.size() > 0 ) { //If epistasis has been initiated
					spSingleGene.setEpiPosVecG(epiPosVecMultiG[i]);
					spSingleGene.setEpiMagVecG(epiMagVecMultiG[i]);
				}
			}

			spSingleGene.setDim(dim);
			spSingleGene.setNumGenes(1);
			spSingleGene.setOutputFile(&outputFile);
			if(spPH->getUseExternalReadStatistic() == true) {
				spSingleGene.setUseExternalReadStatistic(true);
			}
			if(spPH->getFilterHaplotypesMethod() > 0) {
				cout << "Using filter method " << spPH->getFilterHaplotypesMethod() << "\n";
				spSingleGene.setFilterHaplotypesMethod(spPH->getFilterHaplotypesMethod());
			}
			spSingleGene.setRemoveMonomorphicSim(spPH->getRemoveMonomorphicSim());
			cout << "Before simulating data in dataPH\n" << flush;
			g->simulateData(&spSingleGene);
			outputFile.close();

			//Gene is already added if using external read information
			if(spPH->getUseExternalReadStatistic() == false) {
				data.addGene(*g);

				//Check if any data in gene after simulation
				//Update: I think this is no longer necessary as dataPH::simulateData ensures
				//that the resulting dataset is as requested (i.e. covers a specific set of loci)
				if(g->getNumOfDatasets() > 0) {
					genePresent.push_back(true); //Indicate that gene contains data
				} else {
					genePresent.push_back(false);
				}
			}

		
	
		} else { //No data available for said gene
			DataPH g;

			//Gene is already added if using external read information
			if(spPH->getUseExternalReadStatistic() == false) {
				data.addGene(g); //Add empty gene
				genePresent.push_back(false); //Indicate that gene is empty
			}


		}
	}
}

//To use by dataPHMG only. Generates physical pos and sCDS to be used in simulateDataRealistic(sp,physPos,sCDS) method
void DataPHMG::simulateDataRealistic(SimParam* sp) {

	//Setup rng
	 gsl_rng * rng;
	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	rng = gsl_rng_alloc (T);
	gsl_rng_set(rng, sp->getSeed());


	//Generate physical positions
	vector<Sequence> selMulti = *(sp->getSelMulti());
//	vector<string> genes {"HA", "MP", "NA", "NP", "NS", "PA", "PB1", "PB2"};
	vector<int> geneLengths;
	int geneLength = sp->getGeneLength();
	if(geneLength != -1) { //Gene length defined externally

		//All genes share the same gene length
		geneLengths = {geneLength, geneLength, geneLength, geneLength, geneLength, geneLength, geneLength, geneLength};

	} else { //Use standard flu gene lengths
        	geneLengths = {1778,1027,1413,1565,890,2233,2341,2341}; //http://www.rapidreferenceinfluenza.com/chapter/B978-0-7234-3433-7.50009-8/aim/influenza-virus-structure
	}
	

	vector<vector<int > > physPos;

	for(unsigned int g=0; g<selMulti.size(); g++) {

		int maxPos = geneLengths[g];
		vector<int> physPosGene;

		for(int i=0; i<selMulti[g].getLength(); i++) { //Loop over loci in gene g

			int newLocus = -1;
			bool uniqueLocus = false;
			while(uniqueLocus == false) {
				
				uniqueLocus = true;
				newLocus = gsl_rng_uniform_int(rng, maxPos) +1; // Generate between 1 and N
				for(unsigned int j=0; j<physPosGene.size(); j++) {

					if(newLocus == physPosGene[j]) { uniqueLocus = false; break; }
				}

			}

			physPosGene.push_back(newLocus);
		}

		sort(physPosGene.begin(), physPosGene.end());
		physPos.push_back(physPosGene);
	}

	//Generate semi-constrained diploid sequence
	vector<DiploidSequence> sCDS;
	for(unsigned int g=0; g<selMulti.size(); g++) {
	
		Sequence minorEmpty = Sequence(selMulti[g].getLength());
		DiploidSequence sCDSgene = DiploidSequence(selMulti[g], minorEmpty);
		sCDS.push_back(sCDSgene);
	}

	//Run simulateDataRealistic(sp, physPos, sCDS)
	sp->setSeed(gsl_rng_uniform_int(rng,1000000)); //Generate new seed for future use
	simulateDataRealistic(sp, physPos, sCDS);
}
void DataPHMG::simulateDataRealistic(SimParam* sp, vector<vector<int> >& physicalPos, vector<DiploidSequence>& sCDS) { //sCDS = semiConstrainedDiploidSequence
 
	cout << "Simulating (realistic) full haplotype data for multiple genes...\n";
	SimParamPH* spPH = dynamic_cast<SimParamPH*>(sp);
	int dim = spPH->getDim();
	if(dim<3) { //Could in theory have dim==2...
		cout << "Error, the dimension must be above 2. Aborting simulation.";
		return;
	}
	double C = spPH->getC();
	int Nt = spPH->getNt();
	int numGenes= spPH->getNumGenes(); //Number of genes, default 8
	int numGenesWithData = spPH->getNumGenesWithData(); //Number of genes with actual data, e.g.3 
	int seed = spPH->getSeed();
	vector<Sequence> selMulti = *(spPH->getSelMulti()); //To be replaced by sCDS
	vector<vector<double> > selMagVecMulti = *(spPH->getSelMagVecMulti());
	vector<vector<vector<int> > > epiPosVecMulti = *(spPH->getEpiPosVecMulti());
	vector<vector<double> > epiMagVecMulti = *(spPH->getEpiMagVecMulti());

//	vector<string> genes {"HA", "MP", "NA", "NP", "NS", "PA", "PB1", "PB2"};
	vector<int> geneLengths;
	int geneLength = sp->getGeneLength();
	if(geneLength != -1) { //Gene length defined externally

		//All genes share the same gene length
		geneLengths = {geneLength, geneLength, geneLength, geneLength, geneLength, geneLength, geneLength, geneLength};

	} else { //Use standard flu gene lengths
        	geneLengths = {1778,1027,1413,1565,890,2233,2341,2341}; //http://www.rapidreferenceinfluenza.com/chapter/B978-0-7234-3433-7.50009-8/aim/influenza-virus-structure
	}
	



	//Set up RNG
	gsl_rng * rng;
	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	rng = gsl_rng_alloc (T);
	gsl_rng_set(rng, seed);

	//Reset data container
	data.clearData();
    
	for(int i=0; i<numGenes; i++) { //Loop over numGenes genes of which only numGenesWithData genes contain data
        
		if(i<numGenesWithData) { //If gene i contains actual data, simulate it
			cout << "Simulating full haplotype data for gene #" << i << ".\n";
			DataPH dPH;

			int numLoci = sCDS[i].getLength();
			cout << "Num loci for gene " << i << ": " << numLoci << "\n";
	
			//Set up output folder for simulated data (for possible subsequent manual analysis)
			ofstream outputFile;
			string filePathFolder = spPH->getPathToFolder();	
			if(spPH->getRepScenario() == true) {
				string replicateID = spPH->getRepID();
				stringstream ss;
				ss << filePathFolder << "SimulatedDataRep_" << replicateID << "_Gene_" << i << ".dat";
				outputFile.open((ss.str()).c_str()); 
			
			} else { //Non-replicate scenario

				stringstream ss;
				ss << filePathFolder << "SimulatedData_Gene_" << i << ".dat";
				outputFile.open((ss.str()).c_str());
			}        

			SimParamPH spSingleGene; //Could in theory copy spPH
			spSingleGene.setGeneIndex(i);
			spSingleGene.setPathToFolder(filePathFolder);
			spSingleGene.setC(C);
			spSingleGene.setNumLoci(numLoci);
			spSingleGene.setNt(Nt);
			spSingleGene.setFormat(spPH->getFormat());
			spSingleGene.setGrowthFactor(spPH->getGrowthFactor());
			spSingleGene.setSeed(gsl_rng_uniform_int(rng,1000000)); //In theory, should use full integer range, but this should suffice. See https://www.gnu.org/software/gsl/manual/html_node/Sampling-from-a-random-number-generator.html for more info.

			//Transmission
			spSingleGene.setSel(selMulti[i]);
			spSingleGene.setSelMagVec(selMagVecMulti[i]);

			//Not a perfect solution, but fine for now
			if(epiPosVecMulti.size() > 0 ) { //If epistasis has been initiated
				spSingleGene.setEpiPosVec(epiPosVecMulti[i]);
				spSingleGene.setEpiMagVec(epiMagVecMulti[i]);
			}

			//Growth
			if(spPH->getWithinHostSelectionPresent() == true) {
			
				//Multi gene parameters
				vector<Sequence> selMultiG = *(spPH->getSelMultiG());
				vector<vector<double> > selMagVecMultiG = *(spPH->getSelMagVecMultiG());
				vector<vector<vector<int> > > epiPosVecMultiG = *(spPH->getEpiPosVecMultiG());
				vector<vector<double> > epiMagVecMultiG = *(spPH->getEpiMagVecMultiG());
	
				//Single gene parameters
				spSingleGene.setWithinHostSelectionPresent(true);
				spSingleGene.setSelG(selMultiG[i]);
				spSingleGene.setSelMagVecG(selMagVecMultiG[i]);

				//Not a perfect solution, but fine for now
				if(epiPosVecMultiG.size() > 0 ) { //If epistasis has been initiated
					spSingleGene.setEpiPosVecG(epiPosVecMultiG[i]);
					spSingleGene.setEpiMagVecG(epiMagVecMultiG[i]);
				}
			}

			spSingleGene.setDim(dim);
			spSingleGene.setNumGenes(1);
			spSingleGene.setOutputFile(&outputFile);
			spSingleGene.setRemoveMonomorphicSim(spPH->getRemoveMonomorphicSim());

			dPH.simulateDataRealistic(&spSingleGene, physicalPos[i], sCDS[i], geneLengths[i]);
			outputFile.close();
			data.addGene(dPH);

			//Check if any data in gene after simulation
			if(dPH.getNumOfDatasets() > 0) { 
				genePresent.push_back(true); //Indicate that gene contains data
			} else {
				genePresent.push_back(false); //Indicate that gene doesn't contains data
			}
	
		} else { //No data available for said gene
			cout << "No data available for gene #" << i << ".\n";
			DataPH dPH;
			data.addGene(dPH); //Add empty gene
			genePresent.push_back(false); //Indicate that gene is empty
		}
	}

}

//Enforce a cutoff such that any after observations with a frequency less than this cutoff are set to zero.
void DataPHMG::applyObsFrequencyCutOffAfter(double cutOff) {

	if(cutOff >= 0.5) {
		cout << "Applied after observations frequency cut-off cannot be >= 0.5. Exiting.\n";
		exit(1);
	}

	for(int i=0; i<data.getNumGenes(); i++) { //Loop over numGenes genes of which only numGenesWithData genes contain data

		cout << "Applying a frequency cut-off of " << cutOff << " for gene " << i << ":\n";
		data.getGene(i)->applyObsFrequencyCutOffAfter(cutOff);
	}
}

int DataPHMG::getNumGenes() {	return data.getNumGenes(); }
DataPH* DataPHMG::getGene(int index) { return data.getGene(index); }

//Methods related to structs
void DataPHMG::phmgData::addGene(DataPH& g) { genes.push_back(g); }
void DataPHMG::phmgData::clearData() { genes.clear(); }
int DataPHMG::phmgData::getNumGenes() { return (int) genes.size(); }
DataPH* DataPHMG::phmgData::getGene(int index) { return &(genes[index]); }
vector<vector<int> > DataPHMG::getPhysicalPos() { return physicalPos; }
bool DataPHMG::getGenePresent(int g) { return genePresent[g]; }
void DataPHMG::setGenePresent(int g, bool b) { genePresent[g] = b; }
bool DataPHMG::getDataPresent() {

	for(unsigned int g=0; g<genePresent.size(); g++) {

		if(genePresent[g] == true) { return true; }
	}
	return false;
}


void DataPHMG::WHSel::setSelG(Sequence& SELG) { selG = SELG; }
void DataPHMG::WHSel::setSelGmag(vector<double>& SELGMAG) { selGmag = SELGMAG; }
void DataPHMG::WHSel::setEpiG(vector<Model::epistasis>& EPIG) { epiG = EPIG; }
void DataPHMG::WHSel::setEpiGmag(vector<double>& EPIGMAG) { epiGmag = EPIGMAG; }
void DataPHMG::WHSel::setSelPresent(bool& b) { selPresent = b; }
bool DataPHMG::WHSel::getSelPresent() { return selPresent; }
Sequence DataPHMG::WHSel::getSelG() { return selG; }
vector<double> DataPHMG::WHSel::getSelGmag() { return selGmag; }
vector<Model::epistasis> DataPHMG::WHSel::getEpiG() { return epiG; }
vector<double> DataPHMG::WHSel::getEpiGmag() { return epiGmag; }

void DataPHMG::WHSelMG::addWHSel(WHSel& WHSEL) { whVec.push_back(WHSEL); }
void DataPHMG::WHSelMG::clear() { whVec.clear(); }
bool DataPHMG::WHSelMG::getSelPresent() { 

	for(unsigned int g=0; g<whVec.size(); g++) {

		if(whVec[g].getSelPresent() == true) { return true; }
	}

	return false;
}
bool DataPHMG::WHSelMG::getSelPresent(int gene) {

	if((int) whVec.size() >= gene+1) {

		return whVec[gene].getSelPresent();

	} else { return false; }

 }
DataPHMG::WHSel DataPHMG::WHSelMG::getWHSel(int gene) { return whVec[gene]; }
DataPHMG::WHSelMG DataPHMG::getWHSelMG() { return whSelMG; }
bool DataPHMG::getWHSelPresent(int gene) { return  whSelMG.getSelPresent(gene); }


void DataPHMG::loadWithinHostSelection(string& fullPathToFolder) {

	whSelMG.clear(); //Remove any previously loaded within host selection data

	vector<string> genes {"HA", "MP", "NA", "NP", "NS", "PA", "PB1", "PB2"};

	//Add within host selection for each gene in turn
	for(unsigned int i=0; i<genes.size(); i++) {

		WHSel withinHostSel = loadWithinHostSelectionGene(fullPathToFolder, i, genes[i]);

		whSelMG.addWHSel(withinHostSel);
	}

}



//Get information about within host selection from file
//As some data may have been ignored due to low counts or monorphic variants, the within host selection information
//may vary from the data in dataPHMG. As such, we use physical position information to ensure we only keep the selection
//coefficients relevant for the data we have.
DataPHMG::WHSel DataPHMG::loadWithinHostSelectionGene(string & fullPathToFolder, int geneIndex, string & gene) { 

	//File path
	stringstream ssSelFile;
	ssSelFile << fullPathToFolder << "/" << gene << "/Sel_params.out";	
	string selFileString = ssSelFile.str();

	//Selection information containers
	Sequence selG;
	vector<double> selGmag;
	vector<Model::epistasis> epiG;
	vector<double> epiGmag;
	bool selFound = false;
	
	WHSel withinHostSel;
	vector<int> physPosGene = physicalPos[geneIndex];
	
	if(fileExists(selFileString)) {

		cout << "Within host selection file exists for gene " << gene << ".\n";
		cout << "Printing true physical positions: "; printIntVector(physPosGene);

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

					cout << "Error in reading in within host selection information from file. Exiting.\n";
					exit(1);
				}

				bool physPosFound = false;

				for(unsigned int i=0; i<physPosGene.size(); i++) {
				
					if(physPosCurrent == physPosGene[i]) {
						
						physPosFound = true; break;
						cout << "physical position found!\n";
					}
				}

				if(physPosFound == true) { //Current loci is in dataset


					//Some physical pos is not found in the within host selection file
					//i.e. sometimes physPosCurrent will be ahead of the postion in 
					//physicalPos. As such, advance selG and selGmag with neutral data.
					while(physPosCurrent > physPosGene[physPosCounter]) {

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
			
///// We are now limiting selection at the level of haplotype fitness, not selection itself
//
//						double strength = atof(words[1].c_str());
//						if(strength > 10) {
//							cout << "Selection strength for " << physPosCurrent << " is " << strength << " so limiting to 10.\n";
//							strength = 10;
//						} else if(strength < -10) {
//							cout << "Selection strength for " << physPosCurrent << " is " << strength << " so limiting to -10.\n";
//							strength = -10;
//						}
//						selGmag.push_back(strength);

						selGmag.push_back(atof(words[1].c_str()));
						selG.addBase(words[0][words[0].length()-1]); //Getting the last char in word
						selFound = true;
					
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
					for(unsigned int j=0; j<physPosGene.size(); j++) {

						if(physPosCurrent == physPosGene[j]) { locusCurrent = j; break; }
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

///// We are now limiting selection at the level of haplotype fitness, not selection itself
//					double strength = atof(words[words.size()-1].c_str());
//					if(strength > 10) {
//						strength = 10;
//					} else if(strength < -10) {
//						cout << "Epistasis strength for " << printIntVectorToString(positions) << " is " << strength << " so limiting to -10.\n";
//						strength = -10;
//					}
//					epiGmag.push_back(strength);

					epiGmag.push_back(atof(words[words.size()-1].c_str()));
				}
				
			}
		}
	
		//Some physical pos is not found in the within host selection file
		//i.e. sometimes selG and selGmag needs to have neutral data added
		//to the end after going through the within host selection file.
		while(physPosCounter < (int) physPosGene.size()) {

			selG.addBase('-');
			selGmag.push_back(0);
			physPosCounter++;
		}					

	} else { //No within host selection file for current gene

		cout << "Within host selection file doesn't exists for gene " << gene << ".\n";
		cout << "Printing true physical positions: "; printIntVector(physPosGene);
		//It is possible that whilst there is no within host selection present for the current gene,
		//there can still be data for this gene. As such, we need to create a neutral within host
		//selection model for such genes.
		for(unsigned int i=0; i<physPosGene.size(); i++) {

			selG.addBase('-');
			selGmag.push_back(0);
		}
		
	
	}

	
	cout << "Priting final selG for gene " << gene << ": ";
	selG.print();
	printDoubleVector(selGmag);
				

	withinHostSel.setSelG(selG);
	withinHostSel.setSelGmag(selGmag);
	withinHostSel.setEpiG(epiG);
	withinHostSel.setEpiGmag(epiGmag);
	withinHostSel.setSelPresent(selFound);

	return withinHostSel;
}

bool DataPHMG::physicalPosPresent() { 

	if(physicalPos.size() > 0) { return true;
 	} else {
		return false;
	}

}

void DataPHMG::updatePhysicalPos(int gene, DataPH* dPH) {

	physicalPos[gene] = dPH->getPhysicalPos();

}

int DataPHMG::getDeltaDays() { return data.getDeltaDays(); }
int DataPHMG::phmgData::getDeltaDays() { return genes[0].getDeltaDays(); } //Assumes at least one gene is loaded!!
