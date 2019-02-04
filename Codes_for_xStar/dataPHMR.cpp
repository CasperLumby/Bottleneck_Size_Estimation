#include "dataPHMR.hpp"
#include <stdio.h>
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "misc.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>

using namespace std;

//Class related methods + constructors
DataPHMR::DataPHMR() { } //Constructor does nothing
DataPHMR::~DataPHMR() {} //Deconstructor

void DataPHMR::readData(Path * p) {
    
	path = * static_cast<PathPHMR *>(p);
    
	//Add data rep by rep
	for(int i=0; i<path.getNumOfPHMGPaths(); i++) {
		cout << "Adding data for replicate " << i << "\n";
		DataPHMG dphmg;
		PathPHMG pphmg = path.getPHMGPath(i);
		dphmg.readData(&pphmg);

		//Store physical positions
		vector<vector<int> > physicalPosCurrentReplicate = dphmg.getPhysicalPos();
		physicalPos.push_back(physicalPosCurrentReplicate);
        
		data.addReplicate(dphmg);
		cout << "Data for replicate " << i << " added: ";
	}
    
}


//Everything in here is hardcoded for simplicity.
//Considers the case where all variants are shared between all replicates.
void DataPHMR::simulateData(SimParam *sp) {

	cout << "Simulating data for multiple replicates.\n";

	//Set up RNG
	gsl_rng * rng;
	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	rng = gsl_rng_alloc (T);
	gsl_rng_set(rng, sp->getSeed());

	int numReps = sp->getNumReps(); //Number of replicates
	int numGenes= sp->getNumGenes(); //Default 8
	int numGenesWithData = sp->getNumGenesWithData(); //Number of genes with actual data, e.g.3
	int numLociPerGene = sp->getNumLociPerGene();

	//Create vector containing number of loci per gene, e.g. {5,5,5,0,0,0,0,0}
	vector<int> lociPerGenePerRep;
	for(int i=0; i<numGenes; i++) {
		if(i<numGenesWithData) { lociPerGenePerRep.push_back(numLociPerGene); }
		else { lociPerGenePerRep.push_back(0); }
	}
	
	//Gene segment lengths
//	vector<string> genes {"HA", "MP", "NA", "NP", "NS", "PA", "PB1", "PB2"};
	vector<int> geneLengths = {1778,1027,1413,1565,890,2233,2341,2341}; //http://www.rapidreferenceinfluenza.com/chapter/B978-0-7234-3433-7.50009-8/aim/influenza-virus-structure
	
	vector<vector<int> > allPos; //Order: gene, loci
	for(int g=0; g<numGenes; g++) { //Loop over genes

		int maxPos = geneLengths[g];
		vector<int> allPosGene;

		if(g==0) { //set up physicalPos during first round
	
			for(int r=0; r<numReps; r++) {
				vector<vector<int> > empty1;
				physicalPos.push_back(empty1);

				for(int g2=0; g2<numGenes; g2++) {
					vector<int> empty2;
					physicalPos[r].push_back(empty2);
				}
			}
		}

		//Find the shared loci (they are all shared in this scenario)
		for(int i=0; i<lociPerGenePerRep[g]; i++) {

			int newLocus = -1;
			bool uniqueLocus = false;

			while(uniqueLocus == false) {

				uniqueLocus = true;
				newLocus = gsl_rng_uniform_int(rng, maxPos) +1; // Generate between 1 and N
				for(unsigned int j=0; j<allPosGene.size(); j++) {
					if(newLocus == allPosGene[j]) { uniqueLocus = false; break; }
				}
			}

			//Add to all pos
			allPosGene.push_back(newLocus);
	
			//Add loci to replicates (identical to allPosGene as all loci shared)		
			for(int r=0; r<numReps; r++) {
				//Add to physical pos of current replicate
				physicalPos[r][g].push_back(newLocus);
			}
		}

		allPos.push_back(allPosGene);
	}


	//Order loci by position (ascending order)
	for(int g=0; g<numGenes; g++) {
		for(int r=0; r<numReps; r++) {

			sort(physicalPos[r][g].begin(), physicalPos[r][g].end());
		}

		sort(allPos[g].begin(), allPos[g].end());	
	}

	//Set up sequencess for each rep and gene
	//E.g. sequences of the form --AT-, -A-T-, AT--- for rep 0,1,2 for gene 0 	
	//WLOG define, for now, shared sequences to have type A/C
	//Can always find a neater way of inputting this later
	//Therefore sequences need to be diploid and not haploid as earlier
	vector<vector<DiploidSequence> > sCDS; //semiConstrainedDiploidSequence
	vector<Sequence> selMulti = *(sp->getSelMulti());
	for(int r=0; r<numReps; r++) {
		
		vector<DiploidSequence> empty;
		sCDS.push_back(empty);

		for(int g=0; g<numGenes; g++) {

			//Create empty sCDS for gene g - to be filled in below
			DiploidSequence sCDSgene;
			if(r==0) { //First time around
				DiploidSequence sCDSgeneEmpty = DiploidSequence(physicalPos[r][g].size()); //Diploid sequence with length = #loci in gene
				sCDSgene = generateRandomSemiConstrainedDiploidSequence(sCDSgeneEmpty, rng); //Defined in misc
			
			} else { //All other reps just need to be the same as rep 0 in this scenario (all loci shared)

				sCDSgene = sCDS[0][g];
			}


			//Correct for fixed alleles/selection
			for(unsigned int i=0; i<physicalPos[r][g].size(); i++) { //Loop over loci and update sCDSgene as appropriate

				//Check if site under selection
				if(selMulti[g].getBase(i) != '-') { //Site under selection

					sCDSgene.setMajor(i,selMulti[g].getBase(i));
					while(sCDSgene.getMajor().getBase(i) == sCDSgene.getMinor().getBase(i)) {

						int temp = (int) gsl_rng_uniform_int(rng,4);
						if(temp==0) {
							sCDSgene.setMinor(i,'A');
						} else if (temp==1) {
							sCDSgene.setMinor(i,'C');
						} else if (temp==2) {
							sCDSgene.setMinor(i,'G');
						} else {
							sCDSgene.setMinor(i,'T');
						}
					}
				}
			}
	
			cout << "Printing physical pos for replicate " << r << " gene " << g << ":\n";
			printIntVector(physicalPos[r][g]);
			cout << "Printing semi constrained diploid sequence for replicate " << r << " gene " << g << ":\n";
			sCDSgene.print();
			sCDS[r].push_back(sCDSgene);
		}
	}

	//So we now how defined how the various replicates relate to each other
	//The next step is to geneate the data for each replicate in turn
	for(int r=0; r<numReps; r++) {

		SimParamPH* spRep = dynamic_cast<SimParamPH*>(sp); //Do we need a SimParamPHMG? 
		spRep->setSeed(gsl_rng_uniform_int(rng,1000000)); //In theory, should use full integer range, but this should suffice. See https://www.gnu.org/software/gsl/manual/html_node/Sampling-from-a-random-number-generator.html for more info.
		spRep->setRepScenario(true); //Indicate to downstream classes that this is a replicate scenario
		spRep->setRepID(to_string(r));
		DataPHMG dPHMG;
		dPHMG.simulateDataRealistic(spRep, physicalPos[r], sCDS[r]); //Pass more information later, e.g. mean and var read length

		data.addReplicate(dPHMG);
	}
    
}

//This is for simulating data based on real data only
void DataPHMR::simulateRealData(SimParam *sp) {

	//Set up RNG
	gsl_rng * rng;
	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	rng = gsl_rng_alloc (T);
	gsl_rng_set(rng, sp->getSeed());


	vector<vector<int> > allPos = findAllPos();
	vector<vector<vector<int> > > mapFromAllPosToReps = findMapFromAllPosToReps(allPos); //Mapping of form vec<gene,rep,pos> with pos e.g. {-1,0,1,-1,2,3,-1,4} if zeroth index in all pos not in rep, first index in all pos corresponds to zeroth index in rep, etc..
	vector<Sequence> selVec = *sp->getSelMulti();
	vector<vector<double> > selMagVec = *sp->getSelMagVecMulti();
	

	//Run simulated for each replicate
	for(unsigned int r=0; r<mapFromAllPosToReps[0].size(); r++) { //WLOG get numReps from gene 0


		//Convert selVec and selMagVec to individual rep versions
		vector<Sequence> selVecRep;
		vector<vector<double> > selMagVecRep;
	
		for(unsigned int g=0; g<allPos.size(); g++) { //Loop over genes

			Sequence selRep;
			vector<double> selMagRep;
			for(unsigned int i=0; i<allPos[g].size(); i++) { //Loop over sites
			
				if(mapFromAllPosToReps[g][r][i] != -1) { //Pos is int replicate

					selRep.addBase(selVec[g].getBase(i));
					selMagRep.push_back(selMagVec[g][i]);

				}
			}	

			selVecRep.push_back(selRep);
			selMagVecRep.push_back(selMagRep);
		}

		//Simulate data for replicate
		//simParam
		//dataPHMG.simulate
	
		SimParamPH* spRep = dynamic_cast<SimParamPH*>(sp); //Do we need a SimParamPHMG? 
		spRep->setSeed(gsl_rng_uniform_int(rng,1000000)); //In theory, should use full integer range, but this should suffice. See https://www.gnu.org/software/gsl/manual/html_node/Sampling-from-a-random-number-generator.html for more info.
		spRep->setRepScenario(true); //Indicate to downstream classes that this is a replicate scenario
		spRep->setRepID(to_string(r)); //What do I use this for?
		spRep->setSelMulti(selVecRep);
		spRep->setSelMagVecMulti(selMagVecRep);
		
		DataPHMG* dPHMG = data.getReplicate(r);
		dPHMG->simulateData(spRep);
		

	}

}



int DataPHMR::getNumReplicates() {	return data.getNumReplicates(); }
DataPHMG* DataPHMR::getReplicate(int index) { return data.getReplicate(index); }
vector<vector<vector<int> > > DataPHMR::getPhysicalPos() { return physicalPos; }
	
//Accumulating physical positions from multiple replicates into a single
//data container of form <gene<pos> >. Positions are sorted.
vector<vector<int> > DataPHMR::findAllPos() { //Of form: <gene<pos> >

	
	vector<vector<int> > allPos;
	
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


//Mapping of form vec<gene,rep,pos> with pos e.g. {-1,0,1,-1,2,3,-1,4} if zeroth index in all pos not in rep, first index in all pos corresponds to zeroth index in rep, etc..
vector<vector<vector<int> > > DataPHMR::findMapFromAllPosToReps(vector<vector<int> >& allPos) {


	vector<vector<vector<int> > > map;

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

//Of form {gene0={true, false, false, true}, gene1={false,false, true,false},...} where true if allPos[gene][pos] is shared between replicates
vector<vector<bool> > DataPHMR::findMapFromAllPosToShared(vector<vector<int> >&allPos, vector<vector<vector<int> > >& mapFromAllPosToReps) {

	vector<vector<bool> > map; 

	//Loop over genes
	for(unsigned int g=0; g<allPos.size(); g++) {
		
		vector<bool> mapGene;
			
		//Loop over pos
		for(unsigned int i=0; i<allPos[g].size(); i++) {

			bool posFoundInAllReps = true;
			for(unsigned int r=0; r<mapFromAllPosToReps[g].size(); r++) {
		
				if(mapFromAllPosToReps[g][r][i] == -1) { //Pos not found in this replicate

					posFoundInAllReps = false;
					break;
				}
			}


			mapGene.push_back(posFoundInAllReps);	
		}

		map.push_back(mapGene);
	}

	return map;
}


//Methods related to structs
void DataPHMR::phmrData::addReplicate(DataPHMG& g) { replicates.push_back(g); }
void DataPHMR::phmrData::clearData() { replicates.clear(); }
int DataPHMR::phmrData::getNumReplicates() { return (int) replicates.size(); }
DataPHMG* DataPHMR::phmrData::getReplicate(int index) { return &(replicates[index]); }



void DataPHMR::WHSelMR::addWHSelMG(DataPHMG::WHSelMG& WHSELMG) { whMGvec.push_back(WHSELMG); }
void DataPHMR::WHSelMR::clear() { whMGvec.clear(); }
bool DataPHMR::WHSelMR::getSelPresent() { 

	for(unsigned int r=0; r<whMGvec.size(); r++) {

		if(whMGvec[r].getSelPresent() == true) { return true; }
	}

	return false;
}
bool DataPHMR::WHSelMR::getSelPresent(int rep, int gene) { return whMGvec[rep].getSelPresent(gene); }
DataPHMG::WHSelMG DataPHMR::WHSelMR::getWHSelMG(int rep) { return whMGvec[rep]; }
DataPHMG::WHSel DataPHMR::WHSelMR::getWHSel(int rep, int gene) { return whMGvec[rep].getWHSel(gene); }
DataPHMR::WHSelMR DataPHMR::getWHSelMR() { return whSelMR; }
bool DataPHMR::getWHSelPresent(int rep, int gene) { return  whSelMR.getSelPresent(rep,gene); }

void DataPHMR::loadWithinHostSelection(string& fullPathToFolder) {

	whSelMR.clear(); //Remove any previously loaded within host selection data

	for(int r=0; r<getNumReplicates(); r++) { //Loop over reps

		DataPHMG* dataCurrentRep = data.getReplicate(r);
		dataCurrentRep->loadWithinHostSelection(fullPathToFolder);	
	
	}


//	vector<string> genes {"HA", "MP", "NA", "NP", "NS", "PA", "PB1", "PB2"};
//
//	//Add within host selection for each gene in turn
//	for(unsigned int i=0; i<genes.size(); i++) {
//
//		WHSel withinHostSel = loadWithinHostSelectionGene(fullPathToFolder, i, genes[i]);
//
//		whSelMG.addWHSel(withinHostSel);
//	}

}

bool DataPHMR::physicalPosPresent() { 

	if(physicalPos.size() > 0) { return true;
 	} else {
		return false;
	}

}

void DataPHMR::updatePhysicalPos(int rep, DataPHMG* dPHMG) {

	physicalPos[rep] = dPHMG->getPhysicalPos();

}

int DataPHMR::getDeltaDays(int repIndex) { return data.getDeltaDays(repIndex); }
int DataPHMR::phmrData::getDeltaDays(int repIndex) { return replicates[repIndex].getDeltaDays(); }



