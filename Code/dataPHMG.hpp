//Include guard
#ifndef DATAPHMG_H
#define DATAPHMG_H


//Forward declared dependencies


//Included dependencies
#include "data.h"
#include "pathPHMG.hpp"
#include "dataPH.hpp"
#include "diploidSequence.h"
#include "model.hpp"
#include <vector>
#include <string>


//Full haplotype multi gene
class DataPHMG : public Data {
private:
	//Internal variables
	PathPHMG path; //Path to data
	std::vector<std::vector<int> > physicalPos; //Only relevant in some cases
	std::vector<bool> genePresent; //Describes if data is present for each gene (e.g. dim=8 for flu)
    
	struct phmgData{ //Data container for multiple genes
	private:
		std::vector<DataPH> genes;
	public:
		void addGene(DataPH & g);
		void clearData();
		int getNumGenes();
		DataPH* getGene(int index);
		int getDeltaDays();
        
	} data;


public:
	struct WHSel { //Data container for within host selection
	private:
		Sequence selG;
		std::vector<double> selGmag;
		std::vector<Model::epistasis> epiG;
		std::vector<double> epiGmag;
		bool selPresent;
	public:
		void setSelG(Sequence& SELG);
		void setSelGmag(std::vector<double>& SELGMAG);
		void setEpiG(std::vector<Model::epistasis>& EPIG);
		void setEpiGmag(std::vector<double>& EPIGMAG);
		void setSelPresent(bool& b);
		bool getSelPresent();
		Sequence getSelG();
		std::vector<double> getSelGmag();
		std::vector<Model::epistasis> getEpiG();
		std::vector<double> getEpiGmag();
		
	};

	struct WHSelMG { //Multi gene version
	private:
		std::vector<WHSel> whVec;
	public:
		void addWHSel(WHSel& WHSEL);
		WHSel getWHSel(int gene);
		void clear();
		bool getSelPresent();
		bool getSelPresent(int gene);

	} whSelMG;
    
    
public:
	DataPHMG(); //Constructor
	~DataPHMG(); //Deconstructor
    
	void readData(Path * p);
	void simulateData(SimParam *sp);
	void simulateDataRealistic(SimParam *sp); //For use by the dataPHMG only
	void simulateDataRealistic(SimParam* sp, std::vector<std::vector<int> >& physicalPos, std::vector<DiploidSequence>& sCDS); //For use by dataPHMR and dataPHMG
	void applyObsFrequencyCutOffAfter(double cutOff);
	int getNumGenes();
	DataPH* getGene(int index);
	std::vector<std::vector<int> > getPhysicalPos();
	bool getGenePresent(int g);
	void setGenePresent(int g, bool b);
	bool getDataPresent();
	void loadWithinHostSelection(std::string& fullPathToFolder);
	WHSel loadWithinHostSelectionGene(std::string& fullPathToFolder, int geneIndex, std::string& gene);
	WHSelMG getWHSelMG();
	bool getWHSelPresent(int gene);
	bool physicalPosPresent();
	void updatePhysicalPos(int gene, DataPH* dPH);
	int getDeltaDays();
    
};

#endif
