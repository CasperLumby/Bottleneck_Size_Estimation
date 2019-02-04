//Include guard
#ifndef DATAPHMR_H
#define DATAPHMR_H


//Forward declared dependencies


//Included dependencies
#include "data.h"
#include "pathPHMR.hpp"
#include "dataPHMG.hpp"
#include <vector>


//Full haplotype multi gene
class DataPHMR : public Data {
private:
    //Internal variables
    PathPHMR path; //Path to data
	std::vector<std::vector<std::vector<int> > > physicalPos;
    
    struct phmrData{ //Data container for multiple replicates
    private:
        std::vector<DataPHMG> replicates;
    public:
        void addReplicate(DataPHMG & r);
        void clearData();
        int getNumReplicates();
        DataPHMG* getReplicate(int index);
	int getDeltaDays(int repIndex);
        
        
    } data;
    
public:

	struct WHSelMR { //Multi rep version
	private:
		std::vector<DataPHMG::WHSelMG> whMGvec;

	public:
		void addWHSelMG(DataPHMG::WHSelMG& WHSELMG);
		DataPHMG::WHSelMG getWHSelMG(int rep);
		DataPHMG::WHSel getWHSel(int rep, int gene);
		void clear();
		bool getSelPresent();
		bool getSelPresent(int rep, int gene);

	} whSelMR;

    
public:
    DataPHMR(); //Constructor
    ~DataPHMR(); //Deconstructor
    
    void readData(Path * p);
    void simulateData(SimParam *sp);
    void simulateRealData(SimParam *sp);
    int getNumReplicates();
    DataPHMG* getReplicate(int index);

	std::vector<std::vector<std::vector<int> > > getPhysicalPos();
	std::vector<std::vector<int> > findAllPos(); //Of form: <gene<pos> >
	std::vector<std::vector<std::vector<int> > > findMapFromAllPosToReps(std::vector<std::vector<int> >& allPos);
	std::vector<std::vector<bool> > findMapFromAllPosToShared(std::vector<std::vector<int> >& allPos, std::vector<std::vector<std::vector<int> > >& mapFromAllPosToReps);
    
	bool getWHSelPresent(int rep, int gene);
	WHSelMR getWHSelMR();
	void loadWithinHostSelection(std::string& fullPathToFolder);
	bool physicalPosPresent();
	void updatePhysicalPos(int rep, DataPHMG* dPHMG);
	int getDeltaDays(int repIndex);
};

#endif
