//
//  dataPH.hpp
//  FluTransmissionProject
//
//  Created by Casper Lumby on 16/11/2015.
//  Copyright © 2015 Casper Lumby. All rights reserved.
//

#ifndef dataPH_hpp
#define dataPH_hpp

//Included dependencies
#include "data.h"
#include "pathPH.h"
#include "sequence.h"
#include "diploidSequence.h"
#include "simParamPH.h"
#include <vector>
#include <iostream>


class DataPH : public Data {
	friend class DataSL; //DataSL is a friend of DataPH so it can more easily access partial hap information

private:
	//Internal variables
	PathPH* path; //Path to data
	std::vector<int> physicalPos; //Only relevant in some cases
	SimParamPH *spPH;
	std::vector<Sequence> importedFullHaps;
	bool importFullHaps, importHapsMahan;
	std::vector<double> importedqB, importedqA;
	int deltaDays; //DeltaDays describe the difference between before and after, e.g. if before = day 3 and after = day 5, then deltaDays = 2. Default is 1.


	//True full haplotype variables (only existing if data are simulated)
	std::vector<Sequence> fullHapsTrue;
	std::vector<double> qBtrue;


	struct slTraj {
	private:
		int NbMinor, NbMajor, NaMinor, NaMajor, NbTot, NaTot, pos;
		char major, minor;

	public:
		slTraj(); //Constructor initialises varaibles to zero
		void setNbMinor(int N);
		void setNbMajor(int N);
		void setNaMinor(int N);
		void setNaMajor(int N);
		int getNbMinor();
		int getNbMajor();
		int getNaMinor();
		int getNaMajor();
		void incrementNbMinor(int x = 1);
		void incrementNbMajor(int x = 1);
		void incrementNaMinor(int x = 1);
		void incrementNaMajor(int x = 1);
		void countNbTot();
		void countNaTot();
		int getNbTot();
		int getNaTot();
		void setPos(int p);
		int getPos();
		void setMinor(char c);
		void setMajor(char c);
		char getMinor();
		char getMajor();
		void print();
	};
	std::vector<slTraj> slTrajs;
	
   
    struct partialHap {
    private:
        Sequence seq;
        int Nb;
        int Na;
        std::vector<int> contribs;
        
    public:
        void setSeq(Sequence &s);
        void setNb(int NB);
        void setNa(int NA);
        void addContrib(int c);
        Sequence getSeq();
        int getNb();
        int getNa();
        std::vector<int>* getContribs();
        void clearContributions();
        void print();
		void writeToFile(std::ofstream &outputFile);
		void writeToFileMahan(std::ofstream &outputFile, std::vector<int>& physPosPH);
	void removeLocus(int locus);
	int getNumberOfSites();
    };
    
    struct setOfPartialHaps {
    private:
        std::vector<partialHap> phVec;
        int NbTot;
        int NaTot;
	std::vector<int> lociCovered; //E.g. if --XXX- then loci covered 2,3,4
        
    public:
        void setNb(int index, int Nb);
        void setNa(int index, int Na);
        std::vector<int> getNb();
        std::vector<int> getNa();
        void setNbTot(int N);
        void setNaTot(int N);
        void countNbTot();
        void countNaTot();
        int getNbTot();
        int getNaTot();
        void addPartialHap(partialHap &ph);
	partialHap getPartialHap(int phIndex);
        Sequence getSeq(int phIndex);
        int getNumOfPartialHaps();
        void addContrib(int index, int contrib);
        std::vector<int>* getContribs(int index);
        void clearContributions();
        void print();
        std::vector<Sequence> getSequences();
        std::vector<double> getPHFreqs(std::vector<double> &fhFreqs);
        std::vector<double> getPHhapFit(std::vector<double> &fhFreqs, std::vector<double> &fhHapFit);
        void removePH(int phIndex);
        int getNb(int phIndex);
        int getNa(int phIndex);
	void writeToFile(std::ofstream &outputFile);
	void writeToFileMahan(std::ofstream &outputFile, std::vector<int>& physicalPos);
	void writeToFileHaps(std::string &filePathGeneFolder, std::vector<int>& physicalPos, int sophIndex);
	void removeLocus(int locus);
	std::vector<int> getPosCovered(std::vector<int>& allPos);
	void computeDegreeOfOverlapBefore(int& weightedTotal, int& totalNumberOfReads);
	void computeDegreeOfOverlapAfter(int& weightedTotal, int& totalNumberOfReads);
	void setLociCovered(std::vector<int>& lc);
	std::vector<int> getLociCovered();
	int getNumberOfSites();
    };
    
    struct phDataSingleSim { // Data container for sinlge simulation
    private:
        std::vector<setOfPartialHaps> dataset;
        
    public:
        void addSetOfPartialHaps(setOfPartialHaps& setPH);
        void clear();
	void clearObsOnly();
        
        setOfPartialHaps * getSetOfPartialHaps(int index);
        void print();
        void setNb(int index1, int index2, int Nb);
        std::vector<Sequence> getSequences();
        std::vector<int> getBeforeCounts(int sphIndex);
        std::vector<int> getAfterCounts(int sphIndex);
        int numOfPHSets();int getNumDatasets(); // returns number of datasets //THESE TWO DO TEH SAME! REMOVE ONE OF THEM
        int getNumOfPartialHaps(int datasetIndex);
        void clearContributions();
	void computeContribs(std::vector<Sequence> &fhs);
        std::vector<int>* getContribs(int sphIndex, int phIndex);
        std::vector<int> getNb(int sphIndex);
        std::vector<int> getNa(int sphIndex);
        std::vector<double> getPHFreqs(int sphIndex, std::vector<double> &fhFreqs);
        std::vector<double> getPHhapFit(int sphIndex, std::vector<double> &fhFreqs, std::vector<double> &fhHapFit);
	std::vector<std::vector<Sequence> > findMMS(std::vector<Sequence> & MSfhs, int seed);
	int getNbTot();
	int getNaTot();
	bool decomposeAndCheck(std::vector<Sequence> &fhs);
	void writeToFile(std::ofstream &outputFile);
	void writeToFileMahan(std::ofstream &outputFile, std::vector<int>& physicalPos);
	void writeToFileHaps(std::string &filePathGeneFolder, std::vector<int>& physicalPos);
	void removeLocus(int locus);
	void removeSOPH(int index);
	void mergeHaps();
	void orderHaps();
	double computeDegreeOfOverlapBefore();
	double computeDegreeOfOverlapAfter();
	int getNumberOfSites();
	void ensurePHsInImportedHaps(std::vector<Sequence>& fhs);
    } data;
   
	struct pairedEndRead {
	private:
		int start1, end1, start2, end2, hapID;

	public:
		void generateRead(double meanReadLength, double stDevReadLength, double meanGapLength, double stDevGapLength, int length, gsl_rng* r, double C, std::vector<double>& freqs);
		bool posCovered(int& pos);
		bool posCovered(std::vector<int>& pos);
		void print();
		int getHapID();
	};
    
public:
    DataPH(); //Constructor
    ~DataPH(); //Deconstructor
    
    void readData(Path * p);
    void simulateData(SimParam *sp);
    void simulateDataRealistic(SimParam *sp, std::vector<int>& physicalPos, DiploidSequence& sCDS, int geneLength); //sCDS = semiConstrainedDiploidSequence

	//Decompose functions are static as used for MMS methods on arbitrary phdss objects.
	//In future, if MMS methods removed, remove static tag
    static void decomposeFullHaps(std::vector<Sequence> & fullHaps, phDataSingleSim& phdss, std::vector<Sequence> &truePH); //Need to pass phdss explicitly as static function
    static void decomposeFullHapsFixedSize(std::vector<Sequence> & fullHaps, phDataSingleSim& phdss, int size, bool useExternalReadStatistic);
    static void decomposeFullHapsLevel(std::vector<Sequence> & fullHaps, phDataSingleSim& phdss, int level, bool useExternalReadStatistic);
    void drawObservationsBefore(int Nb, std::vector<double> freqs, double C, gsl_rng *r);
    void drawObservationsAfter(int Na, std::vector<double> freqs, double C, gsl_rng *r);
	void drawObservationsBeforeRealisticWrong(int depth, int length, std::vector<int>& physicalPos, double meanReadLength, double stDevReadLength, double meanGapLength, double stDevGapLength, std::vector<Sequence>& haps, std::vector<double>& freqs, double C, gsl_rng* r);
	void drawObservationsBeforeRealistic(int depth, int length, std::vector<int>& physicalPos, double meanReadLength, double stDevReadLength, double meanGapLength, double stDevGapLength, std::vector<Sequence>& haps, std::vector<double>& freqs, double C, gsl_rng* r);
	void drawObservationsAfterRealisticWrong(int depth, int length, std::vector<int>& physicalPos, double meanReadLength, double stDevReadLength, double meanGapLength, double stDevGapLength, std::vector<Sequence>& haps, std::vector<double>& freqs, double C, gsl_rng* r);
	void drawObservationsAfterRealistic(int depth, int length, std::vector<int>& physicalPos, double meanReadLength, double stDevReadLength, double meanGapLength, double stDevGapLength, std::vector<Sequence>& haps, std::vector<double>& freqs, double C, gsl_rng* r);
	std::vector<Sequence> getSequences();
	void removeLowCountObservations();
	void removeLowCountSamfire();
	void applyObsFrequencyCutOffAfter(double cutOff);
	void removeMonomorphicSites();
	void removeLowReadDepth(std::vector<int> posWithSufficientReadDepth);
    void computeContribs(std::vector<Sequence> &fhs);
    std::vector<int> getBeforeCounts(int sphIndex);
    std::vector<int> getAfterCounts(int sphIndex);
    int numOfPHSets();
    void print();
    void clearContributions();
	void clear(); //Clears everything
    int getNumOfDatasets();
    int getNumOfPartialHaps(int datasetIndex);
    std::vector<int> * getContribs(int sphIndex, int phIndex);
    std::vector<int> getNb(int sphIndex);
    std::vector<int> getNa(int sphIndex);
    std::vector<double> getPHFreqs(int sphIndex, std::vector<double> & fhFreqs);
    std::vector<double> getPHhapFit(int sphIndex, std::vector<double> & fhFreqs, std::vector<double> & fhHapFit);
    void importPHs();
	std::vector<std::vector<Sequence> > findMMS(std::vector<Sequence> & MSfhs, int seed);
	int getNbTot();
	int getNaTot();	
	int getNtot();
	bool decomposeAndCheck(std::vector<Sequence> &fhs);
	void writeToFile(std::ofstream &outputFile);
	void writeToFileMahan(std::ofstream &outputFile, std::vector<int>& physicalPos);
	void writeToFileHaps(std::string &filePathGeneFolder, std::vector<int>& physicalPos);
	std::vector<int> getPhysicalPos();
	void clearPhysicalPos();
	std::vector<Sequence>& getImportedFullHaps();
	void setImportedFullHaps(std::vector<Sequence>& ifh);
	void updateDataWithRespectToImportedFullHaps();
	bool getImportFullHaps();
	int getDeltaDays();

	//Full haplotype data -- for debugging use only - requires simulation of data
	std::vector<Sequence> getFullHapsTrue();
	std::vector<double> getqBtrue();


	std::vector<double> getDonorSLFreqs();
	std::vector<double> getRecipientSLFreqs();
	std::vector<int> getRecipientSLObs();
	std::vector<int> getRecipientSLObsTot();
};

#endif /* dataPH_hpp */
