//
//  dataPHgen.hpp
//  FluTransmissionProject
//
//

#ifndef dataPHgen_hpp
#define dataPHgen_hpp

//Included dependencies
#include "data.h"
#include "pathPH.h"
#include "sequence.h"
#include "diploidSequence.h"
#include "simParamPH.h"
#include <vector>
#include <iostream>


//Generalised version of DataPH. In future, consider porting to IDE and refactoring DataPH -> DataPHspecial (or something similar) and DataPHgen -> DataPH.
class DataPHgen : public Data {

private:
	//Internal variables
	PathPH* path; //Path to data
	std::vector<int> physicalPos; //Only relevant in some cases
	SimParamPH *spPH;

public:
	struct partialHap {
	private:
		Sequence seq;
		std::vector<int> pos;
		std::vector<int> times;
		std::vector<int> obs;
		std::vector<int> contribs;
	
	public:
		void readTraj(std::string &trajString); 
		Sequence getSeq();
		std::vector<int> getPos();
		int getObs(int obsIndex);
		void setSeq(Sequence &s);
		int getNumOfTimePoints();
		void clearContributions();
		void addContrib(int c);
		std::vector<int>* getContribs();
		void print();

	};
   
//    struct partialHap {
//    private:
//        Sequence seq;
//	std::vector<int> obs;
//        std::vector<int> contribs;
//        
//    public:
//        void setSeq(Sequence &s);
//        void setObs(int index, int obs);
//        void addContrib(int c);
//        Sequence getSeq();
//        int getObs(int index);
//        std::vector<int>* getContribs();
//        void clearContributions();
//        void print();
//		void writeToFile(std::ofstream &outputFile);
//	void removeLocus(int locus);
//	int getNumberOfSites();
//    };
//    
	struct setOfPartialHaps {
	private:
		std::vector<partialHap> phVec;
		std::vector<int> obsTot;
		std::vector<int> pos; //E.g. 120, 337, 455
		std::vector<int> lociCovered; //E.g. if --XXX- then loci covered 2,3,4
        
	public:
//        void setObs(int phIndex, int obsIndex, int obs);
//        std::vector<std::vector<int> > getObs();
//        void setObsTot(int phIndex, int N);
		void countObsTot();
//        std::vector<int> getObsTot();
		void addPartialHap(partialHap &ph);
		std::vector<int> getPos();
		void setPos(std::vector<int> &p);
//	partialHap getPartialHap(int phIndex);
//        Sequence getSeq(int phIndex);
	        int getNumOfPartialHaps();
		std::vector<int> getObs(int timeIndex); //Get vector of observations from all PHs in SOPH at time point timeIndex
		void addContrib(int index, int contrib);
	        std::vector<int>* getContribs(int index);
        	void clearContributions();
        	void print();
		std::vector<Sequence> getSequences();
		int getNumOfTimePoints(); //Assumes the times are identical for all PHs
//        std::vector<double> getPHFreqs(std::vector<double> &fhFreqs);
//        std::vector<double> getPHhapFit(std::vector<double> &fhFreqs, std::vector<double> &fhHapFit);
//        void removePH(int phIndex);
//        std::vector<int> getObs(int phIndex);
//	void writeToFile(std::ofstream &outputFile);
//	void removeLocus(int locus);
//	std::vector<int> getPosCovered(std::vector<int>& allPos);
		void setLociCovered(std::vector<int>& lc);
//	std::vector<int> getLociCovered();
//	int getNumberOfSites();
	};
    
private:
	std::vector<setOfPartialHaps> SOPHs;


//    struct phData { // Data container for set of SOPHs
//    private:
//        std::vector<setOfPartialHaps> dataset;
///        
//    public:
//  void addSetOfPartialHaps(setOfPartialHaps& setPH);
///        void clear();
//	void clearObsOnly();
//        
//        setOfPartialHaps * getSetOfPartialHaps(int index);
//        void print();
///        void setObs(int sophIndex, int phIndex, int obsIndex, int obs);
//        std::vector<Sequence> getSequences();
//        std::vector<std::vector<int> > getObs(int sophIndex);
//        int getNumOfSOPHs();
//        int getNumOfPartialHaps(int sophIndex);
//        void clearContributions();
//	void computeContribs(std::vector<Sequence> &fhs);
//        std::vector<int>* getContribs(int sphIndex, int phIndex);
///        std::vector<double> getPHFreqs(int sphIndex, std::vector<double> &fhFreqs);
//        std::vector<double> getPHhapFit(int sphIndex, std::vector<double> &fhFreqs, std::vector<double> &fhHapFit);
//	std::vector<int> getObsTot(); //Accumulated total observations
//	bool decomposeAndCheck(std::vector<Sequence> &fhs);
//	void writeToFile(std::ofstream &outputFile);
//	void removeLocus(int locus);
//	void removeSOPH(int index);
//	void mergeHaps();
//	void orderHaps();
//	int getNumberOfSites();
//	void ensurePHsInImportedHaps(std::vector<Sequence>& fhs);
//    } data;
      
public:
	DataPHgen(); //Constructor
	~DataPHgen(); //Deconstructor
    
	void readData(Path * p);
	void simulateData(SimParam *sp);
	int getNumOfSOPHs();
	int getNumOfTimePoints(); //Assumes the times are identical for all PHs
	setOfPartialHaps* getSOPH(int SOPHindex);

	std::vector<Sequence> getSequences();
//	void removeLowCountObservations();
//	void removeMonomorphicSites();
	void computeContribs(std::vector<Sequence> &fhs);
//    int getNumOfSOPHs();
    void print();
//    void clearContributions();
//	void clear(); //Clears everything
//    int getSizeOfDataset();
//    int getNumOfPartialHaps(int datasetIndex);
//    std::vector<int> * getContribs(int datasetIndex, int sophIndex);
//    std::vector<std::vector<int> > getObs(int datasetIndex);
//    std::vector<double> getPHFreqs(int sphIndex, std::vector<double> & fhFreqs);
//    std::vector<double> getPHhapFit(int sphIndex, std::vector<double> & fhFreqs, std::vector<double> & fhHapFit);
//    void importPHs();
//	std::vector<int> getObsTot();
//	bool decomposeAndCheck(std::vector<Sequence> &fhs);
//	void writeToFile(std::ofstream &outputFile);
//	std::vector<int> getPhysicalPos();
//	void clearPhysicalPos();
//	void updateDataWithRespectToImportedFullHaps();
//	int getDeltaDays();
//

};

#endif /* dataPH_hpp */
