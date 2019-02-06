//Include guard
#ifndef PATHPH_H
#define PATHPH_H


//Forward declared dependencies

//Included dependencies
#include <string>
#include "path.h"


class PathPH : public Path {
private:
    std::vector<std::string> hapPaths;
    std::vector<std::string > contribPaths;
    std::vector<std::string > lociPaths;
	std::string singleLocusTrajFilePath, multiLocusTrajFilePath;
	bool physicalPosMatter;
	std::string physicalPosPath;
	bool importHaps;
	bool importHapsMahan;
	int minReadDepth;
	bool filterData;
	std::string importHapsFreq;
	std::string importHapsFile;
	std::string importHapsMahanFile; //Filepath to Mahan's file containing haplotypes and frequencies
	std::string timesFile;
    
public:
    //	Path(std::string& p);
    PathPH(); //Declare default constructor
    ~PathPH(); //Deconstructor
    
    std::vector<std::string> getHapPaths();
    std::vector<std::string> getContribPaths();
    std::vector<std::string> getLociPaths();
	std::string getMultiLocusTrajFilePath();
	std::string getSingleLocusTrajFilePath();
	bool getPhysicalPosMatter();
	void setPhysicalPosMatterToTrue();
	std::string getPhysicalPosPath();
	void setMinReadDepth(int minRD);
	int getMinReadDepth();
	void setFilterData(bool fd);
	bool getFilterData();
	void setImportHapsFreq(std::string& freq);
	std::string getImportHapsFile();
	void setImportHapsMahanFile(std::string& path);
	std::string getImportHapsMahanFile();
	bool getImportHaps();
	bool getImportHapsMahan();
	std::string getTimesFile();
	
    
    //H5N1 data
    void loadH5N1HAGeneFilePaths(std::string& folder);
    void loadH5N1NAGeneFilePaths(std::string& folder);
    void loadH5N1M1GeneFilePaths(std::string& folder);


	//Flu data
	void loadFluGeneFilePaths(std::string &folder, std::string &gene);
    
};

#endif 
