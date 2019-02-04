//
//  pathPHMG.hpp
//  FluTransmissionProject
//
//  Created by Casper Lumby on 11/12/2015.
//  Copyright © 2015 Casper Lumby. All rights reserved.
//

#ifndef pathPHMG_hpp
#define pathPHMG_hpp


#include "path.h"
#include "pathPH.h"
#include <vector>

class PathPHMG : public Path {
    
private:
    std::vector<PathPH> PHpaths;
	bool physicalPosMatter, importHaps, importHapsMahan, filterData;
	std::string importHapsFreq;
	std::string importHapsMahanPath; //Path to folder containing sub folders for each gene
	std::string importHapsMahanFilename; //Filename for input of Mahan haplotypes and frequencies
	int minReadDepth;
    
    
public:
    
    PathPHMG(); //Declare default constructor
    ~PathPHMG(); //Deconstructor
    
    void addPHPath(PathPH& p);
    PathPH getPHPath(int index);
    int getNumOfPHPaths();
	bool getPhysicalPosMatter();
	void setPhysicalPosMatterToTrue();
	void setImportHapsFreq(std::string& freq);
	void setImportHapsMahan(std::string& path);
	void setImportHapsMahanFilename(std::string& filename);
	void setMinReadDepth(int minRD);
	int getMinReadDepth();
	void setFilterData(bool fd);
	bool getFilterData();
    
    void loadH5N1FilePaths(std::string &folder);
	void loadFluFilePaths(std::string &folder);
   
};

#endif /* pathPHMG_hpp */

