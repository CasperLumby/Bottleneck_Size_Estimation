//
//  pathPHMR.hpp
//  FluTransmissionProject
//
//  Created by Casper Lumby on 11/12/2015.
//  Copyright © 2015 Casper Lumby. All rights reserved.
//

#ifndef pathPHMR_hpp
#define pathPHMR_hpp


#include "path.h"
#include "pathPHMG.hpp"
#include <vector>

class PathPHMR : public Path {
    
private:
	std::vector<PathPHMG> PHMGpaths;
	bool importHaps;
	std::string importHapsFreq;
    
    
public:
    
	PathPHMR(); //Declare default constructor
	~PathPHMR(); //Deconstructor
    
	void addPHMGPath(PathPHMG& p);
	PathPHMG getPHMGPath(int index);
	int getNumOfPHMGPaths();
	void setImportHapsFreq(std::string& freq);
    
	void loadFluReplicateFilePaths(std::vector<std::string> & repFolders);
   
};

#endif /* pathPHMR_hpp */

