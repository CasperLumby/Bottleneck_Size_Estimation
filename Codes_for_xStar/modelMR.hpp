//
//  modelMR.hpp
//  FluTransmissionProject
//
//  Created by Casper Lumby on 27/01/2016.
//  Copyright © 2016 Casper Lumby. All rights reserved.
//

#ifndef modelMR_hpp
#define modelMR_hpp

#include <stdio.h>
#include <vector>
#include "model.hpp"
#include "diploidSequence.h"


class ModelMR : public Model {
    
   
private:

	std::vector<int> NtNewVec;
	std::vector<int> NtBestVec;
	std::vector<std::vector<std::vector<double> > > qBmean; //Levels: rep, gene	std::vector<std::vector<gsl_matrix*> > qBvar;
	std::vector<std::vector<gsl_matrix*> > qBvar;
	std::vector<std::vector<std::vector<int> > > map; //Map from overall to individual replicates. Levels: Gene/Replicate/Corresponding_pos_in_individual_rep (e.g. dim is length of overall pos)
	std::vector<std::vector<int> > allPos; //List of loci in combined system (all replicates). Levels: Gene/Pos
	std::vector<DiploidSequence> collapsedFullHaps; //List of alleles for allPos. Levels: Gene

    
public:
    
    //Constructors
    ModelMR();
    ~ModelMR();

	void setMap(std::vector<std::vector<std::vector<int> > > &map);
	void setAllPos(std::vector<std::vector<int> > &allPos);
	void setCollapsedFullHaps(std::vector<DiploidSequence> &collapsedFullHaps);
	
	std::vector<std::vector<std::vector<int> > > getMap();
	std::vector<std::vector<int> > getAllPos();
	std::vector<DiploidSequence> getCollapsedFullHaps();

	int getNtNew(int index);
	int getNtBest(int index);
	std::vector<int> getNtNew();
	std::vector<int> getNtBest();
	void initialiseAndRandomiseCoefs(gsl_rng *r, std::vector<std::vector<std::vector<Sequence> > > &haps);
	void updateMR(gsl_rng* r, double delta, int maxNt, bool useShareBottleneck);
	void updateNt(gsl_rng *r, double delta, int maxNt, bool useSharedBottleneck);
	Model getModelRep(int repIndex);
	void updateBICbest(double BICnew);
	void print();
	void printNew();

	void setqBmean(std::vector<std::vector<std::vector<double> > >& qBmean);
	void setqBvar(std::vector<std::vector<gsl_matrix*> >& qBvar); //Note this one sets it without using the copy method!!

};

#endif /* modelMR_hpp */
