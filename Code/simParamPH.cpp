
//Forward declared dependencies

//Included dependencies
#include "simParamPH.h"
#include <stdio.h>
#include <iostream>

using namespace std;

//Define constructors
SimParamPH::SimParamPH() : SimParam(), dim(3), C(200), format(""), filterData("Transmission"), repScenario(false), geneIndex(-1), numGenerations(1), deltaDays(1), growthFactor(22) {
    
}

SimParamPH::~SimParamPH() {} //Deconstructor

//Setters
void SimParamPH::setDim(int d) { dim = d; }
void SimParamPH::setC(double c) { C = c; }
void SimParamPH::setPathToFolder(string path) { pathToFolder = path; }
void SimParamPH::setRepID(string ID) { repID = ID; }
void SimParamPH::setRepScenario(bool b) { repScenario = b; }
void SimParamPH::setGeneIndex(int g) { geneIndex = g; }
void SimParamPH::setFormat(string f) { format = f; }
void SimParamPH::setFilterData(string fd) { filterData = fd; }
void SimParamPH::setGrowthFactor(int gf) { growthFactor = gf; }
void SimParamPH::setNumGenerations(int ng) { numGenerations = ng; }
void SimParamPH::setDeltaDays(int dd) { deltaDays = dd; }

//Getters
int SimParamPH::getDim() { return dim; }
double SimParamPH::getC() { return C; }
string SimParamPH::getPathToFolder() { return pathToFolder; }
string SimParamPH::getRepID() { return repID; }
bool SimParamPH::getRepScenario() { return repScenario; }
int SimParamPH::getGeneIndex() { return geneIndex; }
string SimParamPH::getFormat() { return format; }
string SimParamPH::getFilterData() { return filterData; }
int SimParamPH::getGrowthFactor() { return growthFactor; }
int SimParamPH::getNumGenerations() { return numGenerations; }
int SimParamPH::getDeltaDays() { return deltaDays; }

