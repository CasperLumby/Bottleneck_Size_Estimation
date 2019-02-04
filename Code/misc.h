//Include guard
#ifndef MISC_H
#define MISC_H


//Forward declared dependencies

//Included dependencies
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include "sequence.h"
#include "diploidSequence.h"
#include "modelMR.hpp"
#include "model.hpp"
#include "simParamPH.h"
#include "anaParam.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"

void rescale(std::vector<double> &vec);
void rescaleZeroAccepted(std::vector<double> &vec);
void rescaleMin(std::vector<double> &vec, double min);
void printDoubleVector(std:: vector<double> &vec);
std::string printDoubleVectorToString(std:: vector<double> &vec);
void printIntVector(std:: vector<int> &vec);
std::string printIntVectorToString(std:: vector<int> &vec);
std::vector<int> multinomialSampling(int N, std::vector<double> p, const gsl_rng *r);
int sumOfVector(std::vector<int> intVec);
double sumOfVector(std::vector<double> doubleVec);
std::vector<double> subtractVectorsDouble(std::vector<double> &a, std::vector<double> &b);
std::vector<double> addVectorsDouble(std::vector<double> &a, std::vector<double> &b);
double my_gsl_linalg_LU_det (gsl_matrix * LU, int signum);
double determinant(gsl_matrix *A);
double logMultivariateNormalPDF(std::vector<double> &x, std::vector<double> &mu, gsl_matrix *sigma);
double logMultivariateNormalPDFPrint(std::vector<double> &x, std::vector<double> &mu, gsl_matrix *sigma);
void printMatrix(gsl_matrix *A);
void printMatrixMathematica(gsl_matrix *A);
std::string printMatrixMathematicaToString(gsl_matrix* A);
Sequence randomSequence(int length, gsl_rng *r);
std::vector<Sequence> generateRandomFullHaplotypes(int length, int number, gsl_rng *r);
std::vector<Sequence> generateMajorMinorFullHaplotypes(int length, int number, gsl_rng *r);
bool identicalSeqs(Sequence & a, Sequence & b);
bool identicalIntVectors(std::vector<int> & a, std::vector<int> & b);
bool nonZeroIntVector(std::vector<int> & a);
bool nonZeroDoubleVector(std::vector<double> & a);

void setupRNG(gsl_rng *r, unsigned long int seed);

DiploidSequence generateRandomDiploidSequence(int length, gsl_rng *r);
DiploidSequence generateRandomSemiConstrainedDiploidSequence(Sequence constraint, gsl_rng *r);
DiploidSequence generateRandomSemiConstrainedDiploidSequence(DiploidSequence constraint, gsl_rng *r);
std::vector<Sequence> generateFullHapsFromAlleles(DiploidSequence &fullHapAlleles);


void randomiseSequenceVector(std::vector<Sequence> & s, gsl_rng *r);

std::vector<int> pickRandomInts(gsl_rng *r, int max, int length);

void addRandom(std::vector<double>& oldVec, std::vector<double>& newVec, double delta, gsl_rng *r);
void addRandomMin(std::vector<double>& oldVec, std::vector<double>& newVec, double delta, gsl_rng *r, double min);
void addRandom(std::vector<double>& vec, double delta, gsl_rng *r);
double logMultinomialProb(std::vector<double> &freq, std::vector<int> & obs);
void findLogFactorial(std::vector<double> & fact_store, int N);
double logDirMultProbC(double C, std::vector<int> &obs, std::vector<double> &inf, std::vector<double> &fact_store);
double logDirMultProb(std::vector<double>& alpha, std::vector<int> &x, int& n);
double logBetaBinProb(double alpha, double beta, int x, int n);
double BetaBinCDF(double alpha, double beta, double threshold, int n);

bool fileExists(std::string &fileName);
std::vector<int> DirMultSampling(int N, std::vector<double> & freqs, double C, const gsl_rng *r);

//Mean, no var
double computeLSelection(std::vector<int> &obsB, std::vector<int> &obsA, int Nt, AnaParam* ap, std::vector<double> &qBfhFreqs, std::vector<double> &hapFit, std::vector<std::vector<int> >* contribs, bool print);
double computeLNeutral(std::vector<int>& obsB, std::vector<int>& obsA, int Nt, AnaParam* ap, std::vector<double>& qBfhFreqs, std::vector<std::vector<int> >* contribs, bool print);

//Mean and var
double computeLSelTVar(std::vector<int>& obsA, int Nt, AnaParam* ap, std::vector<double>& qBfhFreqs, gsl_matrix* var, std::vector<double>& hapFitT, std::vector<std::vector<int> >* contribs, bool print);
double computeLSelTVarOld(std::vector<int>& obsA, int Nt, AnaParam* ap, std::vector<double>& qBfhFreqs, gsl_matrix* var, std::vector<double>& hapFitT, std::vector<std::vector<int> >* contribs, bool print);
double computeLSelTSelGVar(std::vector<int>& obsA, int Nt, AnaParam* ap, std::vector<double>& qBfhFreqs, gsl_matrix* var, std::vector<double>& hapFitT, std::vector<double> &hapFitG, std::vector<std::vector<int> >* contribs, bool print);
double computeLSelTSelGVarOld(std::vector<int>& obsA, int Nt, AnaParam* ap, std::vector<double>& qBfhFreqs, gsl_matrix* var, std::vector<double>& hapFitT, std::vector<double> &hapFitG, std::vector<std::vector<int> >* contribs, bool print);
double computeLNeutralVar(std::vector<int>& obsA, int Nt, AnaParam* ap, std::vector<double>& qBfhFreqs, gsl_matrix* var, std::vector<std::vector<int> >* contribs, bool print);
double computeLNeutralVarAdvanced(std::vector<int>& obsA, int& deltaDays, int Nt, AnaParam* ap, std::vector<double>& qBfhFreqs, gsl_matrix* var, std::vector<std::vector<int> >* contribs, bool print);
double computeLNeutralVarReduced(std::vector<int>& obsA, int Nt, AnaParam* ap, std::vector<double>& qBfhFreqs, gsl_matrix* var, std::vector<std::vector<int> >* contribs, bool print);
double computeLNeutralSelGVar(std::vector<int>& obsA, int Nt, AnaParam* ap, std::vector<double>& qBfhFreqs, gsl_matrix* var, std::vector<double>& hapFitG, std::vector<std::vector<int> >* contribs, bool print);
double computeLNeutralSelGVarOld(std::vector<int>& obsA, int Nt, AnaParam* ap, std::vector<double>& qBfhFreqs, gsl_matrix* var, std::vector<double>& hapFitG, std::vector<std::vector<int> >* contribs, bool print);


void constructMatrixM(std::vector<double> x, gsl_matrix *M);
void constructMatrixMoriginal(std::vector<double> x, gsl_matrix *M);

std::vector<double> computeHapFit(std::vector<Sequence> &fullHaps, Sequence &selModel, std::vector<double> &selCoefsNew);
std::vector<double> computeSimulationHapFitT(std::vector<Sequence> &fullHaps, SimParamPH *spPH);
std::vector<double> computeSimulationHapFitG(std::vector<Sequence> &fullHaps, SimParamPH *spPH);

std::vector<DiploidSequence> createOneLocusCollapsed(std::vector<std::vector<DiploidSequence> > &oneLocusModels);

std::vector<Model> generateModelsFromPreviousBest(std::vector<DiploidSequence> &collapsedFullHaps, Model &bestModelBefore); 
std::vector<ModelMR> generateModelsFromPreviousBest(std::vector<DiploidSequence> &collapsedFullHaps, ModelMR &bestModelBefore); 

int binomialCoeff(int n, int k);
int logBinCoeffStirling(int n, int k);

int getNextLeft(std::vector<int> alreadyPicked, int max);
std::vector<Sequence> findRemainingSeqs(std::vector<Sequence> subset, std::vector<Sequence> set);

DiploidSequence constructDiploidFromFullHap(std::vector<Sequence> &fhs, std::vector<double> &freqs);
DiploidSequence constructDiploidFromFullHap(std::vector<Sequence> &fhs);
bool noEntryTheSameExceptZero(std::vector<double> &v);

std::vector<std::vector<int> > getAllCombs(int N, int k);

std::vector<std::vector<int> > getLinkedPairs(std::vector<Sequence> &haps);

std::vector<std::string> split(const std::string& s, char delim);

//Related to gaining information on memory usage
int parseLine(char* line);
int getValue();

bool stringToBool(std::string& s);
bool compareBySize(const std::vector<int>& a, const std::vector<int>& b);

//gamma and delta functions for compound solutions
double gamma_n(int n, double lambda, int Nt);
double delta_n(int n, double lambda, int Nt);
double qPochhammer(double a, double q, int n);

#endif
