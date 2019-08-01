#include <algorithm>
#include <iostream>
#include <list>
#include <string>
#include <array>
#include <fstream>
#include <math.h> 
#include <vector>
#include <sstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <numeric>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_linalg.h>
#include <time.h> 
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>

using namespace std;

//Forward declaritions
void printVector(vector<double> &vec);

//Normalise the frequency of haplotypes to 1
void NormaliseFreqs(vector<double> &init_freqs)
{
	double tot = 0;
	for (unsigned int i = 0; i<init_freqs.size(); i++) {
		tot = tot + init_freqs[i];
	}
	for (unsigned int i = 0; i<init_freqs.size(); i++) {
		init_freqs[i] = init_freqs[i] / tot;
	}
}
//Check if a file exists. Give full path to file
bool fileExists (string& fileName) 
{
	ifstream f(fileName.c_str());
	if (f.good()) 
	{
		f.close();
		return true;
	}
	else 
	{
		f.close();
		return false;
	}
}

//Takes a vector of frequencies and returns a vector of swaps
vector<vector<int> > findSwaps(std::vector<int>& extinctions, int dimensions) {

	vector<vector<int> > swaps;

	//Go through extinctions one at a time (largest first), then swap with furhest to the right non-extinction
	//Example: All sites = {0,1,2,3,4,5,6}, extinction sites = {5,3,1}.
	//Then first swap 5 with 6, so get order {0,1,2,3,4,6,5} and swap pair {5,6}
	//Then swap 3 with 5, i.e. get {0,1,2,6,4,3,5} and swap pairs {5,6}, {3,5}
	//Finally swap 1 with 4, i.e. get {0,4,2,6,1,3,5} and swap pairs {5,6}, {3,5}, {1,4}.
	//We could further swap {0,4,2,6} to get them in order, but it doesn't matter much/is time consuming to do additional swaps
	for(unsigned int i=0; i<extinctions.size(); i++) { //Loop over extinction entries

		//Check if extinction is too far to the left, if so, move it to right position.
		//E.g. if extinction[i] = 5 and dimensions = 7 and (i+1) = 1, then 5 < 7-1 is true, so make swap with site 7-1=6
		//Similarly, if extinction[i] = 3 and dimensions = 7 and (i+1) = 2, then 3 < 7-2 is true, so make swap with position 7-2=5
		if(extinctions[i] < (int)(dimensions-(i+1))) {
			
			int swapPartner = dimensions - (i+1);
			vector<int> swap {extinctions[i], swapPartner};
			swaps.push_back(swap);
		}
	}

	return swaps;

}

//Perform a set of swaps (e.g. swaps = {{5,6}, {3,5}, {1,4}}) for a matrix mat using gsl methods
void swapMatrix(gsl_matrix *mat, vector<vector<int> > swaps) {


	for(unsigned int i=0; i<swaps.size(); i++) {

		int swap1 = swaps[i][0]; //First entry for swapping, e.g. 5
		int swap2 = swaps[i][1]; //Second enttry for swappging, e.g. 6

		gsl_matrix_swap_rows(mat, swap1, swap2);
		gsl_matrix_swap_columns(mat, swap1, swap2);
		
	}
}

//Perform a set of swaps (e.g. swaps = {{5,6}, {3,5}, {1,4}}) for a double vector (e.g. {a,b,c,d,e,f,g} becomes {a,e,c,g,b,d,f}
void swapVector(vector<double>& v, vector<vector<int> >& swaps) {

//	cout << "Swapping vector: "; printVector(v);
//	cout << "Swap vector is: ";
//	for(unsigned int i=0; i<swaps.size(); i++) {
//
//		for(unsigned int j=0; j<swaps[i].size(); j++) {
//			cout << swaps[i][j] << " ";
//		}
//		cout << "\n";
//	}


	for(unsigned int i=0; i<swaps.size(); i++) {

		int swap1 = swaps[i][0]; //First entry for swapping, e.g. 5
		int swap2 = swaps[i][1]; //Second enttry for swappging, e.g. 6

		iter_swap(v.begin() + swap1, v.begin() + swap2); //Defined in algorithm
//		cout << "Vector after swap " << i << ": "; printVector(v);
	}
}

//Construct the matrix M = diag(q) - q*q^transpose
void constructMatrixM(vector<double> x, gsl_matrix *M) {

    //Create subset of vec
    int dim = (int) x.size();

	//Hardcode the very simple cases as it gives a considerable speed-up. See mMatrixChecker.cpp.
	if(dim==1) { //Very simple case, if x={x1} then M = {{ x1-x1^2 }}, hardcode for efficiency. About 1.5 times as fast

		double value = x[0] - x[0]*x[0];
		gsl_matrix *mat = gsl_matrix_alloc(dim,dim);
		gsl_matrix_set(mat,0,0,value);
		gsl_matrix_memcpy(M,mat); //Copy mat into M, entry by entry
		gsl_matrix_free(mat); //Clean up
		return;

	} else if(dim==2) { //Also simple case. Here M[{x1,x2}] = {{x1 - x1^2, -x1 x2}, {-x1 x2, x2 - x2^2}}. About 1.75 times as fast.

		gsl_matrix *mat = gsl_matrix_alloc(dim,dim);
		gsl_matrix_set(mat,0,0,x[0]-x[0]*x[0]);
		gsl_matrix_set(mat,1,0,-x[0]*x[1]);
		gsl_matrix_set(mat,0,1,-x[0]*x[1]);
		gsl_matrix_set(mat,1,1,x[1]-x[1]*x[1]);
		gsl_matrix_memcpy(M,mat); //Copy mat into M, entry by entry
		gsl_matrix_free(mat); //Clean up
		return;

	}
    
    //Create diagonal matrix with elements of x
    gsl_matrix *diagX = gsl_matrix_alloc(dim,dim);
    for(int j=0; j<dim; j++) {
        for(int k=0; k<dim; k++) {
            
            if(j==k) {
                gsl_matrix_set(diagX,j,k,x[j]);
            } else {
                gsl_matrix_set(diagX,j,k,0);
            }
        }
    }
    
   /*
     / Create outer product of x*x^transpose.
     / This method does not exist in CBLAS library, but can be obtained otherwise.
     / In particular, Outer(a,b) = Matrix(columns=a) * DiagonalMatrix(b).
     */
    //Create matrices
    gsl_matrix *A = gsl_matrix_alloc(dim,dim);
    for(int j=0; j<dim; j++) {
        for(int k=0; k<dim; k++) {
            gsl_matrix_set(A,j,k,x[j]);
        }
    }
    
    //Multiply them and store in XX
    gsl_matrix *XX=gsl_matrix_alloc(dim,dim);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                    1.0, A, diagX,
                    0.0, XX);
    
    //Create M=diag(x) -xxp^t
    gsl_matrix_sub(diagX,XX); //diagX=diagX-XX, i.e. store result in diagPb
    
    //Clean up
    gsl_matrix_memcpy(M,diagX); //Copy to M, entry by entry
    gsl_matrix_free(diagX);
    gsl_matrix_free(A);
    gsl_matrix_free(XX);   
}

//Substract two double vectors
vector<double> subtractVectors(vector<double> &a, vector<double> &b) {
    
	vector<double> result;
	//Assume equal length
	for(unsigned int i=0; i<a.size(); i++) {
		result.push_back(a[i]-b[i]);
	}
	return result;
}

//Compute log likelihood for multivariate normal distribution in x with mean mu and variance sigma
double logMultivariateNormalPDF(vector<double> &x, vector<double> &mu, gsl_matrix *sigma) {
    
    //Set up permutation and matrix copy
    int dim= (int) sigma->size1;
    vector<double> xMinusMu = subtractVectors(x,mu);
    gsl_permutation *p = gsl_permutation_calloc(dim);
    int signum;
    gsl_matrix * tmp_ptr = gsl_matrix_calloc(dim,dim);
    gsl_matrix_memcpy(tmp_ptr,sigma); //Copy elements of sigma into tmp_ptr
    
    //Get LU decomposition
    gsl_linalg_LU_decomp(tmp_ptr,p,&signum);
    
    //Get determinant
    double det = gsl_linalg_LU_det(tmp_ptr, signum);
//	cout <<  "Determinant: " << det << "\n";
//	if(det==0) { det=1e-10; }
    
    //Get inverse
    gsl_matrix *sigmaInv = gsl_matrix_alloc(dim,dim);
    gsl_set_error_handler_off();
    int status = gsl_linalg_LU_invert(tmp_ptr,p,sigmaInv);
    if (status) {
        
	//Clean up, then return negative infinity
	gsl_matrix_free(sigmaInv);
	gsl_matrix_free(tmp_ptr);
	gsl_permutation_free(p);
    

        cout << "Matrix not positive definite. Returning probability of neg infinity.\n";
//	cout << "Determinant should be 0 in this case. Determinant was: " << det << "\n";
        return  -numeric_limits<double>::max();
    }
    
    
    //double logPrefactor= -0.5*log(pow(2*M_PI,dim)*det);
    double logPrefactor= -0.5*dim*log(2*M_PI) -0.5*log(det);
    
    //Convert xMinusMu to a gsl_vector
    gsl_vector *xMinusMuGSL = gsl_vector_alloc(dim);
    for(unsigned int i=0;i<xMinusMu.size();i++) {
        gsl_vector_set(xMinusMuGSL,i,xMinusMu[i]);
    }
    
    //Perform matrix*vector multiplication
    gsl_vector *sigmaInvXminusMu = gsl_vector_alloc(dim);
    gsl_blas_dgemv(CblasNoTrans,1.0,sigmaInv,xMinusMuGSL,0,sigmaInvXminusMu);
    
    //Perform vector*vector multiplication
    double dotProd;
    gsl_blas_ddot(xMinusMuGSL,sigmaInvXminusMu,&dotProd);
    
	//Clean up
	gsl_vector_free(xMinusMuGSL);
	gsl_vector_free(sigmaInvXminusMu);
	gsl_matrix_free(sigmaInv);
	gsl_matrix_free(tmp_ptr);
	gsl_permutation_free(p);
    
	return logPrefactor -0.5*dotProd;
}



//Use the fully continuous method for computing a likelihood
//This method alters the mean and the observations x!
double computeLikelihoodCont(vector<double> x, vector<double>& mean, gsl_matrix* var) {

	//Computing likelihoods depend on dimensionality
	double logLikelihood = -numeric_limits<double>::max(); //Set to neg infinity

	int dim = (int) mean.size();
	if(dim > 2) { //Dimension larger than 2


		/*
		* Reduce dimensionality by 1 to ensure non-degeneracy
		*/
		mean.pop_back(); //WLOG remove last haplotype
		x.pop_back(); //WLOG remove last haplotype
		gsl_matrix_view varReducedView = gsl_matrix_submatrix(var, 0,0, dim-1, dim-1); //Remove last row and column
		gsl_matrix* varReduced = gsl_matrix_alloc(dim-1, dim-1);
                gsl_matrix_memcpy(varReduced, &(varReducedView.matrix)); //Copy the reduced matrix view into a "normal" matrix varReduced
		

		//Compute likelihood
		logLikelihood = logMultivariateNormalPDF(x,mean,varReduced);


	} else { //I.e. dim==2 or dim==1 (same computation)


		//In the case where there is a single haplotype in a set (i.e. dim=1),
		//this can be thought of as a case with a single observed  haplotype
		//and one (or several) unobserved haplotypes.
		//If all loci in the haplotype set are variants then there MUST
		//exist at least one other haplotype containing the unobserved alleles,
		//however, this (or these) haplotype(s) are not observed.
		//As such, WLOG we can consider this a case of dim=2, which gets reduced to a one dimensional system
		//under reduction, i.e. similar to the dim==2 case.
		double meanSingle = mean[0];
		double varSingle = gsl_matrix_get(var,0,0);
		double stDevSingle = sqrt(varSingle);	

		logLikelihood = -(x[0]-meanSingle)*(x[0]-meanSingle)/(2*varSingle) - log(stDevSingle*sqrt(2*M_PI)); //log L
	
	}

	return logLikelihood;
}

void printVector(vector<double> &vec) {
        cout.precision(25);
    if(vec.size() > 0) {
        cout << "Vector = { ";
        for(unsigned int i=0;i<vec.size()-1;i++) {
            cout << vec[i] << " , ";
        }
        cout << vec[vec.size()-1] << " }\n";
    } else {
        cout << "ERROR: Vector has zero length.\n";
    }
        cout.precision(6); //Change back
}


//Compute likelihood
double computeLikelihood(vector<double>& qStarB, vector<double>& qStarA, gsl_matrix* varB, gsl_matrix* varA, int Nt) {

	/*
	 * First compute compound solution. It is of form
	 * 
	 * 	var[qA*] = varA + gamma*M(qB*) + delta*varB
	 * 	
	 * where gamma = (Nt+Ng-1)/(Nt*NG) and delta = (Nt*Ng-Nt-Ng+1)/(Nt*NG)
	 *
	 */
	vector<double> x = qStarA;
	vector<double> mean = qStarB;
	gsl_matrix* var = gsl_matrix_alloc(mean.size(), mean.size());

	//Create factors gamma and delta
	double Ng = Nt*22; //22 fold growth
	double NtNg = Nt*Ng;
	double gamma = (Nt+Ng-1)/NtNg;
	double delta = (NtNg-Nt-Ng+1)/NtNg;

	//Create matrix gamma*M(qB*)
	gsl_matrix* gammaMqStarB = gsl_matrix_alloc(mean.size(),mean.size());
        constructMatrixM(mean, gammaMqStarB);
	gsl_matrix_scale(gammaMqStarB, gamma);


	//Create delta*varB
	gsl_matrix* deltaVarB = gsl_matrix_alloc(mean.size(), mean.size());
	gsl_matrix_memcpy(deltaVarB, varB);
	gsl_matrix_scale(deltaVarB, delta);

	//Copy varA into var
	gsl_matrix_memcpy(var, varA);
	
	//Add gamma*M(qB*) to var
	gsl_matrix_add(var, gammaMqStarB);

	//Add delta*varB to var
	gsl_matrix_add(var, deltaVarB);


	//Free up unused variables
	gsl_matrix_free(gammaMqStarB);
	gsl_matrix_free(deltaVarB);

	//Check if there are any extinctions
	double extinctionThreshold = 1e-6; //Below this qualfies as extinction
	int n = 0; //Number of extinctions
	vector<int> extinctions;

	for(int i=(int)x.size()-1; i>=0; i--) { //Get extinction entries in descending order, e.g. {5,3,1}

                if(x[i] <= extinctionThreshold) {

                        extinctions.push_back(i);
			n++;
                }
        }
	cout << "There were " << n << " extinctions.\n";

	if(n == 0) {

		//Compute likelihood the fully continuous (multivariate normal) way
		return computeLikelihoodCont(x, mean, var);

	} else { //Compute likelihood the continuous-discrete way

		//Start by finding swaps
		int k = x.size(); //Dimensions k
		vector<vector<int> > swaps = findSwaps(extinctions, k);

		
		//Then update x, mean and variances
		swapVector(x, swaps);
		swapVector(mean, swaps);
		swapMatrix(var, swaps);

		//Note that we have altered the original variables. This may not be ideal, consider at a later stage

		/*
		 * We want to compute the likelihood L = P(surviving haplotypes | extinct haplotypes) * P(extinct haplotypes).
		 * We note that the above is written in probability space, in log likelihood space the terms are summed, not multiplied.
		 */

		
		//First we create the conditional distribution for the surviving haplotypes conditional on the extinct haplotypes, see https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Conditional_distributions.
		//Here we have a system of the form:
		// Dimensionality of k = x.size()
		// Number of extinctions n = numExtinctions
		// Inferred frequencies x = q*A = [q1, q2] with q1 being of dim (k-n)x1 and q2 being of dim nx1
		// Mean = [m1, m2] with m1 being of dim (k-n)x1 and m2 being of dim nx1
		// var = {{v11, v12}, {v21, v22}} with v11 being of dim (k-n)x(k-n), v12 with dim (k-n)xn, v21 with dim nx(k-n) and v22 with dim nxn
		// Extinction outcome a = x[k-n, k] (i.e. basically {0,0,0,0...,0} for n entries. Note, this is not entirely true, a some entries can be e.g. 10^-8 and still count as extinction)
		// 
		// Then the conditional distribution has mean
		// meanCond = m1 + v12 . (v22)^-1 . (a - m2)   (here we use '.' to indicate matrix multiplication)
		//
		// And it has variance
		// varCond = v11 - v12 . (v22)^-1 . v21
		
		//First create m1, m2, and a:
		vector<double> m1(mean.begin(), mean.begin() + (k-n));
		vector<double> m2(mean.begin() + (k-n), mean.end());
		vector<double> a(x.begin() + (k-n), x.end());
		
		//Then create v11: (which we will call 'varCond' - the reason for this will become clear below)
		gsl_matrix_view v11View = gsl_matrix_submatrix(var, 0,0, (k-n), (k-n)); //Get sub matrix with upper-left element at (0,0), numRows=numColumns = (k-n)
		gsl_matrix * varCond = gsl_matrix_alloc((k-n),(k-n));
		gsl_matrix_memcpy(varCond, &(v11View.matrix)); //Copy the v11 view into a "normal" matrix v11
	
		//Then create v12:
		gsl_matrix_view v12View = gsl_matrix_submatrix(var, 0,(k-n), (k-n), n); //Get submatrix with upper-left element at (0,(k-n)), numRows=(k-n), numCols=n
		gsl_matrix * v12 = gsl_matrix_alloc((k-n),n);
		gsl_matrix_memcpy(v12, &(v12View.matrix)); //Copy the v12 view into a "normal" matrix v12
	
		//Then create v21:
		gsl_matrix_view v21View = gsl_matrix_submatrix(var, (k-n), 0, n, (k-n)); //Get submatrix with upper-left element at ((k-n),0), numRows=n, numCols=(k-n)
		gsl_matrix * v21 = gsl_matrix_alloc(n,(k-n));
		gsl_matrix_memcpy(v21, &(v21View.matrix)); //Copy the v21 view into a "normal" matrix v21
	
		//Then create v22:
		gsl_matrix_view v22View = gsl_matrix_submatrix(var, (k-n), (k-n), n, n); //Get submatrix with upper-left element at ((k-n),(k-n)), numRows=numColumns = n
		gsl_matrix * v22 = gsl_matrix_alloc(n,n);
		gsl_matrix_memcpy(v22, &(v22View.matrix)); //Copy the v22 view into a "normal" matrix v22
	

		/////// Now we create the conditional mean
		//First create a gsl vector representing (a-m2)
		gsl_vector *aMinusM2 = gsl_vector_alloc(n);
		for(int i=0;i<n;i++) {
			gsl_vector_set(aMinusM2,i,a[i] - m2[i]);
		}


		//Then we want to create the inverse of v22:
		//First set up permutation and signum
		gsl_permutation *p = gsl_permutation_alloc(n);
		int signum;
    
		//Then compute the LU decomposition
		gsl_linalg_LU_decomp(v22,p,&signum); //Note that v22 is changed in this process!
    
		//Finally create inverse of v22
		gsl_matrix *v22Inv = gsl_matrix_alloc(n,n);
		gsl_set_error_handler_off(); //Stops the default error handler from exiting the program - instead we deal with the error ourselves
		int status = gsl_linalg_LU_invert(v22,p,v22Inv); //Status of 0 equal success
		if (status) {
        
			//Error occurred, so clean up and exit
			cout << "V22 matrix not positive definite, so cannot invert. Exiting.\n";
			exit(1);
		}


		//Next we do the multiplications, starting with v12 . (v22)^1
		gsl_matrix* v12v22Inv = gsl_matrix_alloc(k-n,n);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, v12, v22Inv, 0.0, v12v22Inv);
    

		//Followed by the matrix-vector multiplication (v12.(v22)^1) . (a-m2)
		//Actually, we can be a bit clever. The matrix-vector multiplication in gsl is as follows: y=alpha*A*x +beta*y.
		//In other words, if we pre-define a vector y, we can do both matrix-vector multiplication of A times x, but also add y.
		//This is exactly what we need to compute meanCond = m1 + v12 . (v22)^-1 . (a - m2), i.e. we can let y=m1, A=(v12.(v22)^-1) and x=(a-m2).
		gsl_vector *y = gsl_vector_alloc(k-n);
		for(int i=0; i<(k-n); i++) {
			gsl_vector_set(y,i,m1[i]);
		}
		gsl_blas_dgemv(CblasNoTrans, 1.0, v12v22Inv, aMinusM2, 1.0, y); //Result is now stored in y
		
		//Get conditional mean as 'normal' C++ vector
		vector<double> meanCond;
		for(int i=0; i<(k-n); i++) {
			meanCond.push_back(gsl_vector_get(y,i));
		}


		//////// Now we create the conditional variance:  varCond = v11 - v12 . (v22)^-1 . v21
		//Here we can be a bit clever. The matrix-matrix multiplication in gsl is as follows: C=alpha*A*B+beta*C.
		//In other words, if we pre-define a matrix C we can do all the computation in one go. To this end, we
		//associate C=v11, alpha=-1, A=v12v22Inv, B=v21, beta=1.0. To avoid renaming matrices (through wasteful copying) we
		//will refer to C=v11 as 'varCond' (we already did this at the defintion of v11 above).
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,-1.0,v12v22Inv,v21, 1.0, varCond);

		//Create conditional observations, x (i.e. q*A subsetted to match the mean)
		vector<double> xCond(x.begin(), x.begin()+(k-n));

		//Compute the continuous likelihood
		double logLcont = computeLikelihoodCont(xCond, meanCond, varCond);

		//Compute discrete likelihood
		//Noticing that x_i=0 for all the extinct haplotypes, the binomial/multinomial likelihood simplifies to
		//P(q_1,q_2,...,q_n) = (1-q_1-q_2-...-q_n)^(Nt) for n extinct haplotypes. As such
		//L(q_1,q_2,...,q_n) = Nt log(1-q_1-q_2-..._q_n) (log likelihood)
		double sumFreqs = 0;
		for(int i = (k-n); i<k; i++) {
			sumFreqs += mean[i];
		}
		double logLdisc = Nt*log(1-sumFreqs);
		cout << "Nt: " << Nt << "\tlogLDisc: " << logLdisc << "\tlogLcont: " << logLcont << "\n";
		
		//Free matrices + views + permutations etc.
		gsl_matrix_free(var);
		gsl_permutation_free(p);
		gsl_matrix_free(varCond);
		gsl_matrix_free(v12);
		gsl_matrix_free(v21);
		gsl_matrix_free(v22);
		gsl_vector_free(aMinusM2);
		gsl_matrix_free(v22Inv);
		gsl_matrix_free(v12v22Inv);
		gsl_vector_free(y);

		//cout << "Nt: " << Nt << "\tlogLcont: " << logLcont << " \tlogLDisc: " << logLdisc << "\n";
		return logLcont + logLdisc;
	}

}



int main(int argc, char* argv[]) {

	//first check if the input format is right
	if(argc != 3) { 
		cout << "Error in the input format: " << endl << "Please provide outcome_1.txt file (i.e. your reconstructed haplotypes) and test_$t.txt files (i.e. concatenated files for inferred frequencies q** for all the 8 segments)";
		return false;
	}
	
	string qStarFilePath = argv[1];
	string qStarStarFilePath =  argv[2];

	cout << "You provided the following input file paths:\n";
	cout << "q* file path: " << qStarFilePath << "\n";
	cout << "q** file path: " << qStarStarFilePath << "\n";
	cout << "Bottleneck inference will now be initiated.\n";
	
	freopen("maximum_likelihood.txt","w",stdout); //This writes the stdout output to maximum_likelihood.txt

	//Progress with inference only if q* file exists, otherwise print error
	if (fileExists(qStarFilePath) == true) {


		/*
		 * This file looks like the following:
		 *
		 * 1	TTA	0.248196	0
		 * 2	ATG	0.00148182	0.0393609
		 * 3	ATA	0.628543	0.959854
		 * 4	TCA	0.0850215	0
		 * 5	ACA	0.0366405	0.000516205
		 *
		 * i.e. #haplotype_number haplotype frequency_1 frequency_2
		 */

		//Save data as a vector of haplotype entries
		vector<string> qStarData;
		ifstream qStarFile;
		qStarFile.open(qStarFilePath); 
		while (!qStarFile.eof()) { //Read one line at a time
			string line;
			getline(qStarFile, line);
			qStarData.push_back(line);
		}
		if (qStarData.empty())
		{
			cout << "File for q* was empty. Exiting.\n";
			return false;
		}
		qStarFile.close();

		//Now convert each entry in the q* data into a vector of strings, i.e. split by whitespace
		vector < vector<string> > qStarDataSplit;
		for (unsigned int i = 0; i < qStarData.size(); i++) { //Loop over each entry

			vector<string> row;
			istringstream iss(qStarData[i]);
			string term;
			while (iss >> term)
			{
				row.push_back(term);
			}
			if (!row.empty())
			{
				qStarDataSplit.push_back(row);
			}
		}
		
		//Now we get the q* frequencies before and after transmission
		vector<double> qStarBefore;
		vector<double> qStarAfter;
		vector<string> haplotype_IDs; // we record the haplotype IDs in outcome_1.txt in order to later match them with their corresponding q** frequency 

		for (unsigned int i = 0; i < qStarDataSplit.size(); i++) { //Loop over entries in q* file

			//Assume file of format: #haplotype_number haplotype frequency_1 frequency_2 )
			
			//Get before and after frequencies as strings
			string qStarBeforeString = qStarDataSplit[i][2];
			string qStarAfterString = qStarDataSplit[i][3];
		
			//Then convert to doubles and add to qStar vectors
			double qStarBeforeDouble = atof(qStarBeforeString.c_str());
			double qStarAfterDouble = atof(qStarAfterString.c_str());
			qStarBefore.push_back(qStarBeforeDouble);
			qStarAfter.push_back(qStarAfterDouble);

			//Also note down the haplotype ID associated with these frequencies
			haplotype_IDs.push_back(qStarDataSplit[i][0]);
		}
		int numHaps = (int) haplotype_IDs.size();
		

		/*
		 * Now we read in the q** inferences. These are concatenated together, e.g. there might be n=100 sets of inferences in this file.
		 * For example for n=2 replicates:
		 *
		 * 0       0.254717        0.0387761
		 * 1       0       0
		 * 2       0.63447 0.95342
		 * 3       0.077882        0
		 * 4       0.0329313       0.00780431
		 * 0       0.244672        0
		 * 1       0       0
		 * 2       0.642081        0.999423
		 * 3       0.0882288       0
		 * 4       0.0250186       0.000577298
		 *
		 * We here read in the data as a vector of strings (entries, e.g. "0       0.254717        0.0387761" )
		 *
		 */
		vector<string> qStarStarData;
		ifstream qStarStarFile;
		qStarStarFile.open(qStarStarFilePath); //get the full list of concatenated test_i.txt files that are currently stored in Transmission1_bottleneck/segment_i/
		while (!qStarStarFile.eof()) { //Read one line at a time

			string line;
			getline(qStarStarFile, line);
			qStarStarData.push_back(line);
		}
		if (qStarStarData.empty())
		{
			cout << "File for q** was empty. Exiting.\n";
			return false;
		}
		qStarStarFile.close();

		//Now convert each entry in the q** data into a vector of strings, i.e. split by whitespace
		vector<vector<string> > qStarStarDataSplit;
		for (unsigned int i = 0; i < qStarStarData.size(); i++) { //Loop over each entry

			vector<string> row;
			istringstream iss(qStarStarData[i]);
			string term;
			while (iss >> term)
			{
				row.push_back(term);
			}
			if (!row.empty())
			{
				qStarStarDataSplit.push_back(row);
			}
		}
		
		//Get dimensions of q** (i.e. how many replicates we did the q** inference over)
		int numReps;
		if((qStarStarDataSplit.size())%(qStarDataSplit.size()) == 0) { //Check if size of q** is a multiple of size of q*

			numReps = (qStarStarDataSplit.size())/(qStarDataSplit.size()); //this gives the total number of q**s generated for each haplotype in outcome_1.txt

		} else { //If not, error

			cout << "Error: The number of entries in q** is not a multiple of the number of haplotypes in q*. Exiting.\n";
			exit(1);
		}
		
		/*
		 * We now calculate the variance in q*. We do this by assuming a diagonal covariance matrix with
		 * entries defined as var_{i,i} = sum_j (q*_j - q**_j)^2/n    where the sum j is from 1 to n.
		 *
		 */ 
		//First allocate matrices varB and varA and set to zero
		gsl_matrix* varB = gsl_matrix_alloc(numHaps,numHaps);
		gsl_matrix* varA = gsl_matrix_alloc(numHaps,numHaps);
		gsl_matrix_set_zero(varB);
		gsl_matrix_set_zero(varA);


		//Then calculate diagonal variances
		vector<double> varBDiagonal(numHaps, 0.0);
		vector<double> varADiagonal(numHaps, 0.0);
		for(int i=0; i<numReps; i++) { //Loop over all the replicates in q** 


			for(int j=0; j<numHaps; j++) { //Loop over the haplotypes in q*

			//	cout << "Entry (" << i << "," << j << "): " << qStarStarData[i*numHaps+j] << "\n";

				//q** data has format: hapID before_frequency after_frequency
				double qStarStarBefore_j = atof(qStarStarDataSplit[i*numHaps+j][1].c_str());
				double qStarStarAfter_j = atof(qStarStarDataSplit[i*numHaps+j][2].c_str());
				
				varBDiagonal[j] += pow((qStarBefore[j] - qStarStarBefore_j),2)/numReps;
				varADiagonal[j] += pow((qStarAfter[j] - qStarStarAfter_j),2)/numReps;
			}
		}


		//Update matrices
		for(int i=0; i<numHaps; i++) {

			gsl_matrix_set(varB, i, i, varBDiagonal[i]);
			gsl_matrix_set(varA, i, i, varADiagonal[i]);
		}


		//Next we compute likelihoods for bottlenecks in range [1,1000]
		vector<double> likelihoods;
		for(int Nt=1; Nt<=1000; Nt++) {

			double likelihood = computeLikelihood(qStarBefore, qStarAfter, varB, varA, Nt);
		
			likelihoods.push_back(likelihood);
		}

		//Next we identify the maximum likelihood bottleneck estimate
		double max = std::numeric_limits<double>::lowest(); //Start at lowest double value
		int bestNt = -1;
		for(unsigned int i=0; i<likelihoods.size(); i++) { //Loop over the likelihoods, index 0 corresponds to bottleneck of 1

			if(likelihoods[i] > max) {

				max = likelihoods[i];
				bestNt = i+1;
			}
		}

		//If maximum likelihood estimate found, print to screen and write likelihoods to file
		if(bestNt != -1) {

			cout << "Maximum likelihood bottleneck size was " << bestNt << " with likelihood of " << max << "\n";

			//Write likelihoods to file
			ofstream outFile;
			outFile.open("likelihood_distribution.txt");
			for(unsigned int i=0; i< likelihoods.size(); i++) {
		
				outFile << likelihoods[i] << "\n";
			}
			outFile.close();
		
		} else {

			cout << "Error, no maximum bottleneck was found. Exiting.\n"; exit(1);
		}		
	
		
		return 0;

	} else {
		//q* file doesn't exist, so can't infer bottlenecks
		cout << "The provided q* file either doesn't exist. Bottleneck inference interrupted. Exiting.\n";
		return false;
	}

} 
