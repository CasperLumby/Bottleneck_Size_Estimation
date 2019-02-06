
//Forward declared dependencies

//Included dependencies
#include "misc.h"
#include <iostream>
#include <gsl/gsl_linalg.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include "fstream"
#include <algorithm>
#include <limits>
#include "cubature.h"

#if defined(PCUBATURE)
#  define cubature pcubature
#else
#  define cubature hcubature
#endif


using namespace std;

//Defines quantities for numerical integration
const double tol = 1e-06; //Tolerance for numerical integration
const size_t maxEval = 1e07; //MaxEval == 0 signifies unlimited # of evaluations
const unsigned int intOutputDim = 1; //Number of output dimensions must be 1 - it has to be a scaler

//Rescales a vector of doubles to have a sum of 1 and no entry less than 10^-10
void rescale(vector<double> &vec) {
    
	//Get scaling factor
	double sum = 0;
	for (unsigned int i=0; i<vec.size();i++) {
		sum += vec[i];
	}
	if(sum > 0) { //If sum is zero, no need to use scalefactor
		double scaleFactor = 1.0/sum;
    
		//rescale
		for (unsigned int i=0;i<vec.size();i++) {
			vec[i] *= scaleFactor;
		}
	}

	//Check that no values are below 10^-11 as a result of rescaling
	bool changeOccurred = true;
	while(changeOccurred == true) {
	
		changeOccurred = false;

		for (unsigned int i=0; i<vec.size(); i++) {
		
			if(vec[i] < 1e-11) {

				//cout << "Vector < 1e-11: ";
				//printDoubleVector(vec);
				//cout << "\n";

				vec[i] = 1e-10;

				//Define vec2 to be subset of vec that excludes the ith entry
				vector<double> vec2;
				for(unsigned int j=0; j<vec.size(); j++) {
					if(j!=i) {
						vec2.push_back(vec[j]);
					}
				}
				
				if(vec2.size() == 1) { //Can't rescale if of length 1, just keep as is, or bump up to min size
					if(vec2[0] < 1e-11) { vec2[0] = 1e-10; }
				} else {
					rescale(vec2);
				}
			
				//Add new values to vec, but keep ith entry intact	
				for(unsigned int j=0; j<vec2.size(); j++) {
					if(j<i) {
						vec[j]=vec2[j];
					} else {
						vec[j+1]=vec2[j];
					}
				}

				rescale(vec);
				changeOccurred = true;
				break; //Leave for loop and start again
			}
		}
	}	
}

void rescaleZeroAccepted(vector<double> &vec) {

	 //Get scaling factor
	double sum = 0;
	for (unsigned int i=0; i<vec.size();i++) {
		sum += vec[i];
	}
	double scaleFactor = 1.0/sum;
    
	//rescale
	for (unsigned int i=0;i<vec.size();i++) {
		vec[i] *= scaleFactor;
	}
}

//Rescales a vector of doubles to have a sum of 1 and no entry less than min
//Min MUST be less than 0.3 (otherwise 1.5 factor is too small)
void rescaleMin(vector<double> &vec, double min) {
    
	//Get scaling factor
	double sum = 0;
	for (unsigned int i=0; i<vec.size();i++) {
		sum += vec[i];
	}

	if(sum == 0) { //If sum is zero, method breaks down. Return (1/n,1/n,...,1/n) with n the dimension of vec

		for(unsigned int i=0; i<vec.size(); i++) {

			vec[i] = 1.0/(double)vec.size();
		}

	} else { //The majority of cases go here
		double scaleFactor = 1.0/sum;
    
		//rescale
		for (unsigned int i=0;i<vec.size();i++) {
			vec[i] *= scaleFactor;
		}
	}


	//Check that no values are below 10^-11 as a result of rescaling
	bool changeOccurred = true;
	while(changeOccurred == true) {
	
		changeOccurred = false;

		for (unsigned int i=0; i<vec.size(); i++) {
		
			if(vec[i] < min) {

				//cout << "Vector < 1e-11: ";
				//printDoubleVector(vec);
				//cout << "\n";

				vec[i] = min*1.5; //Needs to be marginally bigger than min, otherwise iterations break down

				//Define vec2 to be subset of vec that excludes the ith entry
				vector<double> vec2;
				for(unsigned int j=0; j<vec.size(); j++) {
					if(j!=i) {
						vec2.push_back(vec[j]);
					}
				}

				if(vec2.size() == 1) { //Can't rescale if of length 1, just keep as is, or bump up to min size
					if(vec2[0] < min) { vec2[0] = min*1.5; }
				} else {
					rescaleMin(vec2, min);
				}			

				//Add new values to vec, but keep ith entry intact	
				for(unsigned int j=0; j<vec2.size(); j++) {
					if(j<i) {
						vec[j]=vec2[j];
					} else {
						vec[j+1]=vec2[j];
					}
				}

				rescaleMin(vec, min);
				changeOccurred = true;
				break; //Leave for loop and start again
			}
		}
	}	
}

//Method for printing vectors of doubles
void printDoubleVector(vector<double> &vec) {
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

string printDoubleVectorToString(vector<double> &vec) {
    
    string output;
    if(vec.size() > 0) {
        output += "{ ";
        for(unsigned int i=0;i<vec.size()-1;i++) {
            output += to_string(vec[i]);
		output += " , ";
        }
        output += to_string( vec[vec.size()-1]);
	output += " }";
    } else {
        cout << "ERROR: Vector has zero length.\n";
    }
    return output;
}


//Method for printing vectors of ints
void printIntVector(vector<int> &vec) {
    if(vec.size() > 0) {
        cout << "Vector = { ";
        for(unsigned int i=0;i<vec.size()-1;i++) {
            cout << vec[i] << " , ";
        }
        cout << vec[vec.size()-1] << " }\n";
    } else {
        cout << "ERROR: Vector has zero length.\n";
    }
}

string printIntVectorToString(vector<int> &vec) {
    
    string output;
    if(vec.size() > 0) {
        output += "{ ";
        for(unsigned int i=0;i<vec.size()-1;i++) {
            output += to_string(vec[i]);
		output += " , ";
        }
        output += to_string( vec[vec.size()-1]);
	output += " }";
    } else {
        cout << "ERROR: Vector has zero length.\n";
    }
    return output;
}


vector<int> multinomialSampling(int N, vector<double> p, const gsl_rng *r) {
    
    size_t dim = p.size();
    unsigned int n[dim];
    double* pPointer = &p[0];
    
    gsl_ran_multinomial(r,dim,N,pPointer,n);
    
    vector<int> result(n, n + sizeof n / sizeof n[0]);
    
    return(result);
    
}

int sumOfVector(vector<int> intVec) {
    
    int tot = 0;
    for(unsigned int i=0; i<intVec.size(); i++) {
        tot+= intVec[i];
    }
    return tot;
}

double sumOfVector(vector<double> doubleVec) {
    
    double tot = 0;
    for(unsigned int i=0; i<doubleVec.size(); i++) {
        tot+= doubleVec[i];
    }
    return tot;
}

vector<double> subtractVectorsDouble(vector<double> &a, vector<double> &b) {
    
    vector<double> result;
    //Assume equal length
    for(unsigned int i=0; i<a.size(); i++) {
        result.push_back(a[i]-b[i]);
    }
    return result;
}

vector<double> addVectorsDouble(vector<double> &a, vector<double> &b) {
    
    vector<double> result;
    //Assume equal length
    for(unsigned int i=0; i<a.size(); i++) {
        result.push_back(a[i]+b[i]);
    }
    return result;
}


double my_gsl_linalg_LU_det (gsl_matrix * LU, int signum)
{
    size_t i, n = LU->size1;
    
    double det = (double) signum;
    
    for (i = 0; i < n; i++)
    {
        det *= gsl_matrix_get (LU, i, i);
    }
    
    return det;
}

double determinant(gsl_matrix *A) {
    
    size_t dim = A->size1;
    gsl_permutation *p = gsl_permutation_alloc(dim);
    int signum;
    gsl_matrix * tmp_ptr = gsl_matrix_alloc(dim,dim);
    gsl_matrix_memcpy(tmp_ptr,A);
    gsl_linalg_LU_decomp(tmp_ptr,p,&signum);
    return gsl_linalg_LU_det(tmp_ptr,signum);
    
}

double logMultivariateNormalPDF(vector<double> &x, vector<double> &mu, gsl_matrix *sigma) {
    
    //Set up permutation and matrix copy
    int dim= (int) sigma->size1;
    vector<double> xMinusMu = subtractVectorsDouble(x,mu);
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

double logMultivariateNormalPDFPrint(vector<double> &x, vector<double> &mu, gsl_matrix *sigma) {
	
	cout << "Input mu:\n";
	printDoubleVector(mu);
	cout << "Input sigma:\n";
	printMatrixMathematica(sigma);    

    //Set up permutation and matrix copy
    int dim= (int) sigma->size1;
    vector<double> xMinusMu = subtractVectorsDouble(x,mu);
	cout << "xMinusMu: "; printDoubleVector(xMinusMu);

    gsl_permutation *p = gsl_permutation_calloc(dim);
    int signum;
    gsl_matrix * tmp_ptr = gsl_matrix_calloc(dim,dim);
    gsl_matrix_memcpy(tmp_ptr,sigma); //Copy elements of sigma into tmp_ptr
    
    //Get LU decomposition
    gsl_linalg_LU_decomp(tmp_ptr,p,&signum);
    
    //Get determinant
    double det = gsl_linalg_LU_det(tmp_ptr, signum);
	cout << "Determinant: " << det << "\n";
    
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
        return  -numeric_limits<double>::max();
    }
    
	cout << "Sigma inverse: "; printMatrix(sigmaInv);
    
    //double logPrefactor= -log(sqrt(pow(2*M_PI,dim)*det));
    double logPrefactor= -0.5*dim*log(2*M_PI) -0.5*log(det);
	cout << "logPrefactor: " << logPrefactor << "\n";
    
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
	cout << "dotProduct: " << dotProd << "\n";    

    //Clean up
    gsl_vector_free(xMinusMuGSL);
    gsl_vector_free(sigmaInvXminusMu);
    gsl_matrix_free(sigmaInv);
    gsl_matrix_free(tmp_ptr);
    gsl_permutation_free(p);
    
    return logPrefactor -0.5*dotProd;
}

//Integration parameters
struct parameters {

	double mu;
	double var;
	double stDev;
	vector<double> muVec;
	gsl_matrix* Sigma;

};

//Integrand for numerical integration with cubature
int integrandMultiNormal(unsigned dim, const double *x, void *fdata, unsigned int fdim, double *fval) {

//	count++; //Increase number of evaluations
	parameters* params = (parameters*) (fdata); //Cast fval to a parameters object			
	vector<double> xvals(dim);
	for(unsigned int i=0; i<dim; i++) { xvals[i] = x[i]; }

	*fval = exp(logMultivariateNormalPDF(xvals, params->muVec , params->Sigma));
	
	return  0;
}

//Integrand for (one-dimensional) normal numerical integration with cubature
int integrandNormal(unsigned dim, const double *x, void *fdata, unsigned int fdim, double *fval) {

//	count++; //Increase number of evaluations
	parameters* params = (parameters*) (fdata); //Cast fval to a parameters object			


	*fval = exp(-(x[0]-params->mu)*(x[0]-params->mu)/(2*params->var) - log(params->stDev*sqrt(2*M_PI)));

	return  0;
}

void printMatrix(gsl_matrix *A) {
    
    int dim1= (int) A->size1;
    int dim2= (int) A->size2;
    
    for (int i = 0; i < dim1; i++)  /* OUT OF RANGE ERROR */
        for (int j = 0; j < dim2; j++)
            printf ("m(%d,%d) = %g\n", i, j, gsl_matrix_get (A, i, j));
    
    
}

//Outputs matrix to screen in mathematica format A = {{A11,A12,A13},{A21,A22,A23},{A31,A32,A33}}
void printMatrixMathematica(gsl_matrix *A) {
    
    int dim1= (int) A->size1;
    int dim2= (int) A->size2;
    

	cout << "{";
	for (int i = 0; i < dim1; i++) {
		cout << "{";
		for (int j = 0; j < dim2; j++) {
			cout.precision(25);
			cout << gsl_matrix_get(A,i,j);
			cout.precision(6); //Back to default
			if(j<dim2-1) {
				cout << ",";
			}
  		}
		cout << "}";
		if(i<dim1-1) {
			cout << ",";
		}
	}
	cout << "}\n"; 
}

string printMatrixMathematicaToString(gsl_matrix* A) {
    
	stringstream output;
	int dim1= (int) A->size1;
	int dim2= (int) A->size2;   

	output << "{";
	for (int i = 0; i < dim1; i++) {
		output << "{";
		for (int j = 0; j < dim2; j++) {
			output << gsl_matrix_get(A,i,j);
			if(j<dim2-1) {
				output << ",";
			}
   		}
		output << "}";
		if(i<dim1-1) {
			output << ",";
		}
	}
	output << "}";

	return output.str();
}

//Create a total of number distinct haplotypes of length length
vector<Sequence> generateRandomFullHaplotypes(int length, int number, gsl_rng *r) {
    
    vector<Sequence> fhSeqs;
    Sequence temp;
    for(int i=0;i<number;i++) {
        fhSeqs.push_back(temp);
        for(int j=0;j<length;j++) {
            fhSeqs[i].addBase('-');
        }
    }
    if(number > pow(4,length)) {
        cout << "Number ("<<number<<") too large. Can't create that many distinct haplotypes with a length of only "<<length<<".\n";
        return fhSeqs; //returns -1 at all bases
    }
    
    //A=0,C=1,G=2,T=3
    fhSeqs[0]=randomSequence(length,r);
    for(int i=1;i<number;i++) {
        Sequence temp;
        bool unique = false;
        while(unique==false) {
            unique=true;
            temp = randomSequence(length, r);
            for(int j=0;j<i;j++) {
                if(identicalSeqs(temp,fhSeqs[j]) == true) { unique = false; break; }
            }
        }
        fhSeqs[i]=temp;
    }
    
    return fhSeqs;
}


Sequence randomSequence(int length, gsl_rng *r) {
    
    Sequence result;
    for(int i=0;i<length;i++) {
        int temp = (int) gsl_rng_uniform_int(r,4);
        if(temp==0) {
            result.addBase('A');
        } else if (temp==1) {
            result.addBase('C');
        } else if (temp==2) {
            result.addBase('G');
        } else {
            result.addBase('T');
        }
    }
    return result;
}

vector<Sequence> generateMajorMinorFullHaplotypes(int length, int number, gsl_rng *r) {
    
    vector<Sequence> haps;
    for(int i=0;i<number;i++) { haps.push_back(Sequence(length)); }
    
    if(number > pow(2,length)) {
        cout << "Number ("<<number<<") too large. Cannot make that many distinct haplotypes of length "<<length<<".\n";
        return haps;
    }
    
    for(int i=0;i<length;i++) {
        
        //Find major allele for locus i
        int majorInt = (int) gsl_rng_uniform_int(r,4);
        char major;
        if(majorInt==0) { major = 'A'; }
        else if(majorInt==1) { major = 'C'; }
        else if(majorInt==2) { major = 'G'; }
        else { major = 'T'; }
        
        //Find minor allele
        bool uniqueAllele = false;
        char minor = '_';
        while(uniqueAllele == false) {
            int minorInt = (int) gsl_rng_uniform_int(r,4);
            if(minorInt==0) { minor = 'A'; }
            else if(minorInt==1) { minor = 'C'; }
            else if(minorInt==2) { minor = 'G'; }
            else { minor = 'T'; }
            if(minor != major) { uniqueAllele = true; }
        }
        
        double MAF = 0.3;
        for(int j=0;j<number;j++) {
            if(gsl_rng_uniform(r) > MAF) { haps[j].setBase(i,major); }
            else { haps[j].setBase(i,minor); }
        }
    }
    //Printing
    //for(int i=0;i<number;i++) { haps[i].print(); }
    //cout << length << " " << number << "\n";
    
    return haps;
    
}

bool identicalSeqs(Sequence & a, Sequence & b) {
    
    int lengthA = a.getLength();
    if(lengthA != b.getLength()) {
        cout << "ERROR in identicalSeqs: Sequences a and b are of different lengths!!\n";
        cout << "a = ";
        a.print();
        cout << "b = ";
        b.print();
        cout << "Returning false.\n";
        return false;
    } else {
        
        for(int i=0;i<lengthA;i++) {
            if(a.getBase(i) != b.getBase(i)) {
                return false;
            }
        }
        return true;
    }
    
}

bool identicalIntVectors(std::vector<int> & a, std::vector<int> & b) {

	int lengthA = (int) a.size();
	
	if(lengthA != (int) b.size()) {
		cout << "ERROR in identicalIntVectors: vectors a and b are of different lengths!!\n";
		cout << "a = ";
		printIntVector(a);
		cout << "b = ";
		printIntVector(b);
        cout << "Returning false.\n";
        return false;
    } else {
        
		for(int i=0;i<lengthA;i++) {
            if(a[i] != b[i]) {
                return false;
            }
        }
        return true;
    }
   
}

bool nonZeroIntVector(vector<int> &a) {
    
    for(unsigned int i=0;i<a.size();i++) {
        
        if(a[i] == 0) { return false; }
    }
    
    return true;
}

bool nonZeroDoubleVector(vector<double> &a) {
    
    for(unsigned int i=0;i<a.size();i++) {
        
        if(a[i] == 0) { return false; }
    }
    
    return true;
}

void setupRNG(gsl_rng *r, unsigned long int seed) {
    
    const gsl_rng_type * T; //ISSUE HERE: Cannot define T here. Will be deleted on exit
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r,seed);
    
    cout << "RNG set up.n";
}

//Generates diploid sequence where allele1 and allele2 differ from each other at all positions
DiploidSequence generateRandomDiploidSequence(int length, gsl_rng *r) {
    
    Sequence major = randomSequence(length,r);
    Sequence minor;
    
    DiploidSequence result;
    result.setMajor(major);
    
    for(int i=0;i<length;i++) {
        char current = '-';
        bool alleleDifferentFromMajor = false;
        
        while(alleleDifferentFromMajor==false) {
            int temp = (int) gsl_rng_uniform_int(r,4);
            if(temp==0) {
                current = 'A';
            } else if (temp==1) {
                current = 'C';
            } else if (temp==2) {
                current = 'G';
            } else {
                current = 'T';
            }
            
            if(current != major.getBase(i)) {
                alleleDifferentFromMajor = true;
            }
        }
        minor.addBase(current);
    }
    
    result.setMinor(minor);
    
    return result;
    
}

//Generates diploid sequence where allele1 and allele2 differ from each other at all positions
//Constraint in the form of sequence, e.g. A--T-
DiploidSequence generateRandomSemiConstrainedDiploidSequence(Sequence constraint, gsl_rng *r) {
    
	int length = constraint.getLength();
    Sequence major = randomSequence(length,r);
    Sequence minor;
    
    for(int i=0;i<length;i++) { //Loop over sequence
	
	//Check if major allele need to be updated according to constraint
	if(constraint.getBase(i) != '-') {
		major.setBase(i, constraint.getBase(i));
	}

        char current = '-';
        bool alleleDifferentFromMajor = false;
        
        while(alleleDifferentFromMajor==false) {
            int temp = (int) gsl_rng_uniform_int(r,4);
            if(temp==0) {
                current = 'A';
            } else if (temp==1) {
                current = 'C';
            } else if (temp==2) {
                current = 'G';
            } else {
                current = 'T';
            }
            
            if(current != major.getBase(i)) {
                alleleDifferentFromMajor = true;
            }
        }
        minor.addBase(current);
    }
    
 	DiploidSequence result;
	result.setMajor(major);
    	result.setMinor(minor);
    
    return result;
    
}


//Generates diploid sequence where allele1 and allele2 differ from each other at all positions
//Constraint in the form of sequence, e.g. A/C,-/-,-/-,T/G,-/-
DiploidSequence generateRandomSemiConstrainedDiploidSequence(DiploidSequence constraint, gsl_rng *r) {
    
	int length = constraint.getLength();
	Sequence major = randomSequence(length,r); //ACGTACAGAT (random)
	Sequence minor = Sequence(length); //--------- (undefined)
    
	for(int i=0;i<length;i++) { //Loop over sequence
	
		//Check if locus need to be updated according to constraint
		bool updateMinor = false;
		//Fully constrained case (when both major and minor alleles must be specified)
		if(constraint.getMajor().getBase(i) != '-' && constraint.getMinor().getBase(i) != '-') {
			major.setBase(i, constraint.getMajor().getBase(i));
			minor.setBase(i, constraint.getMinor().getBase(i));
	
		} else if(constraint.getMajor().getBase(i) != '-') { //Semi constraint case (when the minor allele doesn't matter)

			major.setBase(i, constraint.getMajor().getBase(i));
			updateMinor = true;
		
		} else { //Unconstrained case
			updateMinor = true;
		}			


		if(updateMinor == true) {//Update minor allele to be random, but different from major allele

			char current = '-';
			bool alleleDifferentFromMajor = false;
        
			while(alleleDifferentFromMajor==false) {
				int temp = (int) gsl_rng_uniform_int(r,4);
				if(temp==0) {
					current = 'A';
				} else if (temp==1) {
					current = 'C';
				} else if (temp==2) {
					current = 'G';
				} else {
					current = 'T';
				}
            
				if(current != major.getBase(i)) {
					alleleDifferentFromMajor = true;
				}
			}
			minor.setBase(i, current);
		}
	}
    
 	DiploidSequence result;
	result.setMajor(major);
    	result.setMinor(minor);
    
	return result;
    
}

vector<Sequence> generateFullHapsFromAlleles(DiploidSequence &fullHapAlleles) {
    
    //Original data
    vector<Sequence> fullHaps;
    int lengthFH = fullHapAlleles.getMajor().getLength();
    Sequence major = fullHapAlleles.getMajor();
    Sequence minor = fullHapAlleles.getMinor();

	cout << "In generateFullHapsFromAlleles. Major: "; major.print();
	cout << "Minor: "; minor.print();
	cout << "lengthFH: " << lengthFH << "\n";
    
    //Restructured data
    vector<vector<char> > fh;
    for(int i=0;i<lengthFH;i++) {
        fh.push_back({major.getBase(i), minor.getBase(i)});
    }

	for(int n=0; n<lengthFH; n++) { //Loop over loci

		if(n==0) {

			//Create the 0th locus alleles
			Sequence seq1, seq2;
			seq1.addBase(fh[0][0]);
			seq2.addBase(fh[0][1]);
			fullHaps.push_back(seq1);
			fullHaps.push_back(seq2);

		} else {

			//For nth locus, duplicate current sequences and add either fh[n][0] or fh[n][1] to the end
			vector<Sequence> fullHapsDup = fullHaps;

			for(unsigned int i=0; i<fullHaps.size(); i++) {

				fullHaps[i].addBase(fh[n][0]);
				fullHapsDup[i].addBase(fh[n][1]);
			}

			//Add fullHapsDup to fullHaps
			for(unsigned int i=0; i<fullHapsDup.size(); i++) {

				fullHaps.push_back(fullHapsDup[i]);
			}
		}


		//cout << "Printing haps at step " << n << " of " << lengthFH-1 << ":\n";
		for(unsigned int i=0; i<fullHaps.size(); i++) {

			fullHaps[i].print();
		}

	}

        return fullHaps;
}


void randomiseSequenceVector(vector<Sequence> & s, gsl_rng *r) {
    
    int length = (int) s.size();
    int numSwaps = length * 10;
    
    for(int i=0; i<numSwaps;i++) {
        
        int swap1 = gsl_rng_uniform_int(r,length);
        int swap2 = gsl_rng_uniform_int(r,length);
        
        Sequence seq1 = s[swap1];
        Sequence seq2 = s[swap2];
        
        s[swap1] = seq2;
        s[swap2] = seq1;
        
    }
    
}

//Picks length random (but not identical) ints in the range [0,max)
vector<int> pickRandomInts(gsl_rng *r, int max, int length) {
    
    vector<int> result;
    
    if(max >= length) {
        for(int i=0;i<length;i++) {
            
            int pick;
            bool unique = false;
            while(unique == false) {
                pick = gsl_rng_uniform_int(r,max);
                
                unique = true; //Assume it is true, then check if it actually is
                for(unsigned int j=0; j<result.size(); j++) {
                    if(result[j] == pick) { unique = false; break; }
                }
            }
            result.push_back(pick);
        }
    } else {
        cout << "ERROR IN pickRandomInts. max < length\n";
    }
    
    return result;
}


//Method for adding a random small delta to a random frequency index
void addRandom(vector<double> &oldVec, vector<double> &newVec, double delta, gsl_rng *r) {
    
    //Set newFreqs to oldFreqs values
    newVec = oldVec;

	//Update index of newFreq
	int index = gsl_rng_uniform_int(r,oldVec.size());
	double newValue = newVec[index] + delta*(gsl_rng_uniform_pos(r) - 0.5);
	if(newValue < 1e-11) {
		newValue = 1e-10;
	}
	newVec[index] = newValue;

	rescale(newVec);
}

//Method for adding a random small delta to a random frequency index
//None of the entries may be less than min
void addRandomMin(vector<double> &oldVec, vector<double> &newVec, double delta, gsl_rng *r, double min) {
    
    //Set newFreqs to oldFreqs values
    newVec = oldVec;

	//Update index of newFreq
	int index = gsl_rng_uniform_int(r,oldVec.size());
	double newValue = newVec[index] + delta*(gsl_rng_uniform_pos(r) - 0.5);
	if(newValue < min) {
		newValue = 1.5*min;
	}
	newVec[index] = newValue;

	rescaleMin(newVec, min);
}

//Method for adding a random small delta to a random frequency index
void addRandom(vector<double> & vec, double delta, gsl_rng *r) {
    
	//Update index of newFreq
	int index = gsl_rng_uniform_int(r,vec.size());
	double newValue = vec[index] + delta*(gsl_rng_uniform_pos(r) - 0.5);
	if(newValue < 1e-11) {
		newValue = 1e-10;
	}
	vec[index] = newValue;

	rescale(vec);
}

double logMultinomialProb(vector<double> & freq, vector<int> & obs) {
    
    double freqArray[(int) freq.size()];
    copy(freq.begin(), freq.end(), freqArray);
    
    unsigned int obsArray[(int) obs.size()];
    copy(obs.begin(), obs.end(), obsArray);
    
    return gsl_ran_multinomial_lnpdf(freq.size(), freqArray, obsArray);
}

//Store the log of the factorials of 0 to N in fact_store
void findLogFactorial(vector<double>& fact_store,int N){
    double logN=0;
    fact_store.push_back(0);
    for (int i=1;i<=N;i++) {
        logN=logN+log(i);
        fact_store.push_back(logN);
    }
}

//Dirichlet multinomial function taking vectors of observations and inferred frequencies as inputs
double logDirMultProbC(double C, vector<int> &obs, vector<double> &inf, vector<double> &fact_store) {
    
    int N = sumOfVector(obs);
    double bin = 0;
    
    if (N>0) {
        bin=fact_store[N]; //logN! == lnGamma(N+1)
        for (unsigned int i=0;i<obs.size();i++) {
            bin=bin-fact_store[obs[i]];
        }
        vector<double> alpha;
        for (unsigned int i=0;i<inf.size();i++) {
            alpha.push_back(C*inf[i]);
        }
        double a=0;
        for (unsigned int i=0;i<alpha.size();i++) {
            a=a+alpha[i];
            bin=bin-gsl_sf_lngamma(alpha[i]);
        }
        bin=bin+gsl_sf_lngamma(a);
        a=0;
        for (unsigned int i=0;i<alpha.size();i++) {
            double b=alpha[i]+obs[i];
            a=a+b;
            bin=bin+gsl_sf_lngamma(b);
        }
        bin=bin-gsl_sf_lngamma(a);
    } else {
        bin=0;
    }
    
    return(bin);
}

//Dirichlet multinomial function taking vectors of observations and inferred frequencies as inputs
double logDirMultProb(vector<double> &alpha, vector<int> &x, int& n) {
    

	double sumAlpha = 0;
	double sum = 0;
	for(unsigned int i=0; i<alpha.size(); i++) {
		sumAlpha += alpha[i];
		sum += gsl_sf_lngamma(x[i]+alpha[i]) - gsl_sf_lngamma(x[i]+1) - gsl_sf_lngamma(alpha[i]);
	}
	double logL = gsl_sf_lngamma(n+1) + gsl_sf_lngamma(sumAlpha) - gsl_sf_lngamma(n+sumAlpha) + sum;

	return(logL);
}

double logBetaBinProb(double alpha, double beta, int x, int n) {


	return(gsl_sf_lngamma(n+1)+gsl_sf_lngamma(x+alpha)+gsl_sf_lngamma(n-x+beta)+gsl_sf_lngamma(alpha+beta)-gsl_sf_lngamma(x+1)-gsl_sf_lngamma(n-x+1)-gsl_sf_lngamma(n+alpha+beta)-gsl_sf_lngamma(alpha)-gsl_sf_lngamma(beta));
}

//Simple computation of CDF for beta binomial. Will be slow for large threshold and n.
double BetaBinCDF(double alpha, double beta, double threshold, int n) {

	if(threshold < 0) { 
		return 0;
	} else if(threshold >= n) {
		 return 1;
	} else {

		double prob = 0;
		for(int x=0; x<=threshold; x++) {

			prob += exp(logBetaBinProb(alpha,beta,x,n));
		}
		return prob;
	}
}


//Check if a file exists. Give full path to file
bool fileExists (string& fileName) {
    ifstream f(fileName.c_str());
    if (f.good()) {
        f.close();
        return true;
    } else {
        f.close();
        return false;
    }
}

//Mean and var are reduced by 1 dimension (namely the last (kth) dimension)
vector<double> computeAlpha(vector<double>& mean, gsl_matrix* var, int& n) {


	//Simple alpha_0 created from the 0,0 component only
	//q_i,j = var_i,j/(mean_i * (delta_i,j - mean_j/n))
//	double q_00 = gsl_matrix_get(var,0,0)/(mean[0]*(1-mean[0]/n));
//	cout << "q_00: " << q_00 << "\n";
//	double alpha_0 = (n-q_00)/(q_00-1);
//	cout << "alpha_0: " << alpha_0 << "\n";

	//Alpha_0 created from the average of all alpha_0_i,j
	double alpha_0_sum = 0;
	int numIncluded = 0;
	for(unsigned int i=0; i<mean.size(); i++) {
		for(unsigned int j=0; j<mean.size(); j++) {

			int KroneckerDelta = 0;
			if(i==j) { KroneckerDelta= 1; }
			double q_ij = gsl_matrix_get(var,i,j)/(mean[i]*(KroneckerDelta-mean[j]/n));
			double alpha_0_ij = (n-q_ij)/(q_ij-1);
//			cout << "alpha_o_" << i << "," << j <<": " << alpha_0_ij << "\n";
			if(alpha_0_ij > 0) {
				alpha_0_sum += alpha_0_ij;
				numIncluded++;
			}
		}
	}
	cout << "alpha_0_sum: " << alpha_0_sum << "\n";
	//double alpha_0 = alpha_0_sum/(mean.size()*mean.size());
	double alpha_0 = alpha_0_sum/numIncluded;
	cout << "alpha_0: " << alpha_0 << "\n";
//	if(alpha_0<=0) { alpha_0 = 0.01; }
	if(alpha_0_sum == 0) {
	
		for(unsigned int i=0; i<mean.size(); i++) {
			for(unsigned int j=0; j<mean.size(); j++) {

				int KroneckerDelta = 0;
				if(i==j) { KroneckerDelta= 1; }
				cout << "var_" << i << "_" << j << ": " << gsl_matrix_get(var,i,j) << "\n";
				cout << "mean_" << i << ": " << mean[i] << "\n";
				cout << "mean_" << j << ": " << mean[j] << "\n";
				double q_ij = gsl_matrix_get(var,i,j)/(mean[i]*(KroneckerDelta-mean[j]/n));
				double alpha_0_ij = (n-q_ij)/(q_ij-1);
				cout << "alpha_o_" << i << "," << j <<": " << alpha_0_ij << "\n";
				if(alpha_0_ij > 0) {
					alpha_0_sum += alpha_0_ij;
					numIncluded++;
				}
			}
		}

	}

	
	//Create full dimensional alpha
	vector<double> alpha;
	double sumMean = 0;
	for(unsigned int i=0; i<mean.size(); i++) {

		alpha.push_back(alpha_0*mean[i]/n);
		sumMean += mean[i];
	}
	//Do kth dimension
	alpha.push_back(alpha_0*(n-sumMean)/n);
	cout << "alpha: "; printDoubleVector(alpha);
	cout << "mean: "; printDoubleVector(mean);
	cout << "sumMean: " << sumMean << " n: " << n << "\n";
	
	return(alpha);

}

vector<int> DirMultSampling(int N, vector<double> &freqs, double C, const gsl_rng *r) {
    
    //Define the prior, alpha, for the Dirichlet distribution
    double alpha[freqs.size()];
    for(unsigned int i=0;i<freqs.size();i++) {
        alpha[i] = C*freqs[i];
    }
    
    //Sample frequencies p from the Dirichlet distribution
    double p[freqs.size()]; //Placeholder
    gsl_ran_dirichlet(r, freqs.size(), alpha, p);
    
    //Convert p to a vector
    vector<double> pVec;
    for(unsigned int i=0;i<freqs.size();i++) {
        pVec.push_back(p[i]);
    }
    
    //Return multinomial sample based on Dirichlet prior
    return multinomialSampling(N,pVec,r);
}



double computeLSelection(vector<int> &obsB, vector<int> &obsA, int Nt, AnaParam* ap, vector<double> &qBfhFreqs, vector<double> &hapFitT, vector<vector<int> >* contribsPHset, bool print) {

	double C = ap->getC();
	vector<double>& factStore = *(ap->getLogFactStore()); //vector is not copied here! Saves much time
    

    
	if(print == true) {
	
		cout << "Printing intermediate steps in computeLselection.\n";
		cout << "First, the input variables:\n";
		cout << "obsB: "; printIntVector(obsB);
		cout << "obsA: "; printIntVector(obsA);
		cout << "Nt: " << Nt << "\n";
		cout << "C: " << C << "\n";
		cout << "qBfhFreqs: "; printDoubleVector(qBfhFreqs);
		cout << "hapFitT: "; printDoubleVector(hapFitT);
		cout << "Contribs: \n";
		for(unsigned int i=0; i<contribsPHset->size(); i++) {
			cout << "PH " << i << ": "; printIntVector((*contribsPHset)[i]);
		}


	}

	//For now, disregard sigmaG
	vector<double> hapFitG(hapFitT.size(), 0); //Zero fitness

	//Dimensions
	int dimFH = (int) qBfhFreqs.size();
	int dimPH = (int) obsB.size();
	

	//Total observations
	int Na = sumOfVector(obsA);
	
	if(print == true) {
		cout << "dimFH: " << dimFH << "\n";
		cout << "dimPH: " << dimPH << "\n";
		cout << "Na: " << Na << "\n";
	}

	//Construct matrix T
	gsl_matrix *T = gsl_matrix_alloc(dimPH,dimFH);
	for(int i=0; i<dimPH; i++) { //Loop over rows

		int counter=0;
		for(int j=0; j<dimFH; j++) { //Loop over columns

			if((*contribsPHset)[i][counter]==j && counter < (int) (*contribsPHset)[i].size()) { //This assumes the contribs are ordered
				gsl_matrix_set(T,i,j,1);
				counter++;
			} else {
				gsl_matrix_set(T,i,j,0);
			}
		}
	}

	if(print == true) {
		cout << "Matrix T: "; printMatrixMathematica(T);
	}


	/*
	* * 
	* * Before computation
	* *
	*/

	/*
	 * First, transform qB into qB_PH=TqB
	 */
	//Turn qB into gsl vector
	gsl_vector *qBGSL = gsl_vector_alloc(dimFH);
	for(int i=0;i<dimFH;i++) {
		gsl_vector_set(qBGSL,i,qBfhFreqs[i]);
	}

	//Create TqB through matrix vector multiplication
	gsl_vector *TqBGSL = gsl_vector_alloc(dimPH);
	gsl_blas_dgemv(CblasNoTrans, 1.0, T, qBGSL, 0.0, TqBGSL);

	//Convert back into C++ vector
	vector<double> TqB;
	for(int i=0; i<dimPH; i++) {
		TqB.push_back(gsl_vector_get(TqBGSL, i));
	}  

	if(print == true) {
		cout << "TqB: "; printDoubleVector(TqB);
	}
	
	
	/*
	* Finally, compute log likelihood from Dirichet multinomial PDF with dispersion parameter C
	*/
	double logPxBph = logDirMultProbC(C, obsB, TqB, factStore);
	
	if(print == true) {
		cout << "logPxBph = " << logPxBph << "\n";
	}



	/*
	* * 
	* * After computation
	* *
	*/


	/*
 	* First we create the mean of xAPH = Na*T*sigmaG[sigmaT[qB]]
 	*/ 
	vector<double> sigmaTqB;
	double totalT = 0;
	for(int i=0; i<dimFH; i++) {
		totalT += qBfhFreqs[i]*exp(hapFitT[i]);
	}
	for(int i=0; i<dimFH; i++) {
		sigmaTqB.push_back(qBfhFreqs[i]*exp(hapFitT[i])/totalT);
	}

	vector<double> sigmaGsigmaTqB;
	double totalG = 0;
	for(int i=0; i<dimFH; i++) {
		totalG += sigmaTqB[i]*exp(hapFitG[i]);
	}
	for(int i=0; i<dimFH; i++) {
		sigmaGsigmaTqB.push_back(sigmaTqB[i]*exp(hapFitG[i])/totalG);
	}

	if(print == true) {
		cout << "sigmaTqB: "; printDoubleVector(sigmaTqB);
		cout << "sigmaGsigmaTqB: "; printDoubleVector(sigmaGsigmaTqB);
	}

	gsl_vector *sigmaGsigmaTqBGSL = gsl_vector_alloc(dimFH);
	for(int i=0;i<dimFH;i++) {
		gsl_vector_set(sigmaGsigmaTqBGSL,i,sigmaGsigmaTqB[i]);
	}

	//Create TsigmaGsigmaTqB through matrix vector multiplication
	gsl_vector *TsigmaGsigmaTqBGSL = gsl_vector_alloc(dimPH);
	gsl_blas_dgemv(CblasNoTrans, 1.0, T, sigmaGsigmaTqBGSL, 0.0, TsigmaGsigmaTqBGSL);  

	//Convert back to C++ vector
	vector<double> meanXaPH;
	vector<double> TsigmaGsigmaTqB; //To be used in computing variance
	for(int i=0; i<dimPH; i++) {
		meanXaPH.push_back(Na*(gsl_vector_get(TsigmaGsigmaTqBGSL, i)));
		TsigmaGsigmaTqB.push_back(gsl_vector_get(TsigmaGsigmaTqBGSL, i));
	}

	if(print == true) {
		cout << "meanXaPH: "; printDoubleVector(meanXaPH);
	}


	/*
	* Next we compute the variance of xAPH
	*/
	//Starting with the first term: alpha*Na*M(TsigmaGsigmaTqB)
	//Constructing the matrix M
	gsl_matrix *MTsigmaGsigmaTqB = gsl_matrix_alloc(dimPH,dimPH);
	constructMatrixM(TsigmaGsigmaTqB, MTsigmaGsigmaTqB);
	
	//Multiplying by the factor alpha*Na:
	double alpha = (Na+C)/((double)(1+C));
	gsl_matrix_scale(MTsigmaGsigmaTqB, Na*alpha);

	if(print==true) {
		cout << "Na*alpha*MTsigmaGsigmaTqB: "; printMatrixMathematica(MTsigmaGsigmaTqB);
	}



	//Second term: Na(Na-alpha)betaTDsigmaGsigmaTqBMsigmaTqB(DsigmaGsigmaTqB)transTtrans
	//Create MsigmaTqB first
	gsl_matrix *MsigmaTqB = gsl_matrix_alloc(dimFH, dimFH);
	constructMatrixM(sigmaTqB, MsigmaTqB);

	if(print==true) {
		cout << "MsigmaTqB: "; printMatrixMathematica(MsigmaTqB);
	}

	//Compute DsigmaGsigmaTqB
	gsl_matrix *DsigmaGsigmaTqB = gsl_matrix_alloc(dimFH,dimFH);
	for(int i=0; i<dimFH; i++) {
		for(int j=0; j<dimFH; j++) {
            
			if(i==j) {
				gsl_matrix_set(DsigmaGsigmaTqB, i,j, exp(hapFitG[i])/totalG - sigmaTqB[i]*exp(2*hapFitG[i])/(totalG*totalG));
			} else {
				gsl_matrix_set(DsigmaGsigmaTqB, i,j, -sigmaTqB[i]*exp(hapFitG[i]+hapFitG[j])/(totalG*totalG));
			}
		}
	}
	
	if(print==true) {
		cout << "DsigmaGsigmaTqB"; printMatrixMathematica(DsigmaGsigmaTqB);
	}


	//Creating constant factor Na(Na-alpha)beta
	int growthFactor = ap->getGrowthFactor();
	int Ng=growthFactor*Nt;
	double NtNgDouble = ((double)Nt)*((double)Ng); //This computation gets out of integer range if Nt>4634, so turn into double
	double beta = (Nt+Ng-1)/NtNgDouble; //Ng not yet defined. Set equal to Nt or Nt*100? Not quite sure.
	double NaNaMinusAlphaBeta = Na*(Na-alpha)*beta;



	//Compute the matrix multiplications
	gsl_matrix *TDsigmaGsigmaTqB = gsl_matrix_alloc(dimPH,dimFH);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, T, DsigmaGsigmaTqB, 0.0, TDsigmaGsigmaTqB);
	gsl_matrix *TDsigmaGsigmaTqBMsigmaTqB = gsl_matrix_alloc(dimPH,dimFH);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, TDsigmaGsigmaTqB, MsigmaTqB, 0.0, TDsigmaGsigmaTqBMsigmaTqB);
	gsl_matrix *TDsigmaGsigmaTqBMsigmaTqBDsigmaGsigmaTqBTransTTrans = gsl_matrix_alloc(dimPH,dimPH);
	gsl_blas_dgemm (CblasNoTrans, CblasTrans, NaNaMinusAlphaBeta, TDsigmaGsigmaTqBMsigmaTqB, TDsigmaGsigmaTqB, 0.0, TDsigmaGsigmaTqBMsigmaTqBDsigmaGsigmaTqBTransTTrans);

	if(print==true) {
		cout << "Na(Na-alpha)*beta*TDsigmaGsigmaTqBMsigmaTqBDsigmaGsigmaTqBTransTtrans"; printMatrixMathematica(TDsigmaGsigmaTqBMsigmaTqBDsigmaGsigmaTqBTransTTrans);
	}


	//Create the final variance matrix (stored in MTsigmaGsigmaTqB with the other matrix unchanged)
	gsl_matrix_add(MTsigmaGsigmaTqB, TDsigmaGsigmaTqBMsigmaTqBDsigmaGsigmaTqBTransTTrans);

	if(print == true) {
		cout << "varXaPH (full dim): "; printMatrixMathematica(MTsigmaGsigmaTqB);
	}



	//Reduce dimensionality by 1 to ensure non-degeneracy
	meanXaPH.pop_back(); //WLOG remove the last PH
	gsl_matrix_view varXaPHreducedView = gsl_matrix_submatrix(MTsigmaGsigmaTqB, 0,0, dimPH-1, dimPH-1); //Remove row k and column k
	gsl_matrix * varXaPH = gsl_matrix_alloc(dimPH-1, dimPH-1);
	gsl_matrix_memcpy(varXaPH, &(varXaPHreducedView.matrix)); //Copy the reduced matrix view into a "normal" matrix varXaPH
	vector<double> obsAreduced;
	for(int i=0; i<dimPH-1; i++) {
		obsAreduced.push_back(obsA[i]);
	}

	if(print == true) {
		cout << "Reduced quantities: \n";
		cout << "ObsAreduced: "; printDoubleVector(obsAreduced);
		cout << "meanXaPHreduced: "; printDoubleVector(meanXaPH);
		cout << "varXaPHreduced: "; printMatrixMathematica(varXaPH);
	}

	double logPxAph = -numeric_limits<double>::max(); //Set to neg infinity as precaution

	if(dimPH > 2) { //Multivariate approach	

		//Compute likelihood
		logPxAph = logMultivariateNormalPDF(obsAreduced,meanXaPH,varXaPH);
    
		
        
	}else if(dimPH == 2) {
	
		//In this case we retrieve the univariate case when reducing dimensionality by one
		double meanXaPHuni = meanXaPH[0];
		double varXaPHuni = gsl_matrix_get(varXaPH,0,0);
		double stDevXaPHuni = sqrt(varXaPHuni);

		logPxAph = -(obsA[0]-meanXaPHuni)*(obsA[0]-meanXaPHuni)/(2*varXaPHuni) - log(stDevXaPHuni*sqrt(2*M_PI)); //log L


	} else { //error
	
		cout << "ERROR in Misc::computeLSelection: dimPH is <2\n";
		logPxAph = -numeric_limits<double>::max();
	}

	if((logPxAph >= 0 || logPxAph == -numeric_limits<double>::max()) && print == false) {
		cout << "ERROR in Misc::computeLSelection: Likelihood is positive or negative infinity. Likelihood is: " << logPxAph  << "\n";

		if(dimPH>2) {
			cout << "Recomputing likelihood whilst printing out stats:\n";
			logMultivariateNormalPDFPrint(obsAreduced,meanXaPH,varXaPH);
		}

		cout << "-------------------------------------------------------\n";
		cout << "Recomputing mean and variance whilst printing out stats:\n";
		computeLSelection(obsB,obsA,Nt,ap,qBfhFreqs,hapFitT, contribsPHset, true);

		//Clean up, then return negative infinity
		gsl_vector_free(qBGSL);
		gsl_vector_free(TqBGSL);
		gsl_vector_free(sigmaGsigmaTqBGSL);
		gsl_vector_free(TsigmaGsigmaTqBGSL);

		gsl_matrix_free(T);
		gsl_matrix_free(MTsigmaGsigmaTqB);
		gsl_matrix_free(MsigmaTqB);
		gsl_matrix_free(DsigmaGsigmaTqB);
		gsl_matrix_free(TDsigmaGsigmaTqB);
		gsl_matrix_free(TDsigmaGsigmaTqBMsigmaTqB);
		gsl_matrix_free(TDsigmaGsigmaTqBMsigmaTqBDsigmaGsigmaTqBTransTTrans);
		gsl_matrix_free(varXaPH);


		cout << "Returning likelihood of negative infinity.\n\n";
		return -numeric_limits<double>::max();
	}

	//Clean up before exit
	gsl_vector_free(qBGSL);
	gsl_vector_free(TqBGSL);
	gsl_vector_free(sigmaGsigmaTqBGSL);
	gsl_vector_free(TsigmaGsigmaTqBGSL);

	gsl_matrix_free(T);
	gsl_matrix_free(MTsigmaGsigmaTqB);
	gsl_matrix_free(MsigmaTqB);
	gsl_matrix_free(DsigmaGsigmaTqB);
	gsl_matrix_free(TDsigmaGsigmaTqB);
	gsl_matrix_free(TDsigmaGsigmaTqBMsigmaTqB);
	gsl_matrix_free(TDsigmaGsigmaTqBMsigmaTqBDsigmaGsigmaTqBTransTTrans);
	gsl_matrix_free(varXaPH);

	//cout << "logPxAph = " << logPxAph << "\n";
	return logPxBph + logPxAph;
}

double computeLSelTVar(vector<int> &obsA, int Nt, AnaParam *ap, vector<double> &qBmean, gsl_matrix* qBvar, vector<double>& hapFitT, vector<vector<int> >* contribsPHset, bool print) {

	double C = ap->getC();
	double logPxAph = -numeric_limits<double>::max(); //Set to neg infinity as precaution

	if(print == true) {
	
		cout << "Printing intermediate steps in computeLselection.\n";
		cout << "First, the input variables:\n";
		cout << "obsA: "; printIntVector(obsA);
		cout << "Nt: " << Nt << "\n";
		cout << "C: " << C << "\n";
		cout << "qBmean: "; printDoubleVector(qBmean);
		cout << "hapFitT: "; printDoubleVector(hapFitT);
		cout << "Contribs: \n";
		for(unsigned int i=0; i<contribsPHset->size(); i++) {
			cout << "PH " << i << ": "; printIntVector((*contribsPHset)[i]);
		}


	}

	//Dimensions
	int dimFH = (int) qBmean.size();
	int dimPH = (int) obsA.size();
	int dimPHminusOne = -1;
	if(dimPH==1) { dimPHminusOne = 1; } //If dimPH==1, this corresponds to at least dimPH==2 with an unobserved second ph
	else if(dimPH>1) { dimPHminusOne = dimPH-1; }
	


	//Total observations
	int Na = sumOfVector(obsA);
	
	if(print == true) {
		cout << "dimFH: " << dimFH << "\n";
		cout << "dimPH: " << dimPH << "\n";
		cout << "Na: " << Na << "\n";
	}


	//Construct matrix T (dim j-1 x k-1)
	gsl_matrix *T = gsl_matrix_alloc(dimPHminusOne,dimFH-1);

	//Find the full haplotype dimension that needs to be omitted fo reduced dimensionality purposes
	int FHdimOmitted = -1;
	if(dimPH == 1) { //Not straightforward in this case

		//WLOG, the first non-negative integer that isn't a contrib of PH0 will be a FH contrib of an unobserved PH1
		for(int i=0; i<9999; i++) {
		
			if(i<(int) (*contribsPHset)[0].size()) {
				if((*contribsPHset)[0][i] != i) { //Assumes contribs ordered. E.g. if contribs of PH0: {0,1,2,4,6} then FHdimOmitted = 3.

					FHdimOmitted = i;
					break;
				}

			} else {
				FHdimOmitted = i; //Situation where e.g. contribs is {0,1,2,3,4}, then FHdimOmitted = 5
				break;
			}
		}
		
	} else { //Simpler, just omit the last PH dimension (j). Then omit the first of the contribs corresponding to this PH dim
		FHdimOmitted = (*contribsPHset)[dimPH-1][0];
	}

	if(print == true) {
		cout << "FHdimOmitted: " << FHdimOmitted << "\n";
	}
	
	for(int i=0; i<dimPHminusOne; i++) { //Loop over rows of T, omit last row

		int counter=0; //Contribs counter
		for(int j=0; j<dimFH; j++) { //Loop over columns of T

			if(j<FHdimOmitted) {
				if((*contribsPHset)[i][counter]==j && counter < (int) (*contribsPHset)[i].size()) { //This assumes the contribs are ordered
					gsl_matrix_set(T,i,j,1);
					counter++;
				} else {
					gsl_matrix_set(T,i,j,0);
				}

			} else if(j>FHdimOmitted) {
				if((*contribsPHset)[i][counter]==j && counter < (int) (*contribsPHset)[i].size()) { //This assumes the contribs are ordered
					gsl_matrix_set(T,i,j-1,1);
					counter++;
				} else {
					gsl_matrix_set(T,i,j-1,0);
				}
			}
		}
	}


	if(print == true) {
		cout << "Matrix T: "; printMatrixMathematica(T);
	}

	//Create reduced quantities
	vector<double> qBmeanReduced;
	for(unsigned int i=0; i<qBmean.size(); i++) {
		if(((int) i)!=FHdimOmitted) {
			qBmeanReduced.push_back(qBmean[i]);
		}
	} 
	vector<double> obsAreduced;
	for(int i=0; i<dimPHminusOne; i++) {
		obsAreduced.push_back((double) obsA[i]);
	}
	//Create reduced covariance matrix
	gsl_matrix *qBvarReduced = gsl_matrix_alloc(dimFH-1,dimFH-1);
	for(int i=0; i<dimFH; i++) {
		for(int j=0; j<dimFH; j++) {
	
			int deltai = 0;
			int deltaj = 0;
			if(i>FHdimOmitted) { deltai = -1; }
			if(j>FHdimOmitted) { deltaj = -1; }

			if(i!= FHdimOmitted && j!=FHdimOmitted) {
			
				gsl_matrix_set(qBvarReduced,i+deltai,j+deltaj, gsl_matrix_get(qBvar,i,j));
			}	
		}
	}



	/*
	* * 
	* * After computation
	* *
	*/


	/*
 	* First we create the mean of xAPH = Na*T*sigmaT[qBmean]
 	*/ 
	vector<double> sigmaTqBmean;
	double totalT = 0;
	for(int i=0; i<dimFH; i++) {
		totalT += qBmean[i]*exp(hapFitT[i]);
	}
	if(print == true) {
		cout << "totalT: " << totalT << "\n";
	}
	for(int i=0; i<dimFH; i++) {
		if(i!=FHdimOmitted) {
			sigmaTqBmean.push_back(qBmean[i]*exp(hapFitT[i])/totalT);
		}
	}


	if(print == true) {
		cout << "sigmaTqBmean: "; printDoubleVector(sigmaTqBmean);
	}

	gsl_vector *sigmaTqBmeanGSL = gsl_vector_alloc(dimFH-1);
	for(unsigned int i=0;i<sigmaTqBmean.size();i++) {
		gsl_vector_set(sigmaTqBmeanGSL,i,sigmaTqBmean[i]);
	}

	//Create TsigmaTqBmean through matrix vector multiplication
	gsl_vector *TsigmaTqBmeanGSL = gsl_vector_alloc(dimPHminusOne);
	gsl_blas_dgemv(CblasNoTrans, 1.0, T, sigmaTqBmeanGSL, 0.0, TsigmaTqBmeanGSL);  

	//Convert back to C++ vector
	vector<double> meanXaPH;
	vector<double> TsigmaTqBmean; //To be used in computing variance
	for(int i=0; i<dimPHminusOne; i++) {
		meanXaPH.push_back(Na*(gsl_vector_get(TsigmaTqBmeanGSL, i)));
		TsigmaTqBmean.push_back(gsl_vector_get(TsigmaTqBmeanGSL, i));
	}

	if(print == true) {
		cout << "meanXaPH: "; printDoubleVector(meanXaPH);
	}


	
	/*
	* Next we compute the variance of xAPH = Na(alpha+(Na-alpha)beta)M(TsigmaTqBmean) + Na(Na-alpha)gammaTDsigmaTvarqBDsigmaTTrans
	*/
	//Constructing the first term Na(alpha+(Na-alpha)beta)*M(TsigmaT[qBmean])
	gsl_matrix *NaAlphaNaMinusAlphaBetaMTsigmaTqBmean = gsl_matrix_alloc(dimPHminusOne,dimPHminusOne);
	constructMatrixM(TsigmaTqBmean, NaAlphaNaMinusAlphaBetaMTsigmaTqBmean);
	
	//Multiplying by the factor Na(alpha+(Na-alpha)beta):
	double alpha = (Na+C)/((double)(1+C));
	int growthFactor = ap->getGrowthFactor();
	int Ng=growthFactor*Nt;
	double NtNgDouble = ((double)Nt)*((double)Ng); //This computation gets out of integer range if Nt>4634, so turn into double
	double beta = (Nt+Ng-1)/NtNgDouble;

	gsl_matrix_scale(NaAlphaNaMinusAlphaBetaMTsigmaTqBmean, Na*(alpha+(Na-alpha)*beta));

	if(print==true) {
		cout << "Na(alpha-(Na-alpha)beta)MTsigmaTqBmean: "; printMatrixMathematica(NaAlphaNaMinusAlphaBetaMTsigmaTqBmean);
	}
	


	//Constructing the second term Na*(Na-alpha)*gamma*T*D[sigmaT]_qBmean * qBvar * (D[sigmaT]_qBmean)^T)T^T
	
	double gamma = (NtNgDouble-Nt-Ng+1)/NtNgDouble;
	double NaNaMinusAlphaGamma = Na*(Na-alpha)*gamma;

	//Computing D[sigmaT]_qBmean
	gsl_matrix* DsigmaT_qBmean = gsl_matrix_alloc(dimFH-1,dimFH-1);
	for(int i=0; i<dimFH; i++) {
		for(int j=0; j<dimFH; j++) {
	
			int deltai=0;
			int deltaj=0;
			if(i>FHdimOmitted) { deltai = -1; }
			if(j>FHdimOmitted) { deltaj = -1; }

			if(i!=FHdimOmitted && j!=FHdimOmitted) {
				if(i==j) {
					gsl_matrix_set(DsigmaT_qBmean, i+deltai, j+deltaj, exp(hapFitT[i])/totalT -qBmean[i]*exp(hapFitT[i]) * (exp(hapFitT[i]) -exp(hapFitT[FHdimOmitted]))/(totalT*totalT));
				} else {
					gsl_matrix_set(DsigmaT_qBmean, i+deltai, j+deltaj, -qBmean[i]*exp(hapFitT[i]) * (exp(hapFitT[j]) - exp(hapFitT[FHdimOmitted]))/(totalT*totalT));
				}
			}
		}
	}

	

	//Compute Na*(Na-alpha)*gamma * T*D[sigmaT]_qBmean * qBvar * (D[sigmaT]_qBmean)^T*T^T
	gsl_matrix* TDsigmaT_qBmean = gsl_matrix_alloc(dimPHminusOne,dimFH-1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, DsigmaT_qBmean, 0.0, TDsigmaT_qBmean);
	gsl_matrix* qBvarDsigmaT_qBmean_Trans_T_Trans = gsl_matrix_alloc(dimFH-1,dimPHminusOne);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, qBvarReduced, TDsigmaT_qBmean, 0.0, qBvarDsigmaT_qBmean_Trans_T_Trans);
	gsl_matrix* NaNaMinusAlphaGammaTDsigmaT_qBmean_qBvarDsigmaT_qBmean_Trans_T_Trans = gsl_matrix_alloc(dimPHminusOne, dimPHminusOne);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, NaNaMinusAlphaGamma, TDsigmaT_qBmean, qBvarDsigmaT_qBmean_Trans_T_Trans, 0.0, NaNaMinusAlphaGammaTDsigmaT_qBmean_qBvarDsigmaT_qBmean_Trans_T_Trans);


	if(print==true) {
		cout << "Na*(Na-alpha)*gamma*T*DsigmaT(qBmean)* qBvar * (D[sigmaT]_qBmean)^T) *  T^T:  "; printMatrixMathematica(NaNaMinusAlphaGammaTDsigmaT_qBmean_qBvarDsigmaT_qBmean_Trans_T_Trans);
	}
	
	//If qBvar is used, add together
	if(ap->getNoVar() == false && ap->getMeanOnly() == false) {

		//Add together
		gsl_matrix_add(NaAlphaNaMinusAlphaBetaMTsigmaTqBmean, NaNaMinusAlphaGammaTDsigmaT_qBmean_qBvarDsigmaT_qBmean_Trans_T_Trans); //Final variance for after case. Output stored in first matrix.

	} else {

		//Don't add together, i.e. NaAlphaNaMinusAlphaBetaMTsigmaTqBmean remains unchanged
	}

	if(print == true) {
		cout << "Final quantities: \n";
		cout << "ObsAreduced: "; printDoubleVector(obsAreduced);
		cout << "meanXaPHreduced: "; printDoubleVector(meanXaPH);
		cout << "varXaPHreduced: "; printMatrixMathematica(NaAlphaNaMinusAlphaBetaMTsigmaTqBmean);
		cout << "Na: " << Na << "\n";
	}
	
	//Computing likelihoods depend on dimensionality

	if(dimPH > 2) { //Multivariate approach

		//Compute likelihood
		logPxAph = logMultivariateNormalPDF(obsAreduced,meanXaPH,NaAlphaNaMinusAlphaBetaMTsigmaTqBmean);
		if(print==true) {
			logMultivariateNormalPDFPrint(obsAreduced,meanXaPH,NaAlphaNaMinusAlphaBetaMTsigmaTqBmean);
		}
		
		if(ap->getUseIntegralApproach() == true) {
		
			double* xmin = (double *) malloc(dimPHminusOne * sizeof(double));
			double* xmax = (double *) malloc(dimPHminusOne * sizeof(double));
			double* err = (double *) malloc(intOutputDim * sizeof(double));
			err[0] = 9999;
			for(int i=0; i<dimPHminusOne; i++) {
		
				xmin[i] = ((double) obsAreduced[i]) - 0.5; //Define hypercube around observation
				xmax[i] = ((double) obsAreduced[i]) + 0.5; //Boundaries may be < 0 or > Na
			}

			parameters params;
			params.muVec = meanXaPH;
			params.Sigma = NaAlphaNaMinusAlphaBetaMTsigmaTqBmean;
			double PxAph = 0; //Not log parameter
		
			//Perform integral. For dim <= 3, pcubature is best. hcubature is better for larger dimensions, although best for <7.
			int intOutcome = -1;
			if(dimPHminusOne <= 3) {
				intOutcome = pcubature(intOutputDim, integrandMultiNormal, &params, dimPHminusOne, xmin, xmax, maxEval, tol, tol, ERROR_INDIVIDUAL, &PxAph, err); 
			} else {
				intOutcome = hcubature(intOutputDim, integrandMultiNormal, &params, dimPHminusOne, xmin, xmax, maxEval, tol, tol, ERROR_INDIVIDUAL, &PxAph, err); 
			}

			bool error = false;
			if(intOutcome != 0) { //Error in integral

				cout << "Error computing numerical integral. Using original (compound) result instead.\n";
				error = true;
			}
			if(err[0] > 0.001) {
				
				cout << "Error in numerical integral is too large. The error is: " << err[0] << ". Using original (compound) result instead.\n";
				error = true;
			}
			double logPxAphInt = log(PxAph);
			if(::isinf(logPxAphInt) != 0) { //integral evaluation is infinity

				cout << "Error when computing numerical integral: log likelihood is infinity. Using original (compound) result instead.\n";
				error = true;
			}

			if(::isnan(logPxAphInt)!=0) { //integral evaluation is nan

				cout << "Error when computing numerical integral: log likelihood is nan. Using original (compound) result instead.\n";
				error = true;
			}
			if(error == false) { //Using integral approach value if no errors were encountered
				logPxAph = logPxAphInt;
				cout << "Using integral result.\n";
			}
			cout << "logPxAph: " << logPxAph << "\n";

		}
		
        
	} else { //I.e. dimPH == 2 or dimPH == 1 (same computation)
	
		//In this case we retrieve the univariate case
		double meanXaPHuni = meanXaPH[0];
		double varXaPHuni = gsl_matrix_get(NaAlphaNaMinusAlphaBetaMTsigmaTqBmean,0,0);
		double stDevXaPHuni = sqrt(varXaPHuni);
		logPxAph = -(obsA[0]-meanXaPHuni)*(obsA[0]-meanXaPHuni)/(2*varXaPHuni) - log(stDevXaPHuni*sqrt(2*M_PI)); //log L
		

	
		if(ap->getUseIntegralApproach() == true) {

			double* xmin = (double *) malloc(1 * sizeof(double));
			double* xmax = (double *) malloc(1 * sizeof(double));
			xmin[0] = ((double) obsAreduced[0]) - 0.5; //Define +/- 0.5 interval around observation
			xmax[0] = ((double) obsAreduced[0]) + 0.5; //Boundaries may be < 0 or > Na
			double* err = (double *) malloc(intOutputDim * sizeof(double));
			err[0] = 9999;

			parameters params;
			params.mu = meanXaPHuni;
			params.var = varXaPHuni;
			params.stDev = stDevXaPHuni;
			double PxAph = 0; //Not log parameter
		
			//Perform integral using pcubature (best for dim <=3)
			int intOutcome = pcubature(intOutputDim, integrandNormal, &params, 1, xmin, xmax, maxEval, tol, tol, ERROR_INDIVIDUAL, &PxAph, err); 
			bool error = false;
			if(intOutcome != 0) { //Error in integral

				cout << "Error computing numerical integral. Using original (compound) result instead.\n";
				error = true;
			}
			if(err[0] > 0.001) {
				
				cout << "Error in numerical integral is too large. The error is: " << err[0] << ". Using original (compound) result instead.\n";
				error = true;
			}
			double logPxAphInt = log(PxAph);
			if(::isinf(logPxAphInt) != 0) { //integral evaluation is infinity

				cout << "Error when computing numerical integral: log likelihood is infinity. Using original (compound) result instead.\n";
				error = true;
			}

			if(::isnan(logPxAphInt)!=0) { //integral evaluation is nan

				cout << "Error when computing numerical integral: log likelihood is nan. Using original (compound) result instead.\n";
				error = true;
			}
			if(error == false) { //Using integral approach value if no errors were encountered
				logPxAph = logPxAphInt;
				cout << "Using integral result.\n";
			}
			cout << "logPxAph: " << logPxAph << "\n";
		}
	}


	//Check of likelihoods are still negative infinity, which should not happen
	if(logPxAph == -numeric_limits<double>::max() && print == false) {
		cout << "ERROR in Misc::computeLSelTVar: Likelihood for xA is negative infinity. LogPxAph: " << logPxAph <<"\n";

		cout << "Recomputing mean and variance whilst printing out stats:\n";
		computeLSelTVar(obsA,Nt,ap,qBmean, qBvar, hapFitT, contribsPHset, true);

		cout << "-------------------------------------------------------\n";
	}

	//Likelihood is nan, so print out
	if(::isnan(logPxAph) != 0 && print == false) {
		cout << "ERROR in Misc::computeLSelTVar: Likelihood for xA is nan. LogPxAph: " << logPxAph <<"\n";

		cout << "Recomputing mean and variance whilst printing out stats:\n";
		computeLSelTVar(obsA,Nt,ap,qBmean, qBvar, hapFitT, contribsPHset, true);

		cout << "-------------------------------------------------------\n";
	}

	//Clean up before exit
	gsl_vector_free(sigmaTqBmeanGSL);
	gsl_vector_free(TsigmaTqBmeanGSL);

	gsl_matrix_free(T);
	gsl_matrix_free(NaAlphaNaMinusAlphaBetaMTsigmaTqBmean);
	gsl_matrix_free(DsigmaT_qBmean);
	gsl_matrix_free(TDsigmaT_qBmean);
	gsl_matrix_free(qBvarDsigmaT_qBmean_Trans_T_Trans);
	gsl_matrix_free(NaNaMinusAlphaGammaTDsigmaT_qBmean_qBvarDsigmaT_qBmean_Trans_T_Trans);
	gsl_matrix_free(qBvarReduced);

	if(print == true) {

		cout << "Returning likelihood " << logPxAph << "\n";
	}

	return logPxAph;
}


double computeLSelTVarOld(vector<int> &obsA, int Nt, AnaParam *ap, vector<double> &qBmean, gsl_matrix* qBvar, vector<double>& hapFitT, vector<vector<int> >* contribsPHset, bool print) {

	double C = ap->getC();

	if(print == true) {
	
		cout << "Printing intermediate steps in computeLselection.\n";
		cout << "First, the input variables:\n";
		cout << "obsA: "; printIntVector(obsA);
		cout << "Nt: " << Nt << "\n";
		cout << "C: " << C << "\n";
		cout << "qBmean: "; printDoubleVector(qBmean);
		cout << "Contribs: \n";
		for(unsigned int i=0; i<contribsPHset->size(); i++) {
			cout << "PH " << i << ": "; printIntVector((*contribsPHset)[i]);
		}


	}

	//Dimensions
	int dimFH = (int) qBmean.size();
	int dimPH = (int) obsA.size();
	

	//Total observations
	int Na = sumOfVector(obsA);
	
	if(print == true) {
		cout << "dimFH: " << dimFH << "\n";
		cout << "dimPH: " << dimPH << "\n";
		cout << "Na: " << Na << "\n";
	}

	//Construct matrix T
	gsl_matrix *T = gsl_matrix_alloc(dimPH,dimFH);
	for(int i=0; i<dimPH; i++) { //Loop over rows

		int counter=0;
		for(int j=0; j<dimFH; j++) { //Loop over columns

			if((*contribsPHset)[i][counter]==j && counter < (int) (*contribsPHset)[i].size()) { //This assumes the contribs are ordered
				gsl_matrix_set(T,i,j,1);
				counter++;
			} else {
				gsl_matrix_set(T,i,j,0);
			}
		}
	}

	if(print == true) {
		cout << "Matrix T: "; printMatrixMathematica(T);
	}


	/*
	* * 
	* * After computation
	* *
	*/


	/*
 	* First we create the mean of xAPH = Na*T*sigmaT[qBmean]
 	*/ 
	vector<double> sigmaTqBmean;
	double totalT = 0;
	for(int i=0; i<dimFH; i++) {
		totalT += qBmean[i]*exp(hapFitT[i]);
	}
	for(int i=0; i<dimFH; i++) {
		sigmaTqBmean.push_back(qBmean[i]*exp(hapFitT[i])/totalT);
	}


	if(print == true) {
		cout << "sigmaTqBmean: "; printDoubleVector(sigmaTqBmean);
	}

	gsl_vector *sigmaTqBmeanGSL = gsl_vector_alloc(dimFH);
	for(int i=0;i<dimFH;i++) {
		gsl_vector_set(sigmaTqBmeanGSL,i,sigmaTqBmean[i]);
	}

	//Create TsigmaTqBmean through matrix vector multiplication
	gsl_vector *TsigmaTqBmeanGSL = gsl_vector_alloc(dimPH);
	gsl_blas_dgemv(CblasNoTrans, 1.0, T, sigmaTqBmeanGSL, 0.0, TsigmaTqBmeanGSL);  

	//Convert back to C++ vector
	vector<double> meanXaPH;
	vector<double> TsigmaTqBmean; //To be used in computing variance
	for(int i=0; i<dimPH; i++) {
		meanXaPH.push_back(Na*(gsl_vector_get(TsigmaTqBmeanGSL, i)));
		TsigmaTqBmean.push_back(gsl_vector_get(TsigmaTqBmeanGSL, i));
	}

	if(print == true) {
		cout << "meanXaPH: "; printDoubleVector(meanXaPH);
	}


	
	/*
	* Next we compute the variance of xAPH = Na(alpha+(Na-alpha)beta)M(TsigmaTqBmean) + Na(Na-alpha)gammaTDsigmaTvarqBDsigmaTTrans
	*/
	//Constructing the first term Na(alpha+(Na-alpha)beta)*M(TsigmaT[qBmean])
	gsl_matrix *NaAlphaNaMinusAlphaBetaMTsigmaTqBmean = gsl_matrix_alloc(dimPH,dimPH);
	constructMatrixM(TsigmaTqBmean, NaAlphaNaMinusAlphaBetaMTsigmaTqBmean);
	
	//Multiplying by the factor Na(alpha+(Na-alpha)beta):
	double alpha = (Na+C)/((double)(1+C));
	int growthFactor = ap->getGrowthFactor();
	int Ng=growthFactor*Nt;
	double NtNgDouble = ((double)Nt)*((double)Ng); //This computation gets out of integer range if Nt>4634, so turn into double
	double beta = (Nt+Ng-1)/NtNgDouble;

	gsl_matrix_scale(NaAlphaNaMinusAlphaBetaMTsigmaTqBmean, Na*(alpha+(Na-alpha)*beta));

	if(print==true) {
		cout << "Na(alpha-(Na-alpha)beta)MTsigmaTqBmean: "; printMatrixMathematica(NaAlphaNaMinusAlphaBetaMTsigmaTqBmean);
	}
	


	//Constructing the second term Na*(Na-alpha)*gamma*T*D[sigmaT]_qBmean * qBvar * (D[sigmaT]_qBmean)^T)T^T
	
	double gamma = (NtNgDouble-Nt-Ng+1)/NtNgDouble;
	double NaNaMinusAlphaGamma = Na*(Na-alpha)*gamma;


	//Compute D[sigmaT]_qBmean
	gsl_matrix* DsigmaT_qBmean = gsl_matrix_alloc(dimFH,dimFH);
	for(int i=0; i<dimFH; i++) {
		for(int j=0; j<dimFH; j++) {
	
			if(i==j) {
				gsl_matrix_set(DsigmaT_qBmean, i, j, exp(hapFitT[i])/totalT - qBmean[i]*exp(2*hapFitT[i])/(totalT*totalT));
			} else {
				gsl_matrix_set(DsigmaT_qBmean, i, j, -qBmean[i]*exp(hapFitT[i]+hapFitT[j])/(totalT*totalT));
			}
		}
	}	

	//Compute Na*(Na-alpha)*gamma * T*D[sigmaT]_qBmean * qBvar * (D[sigmaT]_qBmean)^T*T^T
	gsl_matrix* TDsigmaT_qBmean = gsl_matrix_alloc(dimPH,dimFH);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, DsigmaT_qBmean, 0.0, TDsigmaT_qBmean);
	gsl_matrix* qBvarDsigmaT_qBmean_Trans_T_Trans = gsl_matrix_alloc(dimFH,dimPH);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, qBvar, TDsigmaT_qBmean, 0.0, qBvarDsigmaT_qBmean_Trans_T_Trans);
	gsl_matrix* NaNaMinusAlphaGammaTDsigmaT_qBmean_qBvarDsigmaT_qBmean_Trans_T_Trans = gsl_matrix_alloc(dimPH, dimPH);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, NaNaMinusAlphaGamma, TDsigmaT_qBmean, qBvarDsigmaT_qBmean_Trans_T_Trans, 0.0, NaNaMinusAlphaGammaTDsigmaT_qBmean_qBvarDsigmaT_qBmean_Trans_T_Trans);


	if(print==true) {
		cout << "Na*(Na-alpha)*gamma*T*DsigmaT(qBmean)* qBvar * (D[sigmaT]_qBmean)^T) *  T^T:  "; printMatrixMathematica(NaNaMinusAlphaGammaTDsigmaT_qBmean_qBvarDsigmaT_qBmean_Trans_T_Trans);
	}
	
	//Add together
	gsl_matrix_add(NaAlphaNaMinusAlphaBetaMTsigmaTqBmean, NaNaMinusAlphaGammaTDsigmaT_qBmean_qBvarDsigmaT_qBmean_Trans_T_Trans); //Final variance for after case. Output stored in first matrix.

	
	//Computing likelihoods depend on dimensionality
	double logPxAph = -numeric_limits<double>::max(); //Set to neg infinity as precaution

	if(dimPH > 1) { //Most common case, i.e. multiple partial haps in a set

		/*
		 * Reduce dimensionality by 1 to ensure non-degeneracy
		 */
		//After case
		meanXaPH.pop_back(); //WLOG remove the last PH
		gsl_matrix_view varXaPHreducedView = gsl_matrix_submatrix(NaAlphaNaMinusAlphaBetaMTsigmaTqBmean, 0,0, dimPH-1, dimPH-1); //Remove row k and column k
		gsl_matrix * varXaPH = gsl_matrix_alloc(dimPH-1, dimPH-1);
		gsl_matrix_memcpy(varXaPH, &(varXaPHreducedView.matrix)); //Copy the reduced matrix view into a "normal" matrix varXaPH
		vector<double> obsAreduced;
		for(int i=0; i<dimPH-1; i++) {
			obsAreduced.push_back(obsA[i]);
		}

		if(print == true) {
			cout << "Reduced quantities: \n";
			cout << "ObsAreduced: "; printDoubleVector(obsAreduced);
			cout << "meanXaPHreduced: "; printDoubleVector(meanXaPH);
			cout << "varXaPHreduced: "; printMatrixMathematica(varXaPH);
		}

		if(dimPH > 2) { //Multivariate approach	

			//Compute likelihood
			logPxAph = logMultivariateNormalPDF(obsAreduced,meanXaPH,varXaPH);
		
        
		}else { //I.e. dimPH==2
	
			//In this case we retrieve the univariate case when reducing dimensionality by one
			double meanXaPHuni = meanXaPH[0];
			double varXaPHuni = gsl_matrix_get(varXaPH,0,0);
			double stDevXaPHuni = sqrt(varXaPHuni);

			logPxAph = -(obsA[0]-meanXaPHuni)*(obsA[0]-meanXaPHuni)/(2*varXaPHuni) - log(stDevXaPHuni*sqrt(2*M_PI)); //log L
		}

		//Deallocate variables created in this scope
		gsl_matrix_free(varXaPH);

	//In theory, can combine dimPH==2 from above with this one, would be slightly faster, although more clear this way
	} else if (dimPH==1) {
		//Case where there is a single partial haplotype in a set.
		//This can be thought of as a case with a single observed partial haplotype
		//and one (or several) unobserved haplotypes.
		//If all loci in the partial haplotype set are variants (i.e. not monomorphic)
		//then there MUST exists at least one other haplotype containing the unobserved alleles,
		//however, this (or these) haplotype(s) are not observed.
		//As such, WLOG we can consider this a case of dimPH=2, which gets reduced to a one dimensional system
		//under reduction, i.e. similar to the dimPH==2 case of above.
				
		double meanXaPHuni = meanXaPH[0];
		double varXaPHuni = gsl_matrix_get(NaAlphaNaMinusAlphaBetaMTsigmaTqBmean,0,0);
		double stDevXaPHuni = sqrt(varXaPHuni);	

		logPxAph = -(obsA[0]-meanXaPHuni)*(obsA[0]-meanXaPHuni)/(2*varXaPHuni) - log(stDevXaPHuni*sqrt(2*M_PI)); //log L

	} else { //error
	
		cout << "ERROR in computeLSelT: dimPH is <1\n";
	}

	//Check of likelihoods are still negative infinity, which should not happen
	if(logPxAph == -numeric_limits<double>::max() && print == false) {
		cout << "ERROR in Misc::computeLSelTVar: Likelihood for xA is negative infinity. LogPxAph: " << logPxAph <<"\n";

		cout << "Recomputing mean and variance whilst printing out stats:\n";
		computeLSelTVar(obsA,Nt,ap,qBmean, qBvar, hapFitT, contribsPHset, true);

		cout << "-------------------------------------------------------\n";
	}


	//Clean up before exit
	gsl_vector_free(sigmaTqBmeanGSL);
	gsl_vector_free(TsigmaTqBmeanGSL);

	gsl_matrix_free(T);
	gsl_matrix_free(NaAlphaNaMinusAlphaBetaMTsigmaTqBmean);
	gsl_matrix_free(DsigmaT_qBmean);
	gsl_matrix_free(TDsigmaT_qBmean);
	gsl_matrix_free(qBvarDsigmaT_qBmean_Trans_T_Trans);
	gsl_matrix_free(NaNaMinusAlphaGammaTDsigmaT_qBmean_qBvarDsigmaT_qBmean_Trans_T_Trans);


	return logPxAph;
}




double computeLSelTSelGVar(vector<int> &obsA, int Nt, AnaParam *ap, vector<double> &qBmean, gsl_matrix* qBvar, vector<double>& hapFitT, vector<double>& hapFitG, vector<vector<int> >* contribsPHset, bool print) {

	double C = ap->getC();
	double logPxAph = -numeric_limits<double>::max(); //Set to neg infinity as precaution

	if(print == true) {
	
		cout << "Printing intermediate steps in computeLselection.\n";
		cout << "First, the input variables:\n";
		cout << "obsA: "; printIntVector(obsA);
		cout << "Nt: " << Nt << "\n";
		cout << "C: " << C << "\n";
		cout << "qBmean: "; printDoubleVector(qBmean);
		cout << "qBvar: "; printMatrixMathematica(qBvar);
		cout << "hapFitT: "; printDoubleVector(hapFitT);
		cout << "hapFitG: "; printDoubleVector(hapFitG);
		vector<double> hapFitTPlusHapFitG = addVectorsDouble(hapFitT, hapFitG);
		cout << "hapFitT + hapFitG = "; printDoubleVector(hapFitTPlusHapFitG);
		cout << "Contribs: \n";
		for(unsigned int i=0; i<contribsPHset->size(); i++) {
			cout << "PH " << i << ": "; printIntVector((*contribsPHset)[i]);
		}


	}

	//Dimensions
	int dimFH = (int) qBmean.size();
	int dimPH = (int) obsA.size();
	int dimPHminusOne = -1;
	if(dimPH==1) { dimPHminusOne = 1; } //If dimPH==1, this corresponds to at least dimPH==2 with an unobserved second ph
	else if(dimPH>1) { dimPHminusOne = dimPH-1; }


	//Total observations
	int Na = sumOfVector(obsA);
	
	if(print == true) {
		cout << "dimFH: " << dimFH << "\n";
		cout << "dimPH: " << dimPH << "\n";
		cout << "Na: " << Na << "\n";
	}

	//Construct matrix T (dim j-1 x k-1)
	gsl_matrix *T = gsl_matrix_alloc(dimPHminusOne,dimFH-1);

	//Find the full haplotype dimension that needs to be omitted fo reduced dimensionality purposes
	int FHdimOmitted = -1;
	if(dimPH == 1) { //Not straightforward in this case

		//WLOG, the first non-negative integer that isn't a contrib of PH0 will be a FH contrib of an unobserved PH1
		//cout << "contribsPHset[0]: "; printIntVector((*contribsPHset)[0]);
		for(int i=0; i<9999; i++) {
		
			if(i<(int)(*contribsPHset)[0].size()) {
				if((*contribsPHset)[0][i] != i) { //Assumes contribs ordered. E.g. if contribs of PH0: {0,1,2,4,6} then FHdimOmitted = 3.

					FHdimOmitted = i;
					break;
				}

			} else {
				FHdimOmitted = i; //Situation where e.g. contribs is {0,1,2,3,4}, then FHdimOmitted = 5
				break;
			}
		}
		
	} else { //Simpler, just omit the last PH dimension (j). Then omit the first of the contribs corresponding to this PH dim
		FHdimOmitted = (*contribsPHset)[dimPH-1][0];
	}

	if(print == true) {
		cout << "FHdimOmitted: " << FHdimOmitted << "\n";
	}
	
	for(int i=0; i<dimPHminusOne; i++) { //Loop over rows of T, omit last row

		int counter=0; //Contribs counter
		for(int j=0; j<dimFH; j++) { //Loop over columns of T

			if(j<FHdimOmitted) {
				if((*contribsPHset)[i][counter]==j && counter < (int) (*contribsPHset)[i].size()) { //This assumes the contribs are ordered
					gsl_matrix_set(T,i,j,1);
					counter++;
				} else {
					gsl_matrix_set(T,i,j,0);
				}

			} else if(j>FHdimOmitted) {
				if((*contribsPHset)[i][counter]==j && counter < (int) (*contribsPHset)[i].size()) { //This assumes the contribs are ordered
					gsl_matrix_set(T,i,j-1,1);
					counter++;
				} else {
					gsl_matrix_set(T,i,j-1,0);
				}
			}
		}
	}


	if(print == true) {
		cout << "Matrix T: "; printMatrixMathematica(T);
	}

	//Create reduced quantities
	//TODO: Probably faster to copy qBmean and then erasing specific index rather than pushing back a large number of entries
	vector<double> qBmeanReduced;
	for(unsigned int i=0; i<qBmean.size(); i++) {
		if(((int) i)!=FHdimOmitted) {
			qBmeanReduced.push_back(qBmean[i]);
		}
	} 
	vector<double> obsAreduced;
	for(int i=0; i<dimPHminusOne; i++) {
		obsAreduced.push_back((double) obsA[i]);
	}
	//Create reduced covariance matrix
	gsl_matrix *qBvarReduced = gsl_matrix_alloc(dimFH-1,dimFH-1);
	for(int i=0; i<dimFH; i++) {
		for(int j=0; j<dimFH; j++) {
	
			int deltai = 0;
			int deltaj = 0;
			if(i>FHdimOmitted) { deltai = -1; }
			if(j>FHdimOmitted) { deltaj = -1; }

			if(i!= FHdimOmitted && j!=FHdimOmitted) {
			
				gsl_matrix_set(qBvarReduced,i+deltai,j+deltaj, gsl_matrix_get(qBvar,i,j));
			}	
		}
	}




	/*
	* * 
	* * After computation
	* *
	*/


	/*
 	* First we create the mean of xAPH = Na*T*sigmaG[sigmaT[qBmean]]
 	*/ 
	vector<double> sigmaTqBmean;
	vector<double> sigmaTqBmeanFullDim;
	double totalT = 0;
	for(int i=0; i<dimFH; i++) {
		totalT += qBmean[i]*exp(hapFitT[i]);
	}
	if(print==true) {
		cout << "totalT: " << totalT << "\n";
	}
	for(int i=0; i<dimFH; i++) {
		if(i!=FHdimOmitted) {
			sigmaTqBmean.push_back(qBmean[i]*exp(hapFitT[i])/totalT);
		}
		sigmaTqBmeanFullDim.push_back(qBmean[i]*exp(hapFitT[i])/totalT);
	}

	vector<double> sigmaGsigmaTqBmean;
	double totalG = 0;
	for(int i=0; i<dimFH; i++) {
		totalG += sigmaTqBmeanFullDim[i]*exp(hapFitG[i]);
	}
	if(print==true) {
		cout << "totalG: " << totalG << "\n";
	}
	for(int i=0; i<dimFH; i++) {
		if(i!=FHdimOmitted) {
			sigmaGsigmaTqBmean.push_back(sigmaTqBmeanFullDim[i]*exp(hapFitG[i])/totalG);
		}
	}

	//Rescale sigmaGsigmaTqBmean to avoid issues with frequencies < 10^-10 arising from very large within-host selection coefficients
	rescale(sigmaGsigmaTqBmean);

	if(print == true) {
		cout << "qBmean: "; printDoubleVector(qBmean);
		cout << "sigmaTqBmean: "; printDoubleVector(sigmaTqBmean);
		cout << "sigmaGsigmaTqBmean: "; printDoubleVector(sigmaGsigmaTqBmean);
	}

	gsl_vector *sigmaGsigmaTqBmeanGSL = gsl_vector_alloc(dimFH-1);
	for(unsigned int i=0;i<sigmaGsigmaTqBmean.size();i++) {
		gsl_vector_set(sigmaGsigmaTqBmeanGSL,i,sigmaGsigmaTqBmean[i]);
	}

	//Create TsigmaGsigmaTqBmean through matrix vector multiplication
	gsl_vector *TsigmaGsigmaTqBmeanGSL = gsl_vector_alloc(dimPHminusOne);
	gsl_blas_dgemv(CblasNoTrans, 1.0, T, sigmaGsigmaTqBmeanGSL, 0.0, TsigmaGsigmaTqBmeanGSL);  

	//Convert back to C++ vector
	vector<double> meanXaPH;
	vector<double> TsigmaGsigmaTqBmean; //To be used in computing variance
	for(int i=0; i<dimPHminusOne; i++) {
		meanXaPH.push_back(Na*(gsl_vector_get(TsigmaGsigmaTqBmeanGSL, i)));
		TsigmaGsigmaTqBmean.push_back(gsl_vector_get(TsigmaGsigmaTqBmeanGSL, i));
	}

	if(print == true) {
		cout << "meanXaPH: "; printDoubleVector(meanXaPH);
	}


	
	/*
	* Next we compute the variance of xAPH = alpha*Na*M(TsigmaG[sigmaT[qBmean]]) + Na*(Na-alpha)*T*D[sigmaG]_sigmaT(qBmean) * (beta*M(sigmaT[qBmean]) + gamma*D[sigmaT]_qBmean * qBvar * (D[sigmaT]_qBmean)^T) * (D[sigmaG]_sigmaT[qBmean])^T * T^T
	*/
	//Constructing the first term alpha*Na*M(TsigmaG[sigmaT[qBmean]])
	gsl_matrix *alphaNaMTsigmaGsigmaTqBmean = gsl_matrix_alloc(dimPHminusOne,dimPHminusOne);
	constructMatrixM(TsigmaGsigmaTqBmean, alphaNaMTsigmaGsigmaTqBmean);
	
	//Multiplying by the factor alpha*Na:
	double alpha = (Na+C)/((double)(1+C));
	gsl_matrix_scale(alphaNaMTsigmaGsigmaTqBmean, alpha*Na);

	if(print==true) {
		cout << "alpha*Na*MTsigmaGsigmaTqBmean: "; printMatrixMathematica(alphaNaMTsigmaGsigmaTqBmean);
	}
	


	//For constructing the second term Na*(Na-alpha)*T*D[sigmaG]_sigmaT(qBmean) * (beta*M(sigmaT[qBmean]) + gamma*D[sigmaT]_qBmean * qBvar * (D[sigmaT]_qBmean)^T) * (D[sigmaG]_sigmaT[qBmean])^T * T^T
	//we break it up into 4 parts:
	//1) Compute beta*M(sigmaT[qBmean]) + gamma*D[sigmaT]_qBmean * qBvar * (D[sigmaT]_qBmean)^T
	//2) Compute T*D[sigmaG]_sigmaT[qBmean]
	//3) Compute 1*2^T
	//4) Compute Na(Na-alpha)*2*(1*2^T)

	//Doing 1) first
	//1a) Compute betaMsigmaTqBmean
	gsl_matrix *betaMsigmaTqBmean = gsl_matrix_alloc(dimFH-1, dimFH-1);
	constructMatrixM(sigmaTqBmean, betaMsigmaTqBmean);
	int growthFactor = ap->getGrowthFactor();
	int Ng=growthFactor*Nt;
	double NtNgDouble = ((double)Nt)*((double)Ng); //This computation gets out of integer range if Nt>4634, so turn into double
	double beta = (Nt+Ng-1)/NtNgDouble;
	double gamma = (NtNgDouble-Nt-Ng+1)/NtNgDouble;
	gsl_matrix_scale(betaMsigmaTqBmean, beta);

	if(print == true) {
		cout << "betaMsigmaTqBmean: "; printMatrixMathematica(betaMsigmaTqBmean);
	}

	//1b) Compute D[sigmaT]_qBmean
	gsl_matrix* DsigmaT_qBmean = gsl_matrix_alloc(dimFH-1,dimFH-1);
	for(int i=0; i<dimFH; i++) {
		for(int j=0; j<dimFH; j++) {
	
			int deltai=0;
			int deltaj=0;
			if(i>FHdimOmitted) { deltai = -1; }
			if(j>FHdimOmitted) { deltaj = -1; }

			if(i!=FHdimOmitted && j!=FHdimOmitted) {
				if(i==j) {
					gsl_matrix_set(DsigmaT_qBmean, i+deltai, j+deltaj, exp(hapFitT[i])/totalT -qBmean[i]*exp(hapFitT[i]) * (exp(hapFitT[i]) -exp(hapFitT[FHdimOmitted]))/(totalT*totalT));
				} else {
					gsl_matrix_set(DsigmaT_qBmean, i+deltai, j+deltaj, -qBmean[i]*exp(hapFitT[i]) * (exp(hapFitT[j]) - exp(hapFitT[FHdimOmitted]))/(totalT*totalT));
				}
			}
		}
	}
	

	//1c) Compute gamma*D[sigmaT]_qBmean * qBvar * (D[sigmaT]_qBmean)^T + betaMsigmaTqBmean
	gsl_matrix* qBvarDsigmaT_qBmean_Trans = gsl_matrix_alloc(dimFH-1,dimFH-1);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, qBvarReduced, DsigmaT_qBmean, 0.0, qBvarDsigmaT_qBmean_Trans);
	if(print == true) {
		cout << "qBvarDsigmaT_qBmean_Trans: "; printMatrixMathematica(qBvarDsigmaT_qBmean_Trans);
	}

	gsl_matrix* gammaDsigmaT_qBmean_qBvarDsigmaT_qBmean_Trans = gsl_matrix_alloc(dimFH-1, dimFH-1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, gamma, DsigmaT_qBmean, qBvarDsigmaT_qBmean_Trans, 0.0, gammaDsigmaT_qBmean_qBvarDsigmaT_qBmean_Trans);
	if(print == true) {
		cout << "gammaDsigmaT_qBmean_qBvarDsigmaT_qBmean_Trans: "; printMatrixMathematica(gammaDsigmaT_qBmean_qBvarDsigmaT_qBmean_Trans);
	}

	//If qBvar is used, add together
	if(ap->getNoVar() == false && ap->getMeanOnly() == false) {

		//1d) Perform 1a+1b (resulting variable is 1a)	
		gsl_matrix_add(betaMsigmaTqBmean, gammaDsigmaT_qBmean_qBvarDsigmaT_qBmean_Trans);
		if(print == true) {
			cout << "1a+1b: "; printMatrixMathematica(betaMsigmaTqBmean);
		}
	} else {

		//Don't add together, i.e. betaMsigmaTqBmean remains unchanged
	}

	//Computing term 2: T*D[sigmaG]_sigmaT[qBmean]
	//2A) Computing D[sigmaG]_sigmaT[qBmean]
	gsl_matrix* DsigmaG_sigmaTqBmean = gsl_matrix_alloc(dimFH-1,dimFH-1);
	for(int i=0; i<dimFH; i++) {
		for(int j=0; j<dimFH; j++) {
	
			int deltai=0;
			int deltaj=0;
			if(i>FHdimOmitted) { deltai = -1; }
			if(j>FHdimOmitted) { deltaj = -1; }

			if(i!=FHdimOmitted && j!=FHdimOmitted) {
				if(i==j) {
					gsl_matrix_set(DsigmaG_sigmaTqBmean, i+deltai, j+deltaj, exp(hapFitG[i])/totalG -sigmaTqBmean[i+deltai]*exp(hapFitG[i]) * (exp(hapFitG[i]) -exp(hapFitG[FHdimOmitted]))/(totalG*totalG));
				} else {
					gsl_matrix_set(DsigmaG_sigmaTqBmean, i+deltai, j+deltaj, -sigmaTqBmean[i+deltai]*exp(hapFitG[i]) * (exp(hapFitG[j]) - exp(hapFitG[FHdimOmitted]))/(totalG*totalG));
				}
			}
		}
	}


	//2B)Compute T*D product
	gsl_matrix* TDsigmaG_sigmaTqBmean = gsl_matrix_alloc(dimPHminusOne,dimFH-1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, DsigmaG_sigmaTqBmean, 0.0, TDsigmaG_sigmaTqBmean);
	if(print == true) {
		cout << "TDsigmaG_sigmaTqBmean: "; printMatrixMathematica(TDsigmaG_sigmaTqBmean);
	}

	//3) Computing term 3 = 1*2^T
	gsl_matrix* three = gsl_matrix_alloc(dimFH-1,dimPHminusOne);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, betaMsigmaTqBmean, TDsigmaG_sigmaTqBmean, 0.0, three);
	if(print == true) {
		cout << "3: "; printMatrixMathematica(three);
	}

	//4) Computing term four = Na(Na-alpha)*2*3
	gsl_matrix* four = gsl_matrix_alloc(dimPHminusOne,dimPHminusOne);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, Na*(Na-alpha), TDsigmaG_sigmaTqBmean, three, 0.0, four);


	if(print==true) {
		cout << "Na*(Na-alpha)*T*D[sigmaG]_sigmaT(qBmean) * (beta*M(sigmaT[qBmean]) + gamma*D[sigmaT]_qBmean * qBvar * (D[sigmaT]_qBmean)^T) * (D[sigmaG]_sigmaT[qBmean])^T * T^T:  "; printMatrixMathematica(four);
	}
	
	//Add together
	gsl_matrix_add(alphaNaMTsigmaGsigmaTqBmean, four); //Final variance for after case. Output stored in first matrix.

	if(print == true) {
		cout << "Final quantities: \n";
		cout << "ObsAreduced: "; printDoubleVector(obsAreduced);
		cout << "meanXaPHreduced: "; printDoubleVector(meanXaPH);
		cout << "varXaPHreduced: "; printMatrixMathematica(alphaNaMTsigmaGsigmaTqBmean);
	}
	
	//Computing likelihoods depend on dimensionality
	if(dimPH > 2) { //Multivariate approach	

		//Compute likelihood
		logPxAph = logMultivariateNormalPDF(obsAreduced,meanXaPH,alphaNaMTsigmaGsigmaTqBmean);
		if(print==true) {
			logMultivariateNormalPDFPrint(obsAreduced,meanXaPH,alphaNaMTsigmaGsigmaTqBmean);
		}
	
		
		if(ap->getUseIntegralApproach() == true) {
			double* xmin = (double *) malloc(dimPHminusOne * sizeof(double));
			double* xmax = (double *) malloc(dimPHminusOne * sizeof(double));
			double* err = (double *) malloc(intOutputDim * sizeof(double));
			err[0] = 9999;	
			for(int i=0; i<dimPHminusOne; i++) {
		
				xmin[i] = ((double) obsAreduced[i]) - 0.5; //Define hypercube around observation
				xmax[i] = ((double) obsAreduced[i]) + 0.5; //Boundaries may be < 0 or > Na
			}

			parameters params;
			params.muVec = meanXaPH;
			params.Sigma = alphaNaMTsigmaGsigmaTqBmean;
			double PxAph = 0; //Not log parameter;
		
			//Perform integral. For dim <= 3, pcubature is best. hcubature is better for larger dimensions, although best for <7.
			int intOutcome = -1;
			if(dimPHminusOne <= 3) {
				intOutcome = pcubature(intOutputDim, integrandMultiNormal, &params, dimPHminusOne, xmin, xmax, maxEval, tol, tol, ERROR_INDIVIDUAL, &PxAph, err); 
			} else {
				intOutcome = hcubature(intOutputDim, integrandMultiNormal, &params, dimPHminusOne, xmin, xmax, maxEval, tol, tol, ERROR_INDIVIDUAL, &PxAph, err); 
			}
			bool error = false;
			if(intOutcome != 0) { //Error in integral

				cout << "Error computing numerical integral. Using original (compound) result instead.\n";
				error = true;
			}
			if(err[0] > 0.001) {
				
				cout << "Error in numerical integral is too large. The error is: " << err[0] << ". Using original (compound) result instead.\n";
				error = true;
			}
			double logPxAphInt = log(PxAph);
			if(::isinf(logPxAphInt) != 0) { //integral evaluation is infinity

				cout << "Error when computing numerical integral: log likelihood is infinity. Using original (compound) result instead.\n";
				error = true;
			}
			if(::isnan(logPxAphInt)!=0) { //integral evaluation is nan

				cout << "Error when computing numerical integral: log likelihood is nan. Using original (compound) result instead.\n";
				error = true;
			}

			if(error == false) { //Using integral approach value if no errors were encountered
				logPxAph = logPxAphInt;
				cout << "Using integral result.\n";
			}
			cout << "logPxAph: " << logPxAph << "\n";

		}
		
        
	}else { //I.e. dimPH == 2 or dimPH == 1 (same computation)
	
		//In this case we retrieve the univariate case
		double meanXaPHuni = meanXaPH[0];
		double varXaPHuni = gsl_matrix_get(alphaNaMTsigmaGsigmaTqBmean,0,0);
		double stDevXaPHuni = sqrt(varXaPHuni);
		logPxAph = -(obsA[0]-meanXaPHuni)*(obsA[0]-meanXaPHuni)/(2*varXaPHuni) - log(stDevXaPHuni*sqrt(2*M_PI)); //log L
		

		if(ap->getUseIntegralApproach() == true) {

			double* xmin = (double *) malloc(1 * sizeof(double));
			double* xmax = (double *) malloc(1 * sizeof(double));
			xmin[0] = ((double) obsAreduced[0]) - 0.5; //Define +/- 0.5 interval around observation
			xmax[0] = ((double) obsAreduced[0]) + 0.5; //Boundaries may be < 0 or > Na
			double* err = (double *) malloc(intOutputDim * sizeof(double));
			err[0] = 9999;	

			parameters params;
			params.mu = meanXaPHuni;
			params.var = varXaPHuni;
			params.stDev = stDevXaPHuni;
			double PxAph = 0; //Not log parameter
		
			//Perform integral using pcubature (best for dim <=3)
			int intOutcome = pcubature(intOutputDim, integrandNormal, &params, 1, xmin, xmax, maxEval, tol, tol, ERROR_INDIVIDUAL, &PxAph, err); 
			bool error = false;
			if(intOutcome != 0) { //Error in integral

				cout << "Error computing numerical integral. Using original (compound) result instead.\n";
				error = true;
			}
			if(err[0] > 0.001) {
				
				cout << "Error in numerical integral is too large. The error is: " << err[0] << ". Using original (compound) result instead.\n";
				error = true;
			}
			double logPxAphInt = log(PxAph);
			if(::isinf(logPxAphInt) != 0) { //integral evaluation is infinity

				cout << "Error when computing numerical integral: log likelihood is infinity. Using original (compound) result instead.\n";
				error = true;
			}

			if(::isnan(logPxAphInt)!=0) { //integral evaluation is nan

				cout << "Error when computing numerical integral: log likelihood is nan. Using original (compound) result instead.\n";
				error = true;
			}
			if(error == false) { //Using integral approach value if no errors were encountered
				logPxAph = logPxAphInt;
				cout << "Using integral result.\n";
			}
			cout << "logPxAph: " << logPxAph << "\n";

		}
	}

	//Check if likelihoods are still negative infinity, which should not happen
	if(logPxAph == -numeric_limits<double>::max() && print == false) {
		cout << "ERROR in Misc::computeLSelTSelGVar: Likelihood for xA is negative infinity. LogPxAph: " << logPxAph <<"\n";

		cout << "Recomputing mean and variance whilst printing out stats:\n";
		computeLSelTSelGVar(obsA,Nt,ap,qBmean, qBvar, hapFitT, hapFitG, contribsPHset, true);

		cout << "-------------------------------------------------------\n";
	}

	//Likelihood is nan, so print out
	if(::isnan(logPxAph) != 0 && print == false) {
		cout << "ERROR in Misc::computeLSelTSelGVar: Likelihood for xA is nan. LogPxAph: " << logPxAph <<"\n";

		cout << "Recomputing mean and variance whilst printing out stats:\n";
		computeLSelTSelGVar(obsA,Nt,ap,qBmean, qBvar, hapFitT, hapFitG, contribsPHset, true);

		cout << "-------------------------------------------------------\n";
	}

	//Clean up before exit
	gsl_vector_free(sigmaGsigmaTqBmeanGSL);
	gsl_vector_free(TsigmaGsigmaTqBmeanGSL);

	gsl_matrix_free(T);
	gsl_matrix_free(alphaNaMTsigmaGsigmaTqBmean);
	gsl_matrix_free(betaMsigmaTqBmean);
	gsl_matrix_free(DsigmaT_qBmean);
	gsl_matrix_free(qBvarDsigmaT_qBmean_Trans);
	gsl_matrix_free(gammaDsigmaT_qBmean_qBvarDsigmaT_qBmean_Trans);
	gsl_matrix_free(DsigmaG_sigmaTqBmean);
	gsl_matrix_free(TDsigmaG_sigmaTqBmean);
	gsl_matrix_free(three);
	gsl_matrix_free(four);
	gsl_matrix_free(qBvarReduced);

	if(print == true) {
			cout << "Returning likelihood " << logPxAph << "\n";
	}

	return logPxAph;
}




double computeLSelTSelGVarOld(vector<int> &obsA, int Nt, AnaParam *ap, vector<double> &qBmean, gsl_matrix* qBvar, vector<double>& hapFitT, vector<double>& hapFitG, vector<vector<int> >* contribsPHset, bool print) {

	double C = ap->getC();

	if(print == true) {
	
		cout << "Printing intermediate steps in computeLselection.\n";
		cout << "First, the input variables:\n";
		cout << "obsA: "; printIntVector(obsA);
		cout << "Nt: " << Nt << "\n";
		cout << "C: " << C << "\n";
		cout << "qBmean: "; printDoubleVector(qBmean);
		cout << "Contribs: \n";
		for(unsigned int i=0; i<contribsPHset->size(); i++) {
			cout << "PH " << i << ": "; printIntVector((*contribsPHset)[i]);
		}


	}

	//Dimensions
	int dimFH = (int) qBmean.size();
	int dimPH = (int) obsA.size();
	

	//Total observations
	int Na = sumOfVector(obsA);
	
	if(print == true) {
		cout << "dimFH: " << dimFH << "\n";
		cout << "dimPH: " << dimPH << "\n";
		cout << "Na: " << Na << "\n";
	}

	//Construct matrix T
	gsl_matrix *T = gsl_matrix_alloc(dimPH,dimFH);
	for(int i=0; i<dimPH; i++) { //Loop over rows

		int counter=0;
		for(int j=0; j<dimFH; j++) { //Loop over columns

			if((*contribsPHset)[i][counter]==j && counter < (int) (*contribsPHset)[i].size()) { //This assumes the contribs are ordered
				gsl_matrix_set(T,i,j,1);
				counter++;
			} else {
				gsl_matrix_set(T,i,j,0);
			}
		}
	}

	if(print == true) {
		cout << "Matrix T: "; printMatrixMathematica(T);
	}


	/*
	* * 
	* * After computation
	* *
	*/


	/*
 	* First we create the mean of xAPH = Na*T*sigmaG[sigmaT[qBmean]]
 	*/ 
	vector<double> sigmaTqBmean;
	double totalT = 0;
	for(int i=0; i<dimFH; i++) {
		totalT += qBmean[i]*exp(hapFitT[i]);
	}
	for(int i=0; i<dimFH; i++) {
		sigmaTqBmean.push_back(qBmean[i]*exp(hapFitT[i])/totalT);
	}

	vector<double> sigmaGsigmaTqBmean;
	double totalG = 0;
	for(int i=0; i<dimFH; i++) {
		totalG += sigmaTqBmean[i]*exp(hapFitG[i]);
	}
	for(int i=0; i<dimFH; i++) {
		sigmaGsigmaTqBmean.push_back(sigmaTqBmean[i]*exp(hapFitG[i])/totalG);
	}

	if(print == true) {
		cout << "sigmaTqBmean: "; printDoubleVector(sigmaTqBmean);
		cout << "sigmaGsigmaTqBmean: "; printDoubleVector(sigmaGsigmaTqBmean);
	}

	gsl_vector *sigmaGsigmaTqBmeanGSL = gsl_vector_alloc(dimFH);
	for(int i=0;i<dimFH;i++) {
		gsl_vector_set(sigmaGsigmaTqBmeanGSL,i,sigmaGsigmaTqBmean[i]);
	}

	//Create TsigmaGsigmaTqBmean through matrix vector multiplication
	gsl_vector *TsigmaGsigmaTqBmeanGSL = gsl_vector_alloc(dimPH);
	gsl_blas_dgemv(CblasNoTrans, 1.0, T, sigmaGsigmaTqBmeanGSL, 0.0, TsigmaGsigmaTqBmeanGSL);  

	//Convert back to C++ vector
	vector<double> meanXaPH;
	vector<double> TsigmaGsigmaTqBmean; //To be used in computing variance
	for(int i=0; i<dimPH; i++) {
		meanXaPH.push_back(Na*(gsl_vector_get(TsigmaGsigmaTqBmeanGSL, i)));
		TsigmaGsigmaTqBmean.push_back(gsl_vector_get(TsigmaGsigmaTqBmeanGSL, i));
	}

	if(print == true) {
		cout << "meanXaPH: "; printDoubleVector(meanXaPH);
	}


	
	/*
	* Next we compute the variance of xAPH = alpha*Na*M(TsigmaG[sigmaT[qBmean]]) + Na*(Na-alpha)*T*D[sigmaG]_sigmaT(qBmean) * (beta*M(sigmaT[qBmean]) + gamma*D[sigmaT]_qBmean * qBvar * (D[sigmaT]_qBmean)^T) * (D[sigmaG]_sigmaT[qBmean])^T * T^T
	*/
	//Constructing the first term alpha*Na*M(TsigmaG[sigmaT[qBmean]])
	gsl_matrix *alphaNaMTsigmaGsigmaTqBmean = gsl_matrix_alloc(dimPH,dimPH);
	constructMatrixM(TsigmaGsigmaTqBmean, alphaNaMTsigmaGsigmaTqBmean);
	
	//Multiplying by the factor alpha*Na:
	double alpha = (Na+C)/((double)(1+C));
	gsl_matrix_scale(alphaNaMTsigmaGsigmaTqBmean, alpha*Na);

	if(print==true) {
		cout << "alpha*Na*MTsigmaGsigmaTqBmean: "; printMatrixMathematica(alphaNaMTsigmaGsigmaTqBmean);
	}
	


	//For constructing the second term Na*(Na-alpha)*T*D[sigmaG]_sigmaT(qBmean) * (beta*M(sigmaT[qBmean]) + gamma*D[sigmaT]_qBmean * qBvar * (D[sigmaT]_qBmean)^T) * (D[sigmaG]_sigmaT[qBmean])^T * T^T
	//we break it up into 4 parts:
	//1) Compute beta*M(sigmaT[qBmean]) + gamma*D[sigmaT]_qBmean * qBvar * (D[sigmaT]_qBmean)^T
	//2) Compute T*D[sigmaG]_sigmaT[qBmean]
	//3) Compute 1*2^T
	//4) Compute Na(Na-alpha)*2*(1*2^T)

	//Doing 1) first
	//1a) Compute betaMsigmaTqBmean
	gsl_matrix *betaMsigmaTqBmean = gsl_matrix_alloc(dimFH, dimFH);
	constructMatrixM(sigmaTqBmean, betaMsigmaTqBmean);
	int growthFactor = ap->getGrowthFactor();
	int Ng=growthFactor*Nt;
	double NtNgDouble = ((double)Nt)*((double)Ng); //This computation gets out of integer range if Nt>4634, so turn into double

	//As Ng=100Nt we can write beta=(Nt+Ng-1)/(Nt*Ng)
	//as (1-1/(100Nt))/Nt
	//double beta = (Nt+Ng-1)/((double)(Nt*Ng));
	//double beta = ( 1 - 1/((double)(100*Nt)) )/((double)Nt);
	double beta = (Nt+Ng-1)/NtNgDouble;

	//double gamma = (Nt*Ng-Nt-Ng+1)/((double)Nt*Ng);
	double gamma = (NtNgDouble-Nt-Ng+1)/NtNgDouble;
	gsl_matrix_scale(betaMsigmaTqBmean, beta);

	if(print == true) {
		cout << "betaMsigmaTqBmean: "; printMatrixMathematica(betaMsigmaTqBmean);
	}

	//1b) Compute D[sigmaT]_qBmean
	gsl_matrix* DsigmaT_qBmean = gsl_matrix_alloc(dimFH,dimFH);
	for(int i=0; i<dimFH; i++) {
		for(int j=0; j<dimFH; j++) {
	
			if(i==j) {
				gsl_matrix_set(DsigmaT_qBmean, i, j, exp(hapFitT[i])/totalT - qBmean[i]*exp(2*hapFitT[i])/(totalT*totalT));
			} else {
				gsl_matrix_set(DsigmaT_qBmean, i, j, -qBmean[i]*exp(hapFitT[i]+hapFitT[j])/(totalT*totalT));
			}
		}
	}	

	//1c) Compute gamma*D[sigmaT]_qBmean * qBvar * (D[sigmaT]_qBmean)^T + betaMsigmaTqBmean
	gsl_matrix* qBvarDsigmaT_qBmean_Trans = gsl_matrix_alloc(dimFH,dimFH);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, qBvar, DsigmaT_qBmean, 0.0, qBvarDsigmaT_qBmean_Trans);
	gsl_matrix* gammaDsigmaT_qBmean_qBvarDsigmaT_qBmean_Trans = gsl_matrix_alloc(dimFH, dimFH);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, gamma, DsigmaT_qBmean, qBvarDsigmaT_qBmean_Trans, 0.0, gammaDsigmaT_qBmean_qBvarDsigmaT_qBmean_Trans);

	//1d) Perform 1a+1b (resulting variable is 1a)	
	gsl_matrix_add(betaMsigmaTqBmean, gammaDsigmaT_qBmean_qBvarDsigmaT_qBmean_Trans);

	//Computing term 2: T*D[sigmaG]_sigmaT[qBmean]
	//2A) Computing D[sigmaG]_sigmaT[qBmean]
	gsl_matrix* DsigmaG_sigmaTqBmean = gsl_matrix_alloc(dimFH,dimFH);
	for(int i=0; i<dimFH; i++) {
		for(int j=0; j<dimFH; j++) {
	
			if(i==j) {
				gsl_matrix_set(DsigmaG_sigmaTqBmean, i, j, exp(hapFitG[i])/totalG - sigmaTqBmean[i]*exp(2*hapFitG[i])/(totalG*totalG));
			} else {
				gsl_matrix_set(DsigmaG_sigmaTqBmean, i, j, -sigmaTqBmean[i]*exp(hapFitG[i]+hapFitG[j])/(totalG*totalG));
			}
		}
	}	

	//2B)Compute T*D product
	gsl_matrix* TDsigmaG_sigmaTqBmean = gsl_matrix_alloc(dimPH,dimFH);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, DsigmaG_sigmaTqBmean, 0.0, TDsigmaG_sigmaTqBmean);

	//3) Computing term 3 = 1*2^T
	gsl_matrix* three = gsl_matrix_alloc(dimFH,dimPH);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, betaMsigmaTqBmean, TDsigmaG_sigmaTqBmean, 0.0, three);

	//4) Computing term four = Na(Na-alpha)*2*3
	gsl_matrix* four = gsl_matrix_alloc(dimPH,dimPH);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, Na*(Na-alpha), TDsigmaG_sigmaTqBmean, three, 0.0, four);


	if(print==true) {
		cout << "Na*(Na-alpha)*T*D[sigmaG]_sigmaT(qBmean) * (beta*M(sigmaT[qBmean]) + gamma*D[sigmaT]_qBmean * qBvar * (D[sigmaT]_qBmean)^T) * (D[sigmaG]_sigmaT[qBmean])^T * T^T:  "; printMatrixMathematica(four);
	}
	
	//Add together
	gsl_matrix_add(alphaNaMTsigmaGsigmaTqBmean, four); //Final variance for after case. Output stored in first matrix.

	
	//Computing likelihoods depend on dimensionality
	double logPxAph = -numeric_limits<double>::max(); //Set to neg infinity as precaution

	if(dimPH > 1) { //Most common case, i.e. multiple partial haps in a set

		/*
		 * Reduce dimensionality by 1 to ensure non-degeneracy
		 */
		//After case
		meanXaPH.pop_back(); //WLOG remove the last PH
		gsl_matrix_view varXaPHreducedView = gsl_matrix_submatrix(alphaNaMTsigmaGsigmaTqBmean, 0,0, dimPH-1, dimPH-1); //Remove row k and column k
		gsl_matrix * varXaPH = gsl_matrix_alloc(dimPH-1, dimPH-1);
		gsl_matrix_memcpy(varXaPH, &(varXaPHreducedView.matrix)); //Copy the reduced matrix view into a "normal" matrix varXaPH
		vector<double> obsAreduced;
		for(int i=0; i<dimPH-1; i++) {
			obsAreduced.push_back(obsA[i]);
		}

		if(print == true) {
			cout << "Reduced quantities: \n";
			cout << "ObsAreduced: "; printDoubleVector(obsAreduced);
			cout << "meanXaPHreduced: "; printDoubleVector(meanXaPH);
			cout << "varXaPHreduced: "; printMatrixMathematica(varXaPH);
		}

		if(dimPH > 2) { //Multivariate approach	

			//Compute likelihood
			logPxAph = logMultivariateNormalPDF(obsAreduced,meanXaPH,varXaPH);
		
        
		}else { //I.e. dimPH==2
	
			//In this case we retrieve the univariate case when reducing dimensionality by one
			double meanXaPHuni = meanXaPH[0];
			double varXaPHuni = gsl_matrix_get(varXaPH,0,0);
			double stDevXaPHuni = sqrt(varXaPHuni);

			logPxAph = -(obsA[0]-meanXaPHuni)*(obsA[0]-meanXaPHuni)/(2*varXaPHuni) - log(stDevXaPHuni*sqrt(2*M_PI)); //log L
		}

		//Deallocate variables created in this scope
		gsl_matrix_free(varXaPH);

	//In theory, can combine dimPH==2 from above with this one, would be slightly faster, although more clear this way
	} else if (dimPH==1) {
		//Case where there is a single partial haplotype in a set.
		//This can be thought of as a case with a single observed partial haplotype
		//and one (or several) unobserved haplotypes.
		//If all loci in the partial haplotype set are variants (i.e. not monomorphic)
		//then there MUST exists at least one other haplotype containing the unobserved alleles,
		//however, this (or these) haplotype(s) are not observed.
		//As such, WLOG we can consider this a case of dimPH=2, which gets reduced to a one dimensional system
		//under reduction, i.e. similar to the dimPH==2 case of above.
				
		double meanXaPHuni = meanXaPH[0];
		double varXaPHuni = gsl_matrix_get(alphaNaMTsigmaGsigmaTqBmean,0,0);
		double stDevXaPHuni = sqrt(varXaPHuni);	

		logPxAph = -(obsA[0]-meanXaPHuni)*(obsA[0]-meanXaPHuni)/(2*varXaPHuni) - log(stDevXaPHuni*sqrt(2*M_PI)); //log L

	} else { //error
	
		cout << "ERROR in computeLSelTSelG: dimPH is <1\n";
	}

	//Check of likelihoods are still negative infinity, which should not happen
	if(logPxAph == -numeric_limits<double>::max() && print == false) {
		cout << "ERROR in Misc::computeLSelTSelGVar: Likelihood for xA is negative infinity. LogPxAph: " << logPxAph <<"\n";

		cout << "Recomputing mean and variance whilst printing out stats:\n";
		computeLSelTSelGVar(obsA,Nt,ap,qBmean, qBvar, hapFitT, hapFitG, contribsPHset, true);

		cout << "-------------------------------------------------------\n";
	}


	//Clean up before exit
	gsl_vector_free(sigmaGsigmaTqBmeanGSL);
	gsl_vector_free(TsigmaGsigmaTqBmeanGSL);

	gsl_matrix_free(T);
	gsl_matrix_free(alphaNaMTsigmaGsigmaTqBmean);
	gsl_matrix_free(betaMsigmaTqBmean);
	gsl_matrix_free(DsigmaT_qBmean);
	gsl_matrix_free(qBvarDsigmaT_qBmean_Trans);
	gsl_matrix_free(gammaDsigmaT_qBmean_qBvarDsigmaT_qBmean_Trans);
	gsl_matrix_free(DsigmaG_sigmaTqBmean);
	gsl_matrix_free(TDsigmaG_sigmaTqBmean);
	gsl_matrix_free(three);
	gsl_matrix_free(four);


	return logPxAph;
}






double computeLNeutral(vector<int> &obsB, vector<int> &obsA, int Nt, AnaParam *ap, vector<double> &qBfhFreqs, vector<vector<int> >* contribsPHset, bool print) {

	double C = ap->getC();
	vector<double>& factStore = *(ap->getLogFactStore()); //vector is not copied here! Saves much time
    
	if(print == true) {
	
		cout << "Printing intermediate steps in computeLselection.\n";
		cout << "First, the input variables:\n";
		cout << "obsB: "; printIntVector(obsB);
		cout << "obsA: "; printIntVector(obsA);
		cout << "Nt: " << Nt << "\n";
		cout << "C: " << C << "\n";
		cout << "qBfhFreqs: "; printDoubleVector(qBfhFreqs);
		cout << "Contribs: \n";
		for(unsigned int i=0; i<contribsPHset->size(); i++) {
			cout << "PH " << i << ": "; printIntVector((*contribsPHset)[i]);
		}


	}

	//Dimensions
	int dimFH = (int) qBfhFreqs.size();
	int dimPH = (int) obsB.size();
	

	//Total observations
	int Na = sumOfVector(obsA);
	
	if(print == true) {
		cout << "dimFH: " << dimFH << "\n";
		cout << "dimPH: " << dimPH << "\n";
		cout << "Na: " << Na << "\n";
	}

	//Construct matrix T
	gsl_matrix *T = gsl_matrix_alloc(dimPH,dimFH);
	for(int i=0; i<dimPH; i++) { //Loop over rows

		int counter=0;
		for(int j=0; j<dimFH; j++) { //Loop over columns

			if((*contribsPHset)[i][counter]==j && counter < (int) (*contribsPHset)[i].size()) { //This assumes the contribs are ordered
				gsl_matrix_set(T,i,j,1);
				counter++;
			} else {
				gsl_matrix_set(T,i,j,0);
			}
		}
	}

	if(print == true) {
		cout << "Matrix T: "; printMatrixMathematica(T);
	}


	/*
	* * 
	* * Before computation
	* *
	*/

	/*
	 * First, transform qB into qB_PH=TqB
	 */
	//Turn qB into gsl vector
	gsl_vector *qBGSL = gsl_vector_alloc(dimFH);
	for(int i=0;i<dimFH;i++) {
		gsl_vector_set(qBGSL,i,qBfhFreqs[i]);
	}

	//Create TqB through matrix vector multiplication
	gsl_vector *TqBGSL = gsl_vector_alloc(dimPH);
	gsl_blas_dgemv(CblasNoTrans, 1.0, T, qBGSL, 0.0, TqBGSL);

	//Convert back into C++ vector
	vector<double> TqB;
	for(int i=0; i<dimPH; i++) {
		TqB.push_back(gsl_vector_get(TqBGSL, i));
	}  

	if(print == true) {
		cout << "TqB: "; printDoubleVector(TqB);
	}
	
	
	/*
	* Finally, compute log likelihood from Dirichet multinomial PDF with dispersion parameter C
	*/
	double logPxBph = logDirMultProbC(C, obsB, TqB, factStore);
	
	if(print == true) {
		cout << "logPxBph = " << logPxBph << "\n";
	}



	/*
	* * 
	* * After computation
	* *
	*/


	/*
 	* First we create the mean of xAPH = Na*T*qB
 	*/ 
	//Scale vector
	vector<double> meanXaPH;
	for(int i=0; i<dimPH; i++) {

		meanXaPH.push_back(Na*TqB[i]);
	}
	
	/*
	* Next we compute the variance of xAPH = Na(alpha+(Na-alpha)beta)M(TqB)
	*/
	//Constructing the matrix M
	gsl_matrix *MTqB = gsl_matrix_alloc(dimPH,dimPH);
	constructMatrixM(TqB, MTqB);
	
	//Multiplying by the factor Na(alpha+(Na-alpha)beta):
	int growthFactor = ap->getGrowthFactor();
	int Ng=growthFactor*Nt;
	double NtNgDouble = ((double)Nt)*((double)Ng); //This computation gets out of integer range if Nt>4634, so turn into double
	double alpha = (Na+C)/(double)(1+C);
	double beta = (Nt+Ng-1)/NtNgDouble;
	double NaAlphaPlusNaMinusAlphaBeta = Na*(alpha+(Na-alpha)*beta);
	gsl_matrix_scale(MTqB, NaAlphaPlusNaMinusAlphaBeta);

	if(print==true) {
		cout << "Na(alpha+(Na-alpha)beta)MTqB: "; printMatrixMathematica(MTqB);
	}


	//Reduce dimensionality by 1 to ensure non-degeneracy
	meanXaPH.pop_back(); //WLOG remove the last PH
	gsl_matrix_view varXaPHreducedView = gsl_matrix_submatrix(MTqB, 0,0, dimPH-1, dimPH-1); //Remove row k and column k
	gsl_matrix * varXaPH = gsl_matrix_alloc(dimPH-1, dimPH-1);
	gsl_matrix_memcpy(varXaPH, &(varXaPHreducedView.matrix)); //Copy the reduced matrix view into a "normal" matrix varXaPH
	vector<double> obsAreduced;
	for(int i=0; i<dimPH-1; i++) {
		obsAreduced.push_back(obsA[i]);
	}

	if(print == true) {
		cout << "Reduced quantities: \n";
		cout << "ObsAreduced: "; printDoubleVector(obsAreduced);
		cout << "meanXaPHreduced: "; printDoubleVector(meanXaPH);
		cout << "varXaPHreduced: "; printMatrixMathematica(varXaPH);
	}

	double logPxAph = -numeric_limits<double>::max(); //Set to neg infinity as precaution

	if(dimPH > 2) { //Multivariate approach	

		//Compute likelihood
		logPxAph = logMultivariateNormalPDF(obsAreduced,meanXaPH,varXaPH);
    
		
        
	}else if(dimPH == 2) {
	
		//In this case we retrieve the univariate case when reducing dimensionality by one
		double meanXaPHuni = meanXaPH[0];
		double varXaPHuni = gsl_matrix_get(varXaPH,0,0);
		double stDevXaPHuni = sqrt(varXaPHuni);

		logPxAph = -(obsA[0]-meanXaPHuni)*(obsA[0]-meanXaPHuni)/(2*varXaPHuni) - log(stDevXaPHuni*sqrt(2*M_PI)); //log L


	} else { //error
	
		cout << "ERROR in computeLSelection: dimPH is <2\n";
		logPxAph = -numeric_limits<double>::max();
	}

	if((logPxAph >= 0 || logPxAph == -numeric_limits<double>::max()) && print == false) {
		cout << "ERROR in Misc::computeLSelection: Likelihood is positive or negative infinity. Likelihood is: " << logPxAph  << "\n";

		if(dimPH>2) {
			cout << "Recomputing likelihood whilst printing out stats:\n";
			logMultivariateNormalPDFPrint(obsAreduced,meanXaPH,varXaPH);
		}

		cout << "-------------------------------------------------------\n";
		cout << "Recomputing mean and variance whilst printing out stats:\n";
		computeLNeutral(obsB,obsA,Nt,ap,qBfhFreqs, contribsPHset, true);

		//Clean up and return negative infinity
		gsl_vector_free(qBGSL);
		gsl_vector_free(TqBGSL);

		gsl_matrix_free(T);
		gsl_matrix_free(MTqB);
		gsl_matrix_free(varXaPH);

		cout << "Returning likelihood of negative infinity.\n\n";
		return -numeric_limits<double>::max();
	}

	//Clean up before exit
	gsl_vector_free(qBGSL);
	gsl_vector_free(TqBGSL);

	gsl_matrix_free(T);
	gsl_matrix_free(MTqB);
	gsl_matrix_free(varXaPH);

	//cout << "logPxAph = " << logPxAph << "\n";
	return logPxBph + logPxAph;
}

//Computes likelyhood L(Nt,sigma|qB,xA) given qB (represented by mean and var)
double computeLNeutralVar(vector<int> &obsA, int Nt, AnaParam *ap, vector<double> &qBmean, gsl_matrix* qBvar, vector<vector<int> >* contribsPHset, bool print) {

	double C = ap->getC();

	if(print == true) {
	
		cout << "Printing intermediate steps in computeLselection.\n";
		cout << "First, the input variables:\n";
		cout << "obsA: "; printIntVector(obsA);
		cout << "Nt: " << Nt << "\n";
		cout << "C: " << C << "\n";
		cout << "qBmean: "; printDoubleVector(qBmean);
		cout << "qBvar: "; printMatrixMathematica(qBvar);
		cout << "Contribs: \n";
		for(unsigned int i=0; i<contribsPHset->size(); i++) {
			cout << "PH " << i << ": "; printIntVector((*contribsPHset)[i]);
		}


	}

	//Dimensions
	int dimFH = (int) qBmean.size();
	int dimPH = (int) obsA.size();
	

	//Total observations
	int Na = sumOfVector(obsA);
	
	if(print == true) {
		cout << "dimFH: " << dimFH << "\n";
		cout << "dimPH: " << dimPH << "\n";
		cout << "Na: " << Na << "\n";
	}

	//Construct matrix T
	gsl_matrix *T = gsl_matrix_alloc(dimPH,dimFH);
	for(int i=0; i<dimPH; i++) { //Loop over rows

		int counter=0;
		for(int j=0; j<dimFH; j++) { //Loop over columns

			if((*contribsPHset)[i][counter]==j && counter < (int) (*contribsPHset)[i].size()) { //This assumes the contribs are ordered
				gsl_matrix_set(T,i,j,1);
				counter++;
			} else {
				gsl_matrix_set(T,i,j,0);
			}
		}
	}

	if(print == true) {
		cout << "Matrix T: "; printMatrixMathematica(T);
	}

	/*
	* * 
	* * After computation
	* *
	*/


	/*
 	* First we create the mean of xAPH = Na*T*qBmean
 	*/ 

	//Turn qBmean into gsl vector
	gsl_vector *qBmeanGSL = gsl_vector_alloc(dimFH);
	for(int i=0;i<dimFH;i++) {
		gsl_vector_set(qBmeanGSL,i,qBmean[i]);
	}

	//Create TqBmean through matrix vector multiplication
	gsl_vector *TqBmeanGSL = gsl_vector_alloc(dimPH);
	gsl_blas_dgemv(CblasNoTrans, 1.0, T, qBmeanGSL, 0.0, TqBmeanGSL);

	//Convert back into C++ vector
	vector<double> TqBmean; //To be used later
	vector<double> meanXaPH;
	for(int i=0; i<dimPH; i++) {
		TqBmean.push_back(gsl_vector_get(TqBmeanGSL, i));
		meanXaPH.push_back(Na*TqBmean[i]);
	}  

	if(print == true) {
		cout << "TqB: "; printDoubleVector(TqBmean);
	}

	
	/*
	* Next we compute the variance of xAPH = Na(alpha+(Na-alpha)beta)M(TqBmean) + Na(Na-alpha)gammaTqBvarTTrans
	*/
	//Constructing the first term
	gsl_matrix *NaAlphaPlusNaMinusAlphaBetaMTqBmean = gsl_matrix_alloc(dimPH,dimPH);
	constructMatrixM(TqBmean, NaAlphaPlusNaMinusAlphaBetaMTqBmean);

	if(print==true) {
		cout << "M(TqBmean): "; printMatrixMathematica(NaAlphaPlusNaMinusAlphaBetaMTqBmean);
	}
	
	//Multiplying by the factor Na(alpha+(Na-alpha)beta):
	int growthFactor = ap->getGrowthFactor();
	int Ng=growthFactor*Nt;
	if(print == true) {
		cout << "Growth factor: " << growthFactor << "\n";
		cout << "Nt: " << Nt << "\n";
		cout << "Ng: " << Ng << "\n";
	}
	double NtNgDouble = ((double)Nt)*((double)Ng); //This computation gets out of integer range if Nt>4634, so turn into double
	double alpha = (Na+C)/((double)(1+C));

	double beta = (Nt+Ng-1)/NtNgDouble;
	double NaAlphaPlusNaMinusAlphaBeta = Na*(alpha+(Na-alpha)*beta);
	if(print==true) {
		cout << "NaAlphaPlusNaMinusAlphaBeta: " << NaAlphaPlusNaMinusAlphaBeta << "\n";
	}
	gsl_matrix_scale(NaAlphaPlusNaMinusAlphaBetaMTqBmean, NaAlphaPlusNaMinusAlphaBeta);

	if(print==true) {
		cout << "Na(alpha+(Na-alpha)beta)MTqB: "; printMatrixMathematica(NaAlphaPlusNaMinusAlphaBetaMTqBmean);
	}
	
	//Creating TqBvarTTrans
	gsl_matrix* TqBvar = gsl_matrix_alloc(dimPH,dimFH);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, qBvar, 0.0, TqBvar);
	gsl_matrix* TqBvarTTrans = gsl_matrix_alloc(dimPH,dimPH);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, TqBvar, T, 0.0, TqBvarTTrans);

	//Constructing the second term
	double gamma = (NtNgDouble-Nt-Ng+1)/NtNgDouble;
	gsl_matrix* NaNaMinusAlphaGammaTqBvarTTrans = gsl_matrix_alloc(dimPH,dimPH);
	gsl_matrix_memcpy(NaNaMinusAlphaGammaTqBvarTTrans, TqBvarTTrans);
	gsl_matrix_scale(NaNaMinusAlphaGammaTqBvarTTrans, Na*(Na-alpha)*gamma);

	//If qBvar is used, add together
	if(ap->getNoVar() == false && ap->getMeanOnly() == false) {

		//Add together
		gsl_matrix_add(NaAlphaPlusNaMinusAlphaBetaMTqBmean, NaNaMinusAlphaGammaTqBvarTTrans); //Final variance for after case. Output stored in first matrix.
		if(print==true) {
			cout << "Na(alpha+(Na-alpha)beta)M(TqBmean) + Na(Na-alpha)gamma*T*qBvar*T^Trans:"; printMatrixMathematica(NaAlphaPlusNaMinusAlphaBetaMTqBmean);
		}
	} else {

		//Don't add together, i.e. NaAlphaPlusNaMinusAlphaBetaMTqBmean remains unchanged
	}


	//Computing likelihoods depend on dimensionality
	double logPxAph = -numeric_limits<double>::max(); //Set to neg infinity as precaution
	if(dimPH > 2) { //I.e. most common case with multiple partial haps in a set

		/*
		 * Reduce dimensionality by 1 to ensure non-degeneracy
		 */
		//After case
		meanXaPH.pop_back(); //WLOG remove the last PH
		gsl_matrix_view varXaPHreducedView = gsl_matrix_submatrix(NaAlphaPlusNaMinusAlphaBetaMTqBmean, 0,0, dimPH-1, dimPH-1); //Remove row k and column k
		gsl_matrix * varXaPH = gsl_matrix_alloc(dimPH-1, dimPH-1);
		gsl_matrix_memcpy(varXaPH, &(varXaPHreducedView.matrix)); //Copy the reduced matrix view into a "normal" matrix varXaPH
		vector<double> obsAreduced;
		for(int i=0; i<dimPH-1; i++) {
			obsAreduced.push_back(obsA[i]);
		}
	
		if(print == true) {
			cout << "Reduced quantities: \n";

			cout << "ObsAreduced: "; printDoubleVector(obsAreduced);
			cout << "meanXaPHreduced: "; printDoubleVector(meanXaPH);
			cout << "varXaPHreduced: "; printMatrixMathematica(varXaPH);
		}
	
        	//Compute likelihood
		logPxAph = logMultivariateNormalPDF(obsAreduced,meanXaPH,varXaPH);
		if(print==true) {
			logMultivariateNormalPDFPrint(obsAreduced,meanXaPH,varXaPH);
		}

		
		if(ap->getUseIntegralApproach() == true) {
			double* xmin = (double *) malloc((dimPH-1) * sizeof(double));
			double* xmax = (double *) malloc((dimPH-1) * sizeof(double));
			double* err = (double *) malloc(intOutputDim * sizeof(double));
			err[0] = 9999;
			for(int i=0; i<(dimPH-1); i++) {
		
				xmin[i] = ((double) obsAreduced[i]) - 0.5; //Define hypercube around observation
				xmax[i] = ((double) obsAreduced[i]) + 0.5; //Boundaries may be < 0 or > Na
			}

			parameters params;
			params.muVec = meanXaPH;
			params.Sigma = varXaPH;
			double PxAph = 0; //Not log parameter
		
			//Perform integral. For dim <= 3, pcubature is best. hcubature is better for larger dimensions, although best for <7.
			int intOutcome = -1;
			if((dimPH-1) <= 3) {
				intOutcome = pcubature(intOutputDim, integrandMultiNormal, &params, dimPH-1, xmin, xmax, maxEval, tol, tol, ERROR_INDIVIDUAL, &PxAph, err); 
			} else {
				intOutcome = hcubature(intOutputDim, integrandMultiNormal, &params, dimPH-1, xmin, xmax, maxEval, tol, tol, ERROR_INDIVIDUAL, &PxAph, err); 
			}

			bool error = false;
			if(intOutcome != 0) { //Error in integral

				cout << "Error computing numerical integral. Using original (compound) result instead.\n";
				error = true;
			}
			if(err[0] > 0.001) {
				
				cout << "Error in numerical integral is too large. The error is: " << err[0] << ". Using original (compound) result instead.\n";
				error = true;
			}
			double logPxAphInt = log(PxAph);
			if(::isinf(logPxAphInt) != 0) { //integral evaluation is infinity

				cout << "Error when computing numerical integral: log likelihood is infinity. Using original (compound) result instead.\n";
				error = true;
			}

			if(::isnan(logPxAphInt)!=0) { //integral evaluation is nan

				cout << "Error when computing numerical integral: log likelihood is nan. Using original (compound) result instead.\n";
				error = true;
			}
			if(error == false) { //Using integral approach value if no errors were encountered
				logPxAph = logPxAphInt;
				cout << "Using integral result.\n";
			}
			cout << "logPxAph: " << logPxAph << "\n";
			
		}


		//Deallocate variables created in this scope
		gsl_matrix_free(varXaPH);

	} else { //I.e. dimPH==2 or dimPH==1 (same computation)
		//Case where there is a single partial haplotype in a set.
		//This can be thought of as a case with a single observed partial haplotype
		//and one (or several) unobserved haplotypes.
		//If all loci in the partial haplotype set are variants (i.e. not monomorphic)
		//then there MUST exists at least one other haplotype containing the unobserved alleles,
		//however, this (or these) haplotype(s) are not observed.
		//As such, WLOG we can consider this a case of dimPH=2, which gets reduced to a one dimensional system
		//under reduction, i.e. similar to the dimPH==2 case of above.
		double meanXaPHuni = meanXaPH[0];
		double varXaPHuni = gsl_matrix_get(NaAlphaPlusNaMinusAlphaBetaMTqBmean,0,0);
		double stDevXaPHuni = sqrt(varXaPHuni);	

		logPxAph = -(obsA[0]-meanXaPHuni)*(obsA[0]-meanXaPHuni)/(2*varXaPHuni) - log(stDevXaPHuni*sqrt(2*M_PI)); //log L
		if(print==true) {
			cout << "MeanXaPHuni: " << meanXaPHuni << "\n";
			cout << "VarXaPHuni: " << varXaPHuni << "\n";
			cout << "logPxAph: " << logPxAph << "\n";
		}
		

		if(ap->getUseIntegralApproach() == true) {

			double* xmin = (double *) malloc(1 * sizeof(double));
			double* xmax = (double *) malloc(1 * sizeof(double));
			xmin[0] = ((double) obsA[0]) - 0.5; //Define +/- 0.5 interval around observation
			xmax[0] = ((double) obsA[0]) + 0.5; //Boundaries may be < 0 or > Na
			double* err = (double *) malloc(intOutputDim * sizeof(double));
			err[0] = 9999;

			parameters params;
			params.mu = meanXaPHuni;
			params.var = varXaPHuni;
			params.stDev = stDevXaPHuni;
			double PxAph = 0;
		
			//Perform integral using pcubature (best for dim <=3)
			int intOutcome = pcubature(intOutputDim, integrandNormal, &params, 1, xmin, xmax, maxEval, tol, tol, ERROR_INDIVIDUAL, &PxAph, err); 
			bool error = false;
			if(intOutcome != 0) { //Error in integral

				cout << "Error computing numerical integral. Using original (compound) result instead.\n";
				error = true;
			}
			if(err[0] > 0.001) {
				
				cout << "Error in numerical integral is too large. The error is: " << err[0] << ". Using original (compound) result instead.\n";
				error = true;
			}
			double logPxAphInt = log(PxAph);
			if(::isinf(logPxAphInt) != 0) { //integral evaluation is infinity

				cout << "Error when computing numerical integral: log likelihood is infinity. Using original (compound) result instead.\n";
				error = true;
			}

			if(::isnan(logPxAphInt)!=0) { //integral evaluation is nan

				cout << "Error when computing numerical integral: log likelihood is nan. Using original (compound) result instead.\n";
				error = true;
			}
			if(error == false) { //Using integral approach value if no errors were encountered
				logPxAph = logPxAphInt;
				cout << "Using integral result.\n";
			}
			cout << "logPxAph: " << logPxAph << "\n";
		}
	}
	
	if(logPxAph == -numeric_limits<double>::max() && print == false) {
		cout << "ERROR in Misc::computeLNeutralVar: Likelihood for xA is negative infinity. LogPxAph: " << logPxAph <<"\n";

		cout << "Recomputing mean and variance whilst printing out stats:\n";
		computeLNeutralVar(obsA,Nt,ap,qBmean, qBvar, contribsPHset, true);

		cout << "-------------------------------------------------------\n";
	}

	//Likelihood is nan, so print out
	if(::isnan(logPxAph) != 0 && print == false) {
		cout << "ERROR in Misc::computeLNeutralVarOrignal: Likelihood for xA is nan. LogPxAph: " << logPxAph <<"\n";

		cout << "Recomputing mean and variance whilst printing out stats:\n";
		computeLNeutralVar(obsA,Nt,ap,qBmean, qBvar, contribsPHset, true);

		cout << "-------------------------------------------------------\n";
	}

	//Clean up before exit
	gsl_vector_free(qBmeanGSL);
	gsl_vector_free(TqBmeanGSL);

	gsl_matrix_free(T);
	gsl_matrix_free(TqBvar);
	gsl_matrix_free(TqBvarTTrans);
	gsl_matrix_free(NaAlphaPlusNaMinusAlphaBetaMTqBmean);
	gsl_matrix_free(NaNaMinusAlphaGammaTqBvarTTrans);

	return logPxAph;
}



//Computes likelyhood L(Nt,sigma|qB,xA) given qB (represented by mean and var)
double computeLNeutralVarAdvanced(vector<int> &obsA, int& deltaDays, int Nt, AnaParam *ap, vector<double> &qBmean, gsl_matrix* qBvar, vector<vector<int> >* contribsPHset, bool print) {

	double C = ap->getC();

	if(print == true) {
	
		cout << "Printing intermediate steps in computeLselection.\n";
		cout << "First, the input variables:\n";
		cout << "obsA: "; printIntVector(obsA);
		cout << "deltaDays: " << deltaDays << "\n";
		cout << "Nt: " << Nt << "\n";
		cout << "C: " << C << "\n";
		cout << "qBmean: "; printDoubleVector(qBmean);
		cout << "qBvar: "; printMatrixMathematica(qBvar);
		cout << "Contribs: \n";
		for(unsigned int i=0; i<contribsPHset->size(); i++) {
			cout << "PH " << i << ": "; printIntVector((*contribsPHset)[i]);
		}


	}

	//Dimensions
	int dimFH = (int) qBmean.size();
	int dimPH = (int) obsA.size();
	

	//Total observations
	int Na = sumOfVector(obsA);
	
	if(print == true) {
		cout << "dimFH: " << dimFH << "\n";
		cout << "dimPH: " << dimPH << "\n";
		cout << "Na: " << Na << "\n";
	}

	//Construct matrix T
	gsl_matrix *T = gsl_matrix_alloc(dimPH,dimFH);
	for(int i=0; i<dimPH; i++) { //Loop over rows

		int counter=0;
		for(int j=0; j<dimFH; j++) { //Loop over columns

			if((*contribsPHset)[i][counter]==j && counter < (int) (*contribsPHset)[i].size()) { //This assumes the contribs are ordered
				gsl_matrix_set(T,i,j,1);
				counter++;
			} else {
				gsl_matrix_set(T,i,j,0);
			}
		}
	}

	if(print == true) {
		cout << "Matrix T: "; printMatrixMathematica(T);
	}

	/*
	* * 
	* * After computation
	* *
	*/


	/*
 	* First we create the mean of xAPH = Na*T*qBmean
 	*/ 

	//Turn qBmean into gsl vector
	gsl_vector *qBmeanGSL = gsl_vector_alloc(dimFH);
	for(int i=0;i<dimFH;i++) {
		gsl_vector_set(qBmeanGSL,i,qBmean[i]);
	}

	//Create TqBmean through matrix vector multiplication
	gsl_vector *TqBmeanGSL = gsl_vector_alloc(dimPH);
	gsl_blas_dgemv(CblasNoTrans, 1.0, T, qBmeanGSL, 0.0, TqBmeanGSL);

	//Convert back into C++ vector
	vector<double> TqBmean; //To be used later
	vector<double> meanXaPH;
	for(int i=0; i<dimPH; i++) {
		TqBmean.push_back(gsl_vector_get(TqBmeanGSL, i));
		meanXaPH.push_back(Na*TqBmean[i]);
	}  

	if(print == true) {
		cout << "TqB: "; printDoubleVector(TqBmean);
	}

	
	/*
	* Next we compute the variance of xAPH = Na(alpha+(Na-alpha)(gamma_N))M(TqBmean) + Na(Na-alpha)(delta_N)TqBvarTTrans
	*/
	//Constructing the first term
	gsl_matrix *NaAlphaPlusNaMinusAlphaGammaMTqBmean = gsl_matrix_alloc(dimPH,dimPH);
	constructMatrixM(TqBmean, NaAlphaPlusNaMinusAlphaGammaMTqBmean);

	if(print==true) {
		cout << "M(TqBmean): "; printMatrixMathematica(NaAlphaPlusNaMinusAlphaGammaMTqBmean);
	}
	
	//Multiplying by the factor Na(alpha+(Na-alpha)gamma):
	double alpha = (Na+C)/((double)(1+C));
	int growthFactor = ap->getGrowthFactor();
	if(print == true) {
		cout << "Growth factor: " << growthFactor << "\n";
	}

	//gamma_N = gamma_{deltaDays*2}, e.g. if deltaDays=2, then gamma_N = gamma_4
	double gammaN = gamma_n(deltaDays*2, growthFactor, Nt); //Defined in misc
	if(print == true) {
		cout << "gammaN: " << gammaN << "\n";
	}
	double NaAlphaPlusNaMinusAlphaGamma = Na*(alpha+(Na-alpha)*gammaN);
	if(print==true) {
		cout << "NaAlphaPlusNaMinusAlphaGamma: " << NaAlphaPlusNaMinusAlphaGamma << "\n";
	}
	gsl_matrix_scale(NaAlphaPlusNaMinusAlphaGammaMTqBmean, NaAlphaPlusNaMinusAlphaGamma);

	if(print==true) {
		cout << "Na(alpha+(Na-alpha)gamma)MTqB: "; printMatrixMathematica(NaAlphaPlusNaMinusAlphaGammaMTqBmean);
	}
	
	//Creating TqBvarTTrans
	gsl_matrix* TqBvar = gsl_matrix_alloc(dimPH,dimFH);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, qBvar, 0.0, TqBvar);
	gsl_matrix* TqBvarTTrans = gsl_matrix_alloc(dimPH,dimPH);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, TqBvar, T, 0.0, TqBvarTTrans);

	//delta_N = delta_{deltaDays*2}, e.g. if deltaDays=2, then delta_N = delta_4
	double deltaN = delta_n(deltaDays*2, growthFactor, Nt);
	if(print == true) {
		cout << "delta_n: " << deltaN << "\n";
	}


	//Constructing the second term
	gsl_matrix* NaNaMinusAlphaDeltaTqBvarTTrans = gsl_matrix_alloc(dimPH,dimPH);
	gsl_matrix_memcpy(NaNaMinusAlphaDeltaTqBvarTTrans, TqBvarTTrans);
	gsl_matrix_scale(NaNaMinusAlphaDeltaTqBvarTTrans, Na*(Na-alpha)*deltaN);

	//If qBvar is used, add together
	if(ap->getNoVar() == false && ap->getMeanOnly() == false) {

		//Add together
		gsl_matrix_add(NaAlphaPlusNaMinusAlphaGammaMTqBmean, NaNaMinusAlphaDeltaTqBvarTTrans); //Final variance for after case. Output stored in first matrix.
		if(print==true) {
			cout << "Na(alpha+(Na-alpha)gamma)M(TqBmean) + Na(Na-alpha)delta*T*qBvar*T^Trans:"; printMatrixMathematica(NaAlphaPlusNaMinusAlphaGammaMTqBmean);
		}
	} else {

		//Don't add together, i.e. NaAlphaPlusNaMinusAlphaGammaMTqBmean remains unchanged
	}


	//Computing likelihoods depend on dimensionality
	double logPxAph = -numeric_limits<double>::max(); //Set to neg infinity as precaution
	if(dimPH > 2) { //I.e. most common case with multiple partial haps in a set

		/*
		 * Reduce dimensionality by 1 to ensure non-degeneracy
		 */
		//After case
		meanXaPH.pop_back(); //WLOG remove the last PH
		gsl_matrix_view varXaPHreducedView = gsl_matrix_submatrix(NaAlphaPlusNaMinusAlphaGammaMTqBmean, 0,0, dimPH-1, dimPH-1); //Remove row k and column k
		gsl_matrix * varXaPH = gsl_matrix_alloc(dimPH-1, dimPH-1);
		gsl_matrix_memcpy(varXaPH, &(varXaPHreducedView.matrix)); //Copy the reduced matrix view into a "normal" matrix varXaPH
		vector<double> obsAreduced;
		for(int i=0; i<dimPH-1; i++) {
			obsAreduced.push_back(obsA[i]);
		}
	
		if(print == true) {
			cout << "Reduced quantities: \n";

			cout << "ObsAreduced: "; printDoubleVector(obsAreduced);
			cout << "meanXaPHreduced: "; printDoubleVector(meanXaPH);
			cout << "varXaPHreduced: "; printMatrixMathematica(varXaPH);
		}
	
        	//Compute likelihood
		logPxAph = logMultivariateNormalPDF(obsAreduced,meanXaPH,varXaPH);
		if(print==true) {
			logMultivariateNormalPDFPrint(obsAreduced,meanXaPH,varXaPH);
		}

		
		if(ap->getUseIntegralApproach() == true) {
			double* xmin = (double *) malloc((dimPH-1) * sizeof(double));
			double* xmax = (double *) malloc((dimPH-1) * sizeof(double));
			double* err = (double *) malloc(intOutputDim * sizeof(double));
			err[0] = 9999;
			for(int i=0; i<(dimPH-1); i++) {
		
				xmin[i] = ((double) obsAreduced[i]) - 0.5; //Define hypercube around observation
				xmax[i] = ((double) obsAreduced[i]) + 0.5; //Boundaries may be < 0 or > Na
			}

			parameters params;
			params.muVec = meanXaPH;
			params.Sigma = varXaPH;
			double PxAph = 0; //Not log parameter
		
			//Perform integral. For dim <= 3, pcubature is best. hcubature is better for larger dimensions, although best for <7.
			int intOutcome = -1;
			if((dimPH-1) <= 3) {
				intOutcome = pcubature(intOutputDim, integrandMultiNormal, &params, dimPH-1, xmin, xmax, maxEval, tol, tol, ERROR_INDIVIDUAL, &PxAph, err); 
			} else {
				intOutcome = hcubature(intOutputDim, integrandMultiNormal, &params, dimPH-1, xmin, xmax, maxEval, tol, tol, ERROR_INDIVIDUAL, &PxAph, err); 
			}

			bool error = false;
			if(intOutcome != 0) { //Error in integral

				cout << "Error computing numerical integral. Using original (compound) result instead.\n";
				error = true;
			}
			if(err[0] > 0.001) {
				
				cout << "Error in numerical integral is too large. The error is: " << err[0] << ". Using original (compound) result instead.\n";
				error = true;
			}
			double logPxAphInt = log(PxAph);
			if(::isinf(logPxAphInt) != 0) { //integral evaluation is infinity

				cout << "Error when computing numerical integral: log likelihood is infinity. Using original (compound) result instead.\n";
				error = true;
			}

			if(::isnan(logPxAphInt)!=0) { //integral evaluation is nan

				cout << "Error when computing numerical integral: log likelihood is nan. Using original (compound) result instead.\n";
				error = true;
			}
			if(error == false) { //Using integral approach value if no errors were encountered
				logPxAph = logPxAphInt;
				cout << "Using integral result.\n";
			}
			cout << "logPxAph: " << logPxAph << "\n";
			
		}


		//Deallocate variables created in this scope
		gsl_matrix_free(varXaPH);

	} else { //I.e. dimPH==2 or dimPH==1 (same computation)
		//Case where there is a single partial haplotype in a set.
		//This can be thought of as a case with a single observed partial haplotype
		//and one (or several) unobserved haplotypes.
		//If all loci in the partial haplotype set are variants (i.e. not monomorphic)
		//then there MUST exists at least one other haplotype containing the unobserved alleles,
		//however, this (or these) haplotype(s) are not observed.
		//As such, WLOG we can consider this a case of dimPH=2, which gets reduced to a one dimensional system
		//under reduction, i.e. similar to the dimPH==2 case of above.
		double meanXaPHuni = meanXaPH[0];
		double varXaPHuni = gsl_matrix_get(NaAlphaPlusNaMinusAlphaGammaMTqBmean,0,0);
		double stDevXaPHuni = sqrt(varXaPHuni);	

		logPxAph = -(obsA[0]-meanXaPHuni)*(obsA[0]-meanXaPHuni)/(2*varXaPHuni) - log(stDevXaPHuni*sqrt(2*M_PI)); //log L
		if(print==true) {
			cout << "MeanXaPHuni: " << meanXaPHuni << "\n";
			cout << "VarXaPHuni: " << varXaPHuni << "\n";
			cout << "logPxAph: " << logPxAph << "\n";
		}
		

		if(ap->getUseIntegralApproach() == true) {

			double* xmin = (double *) malloc(1 * sizeof(double));
			double* xmax = (double *) malloc(1 * sizeof(double));
			xmin[0] = ((double) obsA[0]) - 0.5; //Define +/- 0.5 interval around observation
			xmax[0] = ((double) obsA[0]) + 0.5; //Boundaries may be < 0 or > Na
			double* err = (double *) malloc(intOutputDim * sizeof(double));
			err[0] = 9999;

			parameters params;
			params.mu = meanXaPHuni;
			params.var = varXaPHuni;
			params.stDev = stDevXaPHuni;
			double PxAph = 0;
		
			//Perform integral using pcubature (best for dim <=3)
			int intOutcome = pcubature(intOutputDim, integrandNormal, &params, 1, xmin, xmax, maxEval, tol, tol, ERROR_INDIVIDUAL, &PxAph, err); 
			bool error = false;
			if(intOutcome != 0) { //Error in integral

				cout << "Error computing numerical integral. Using original (compound) result instead.\n";
				error = true;
			}
			if(err[0] > 0.001) {
				
				cout << "Error in numerical integral is too large. The error is: " << err[0] << ". Using original (compound) result instead.\n";
				error = true;
			}
			double logPxAphInt = log(PxAph);
			if(::isinf(logPxAphInt) != 0) { //integral evaluation is infinity

				cout << "Error when computing numerical integral: log likelihood is infinity. Using original (compound) result instead.\n";
				error = true;
			}

			if(::isnan(logPxAphInt)!=0) { //integral evaluation is nan

				cout << "Error when computing numerical integral: log likelihood is nan. Using original (compound) result instead.\n";
				error = true;
			}
			if(error == false) { //Using integral approach value if no errors were encountered
				logPxAph = logPxAphInt;
				cout << "Using integral result.\n";
			}
			cout << "logPxAph: " << logPxAph << "\n";
		}
	}
	
	if(logPxAph == -numeric_limits<double>::max() && print == false) {
		cout << "ERROR in Misc::computeLNeutralVarAdvanced: Likelihood for xA is negative infinity. LogPxAph: " << logPxAph <<"\n";

		cout << "Recomputing mean and variance whilst printing out stats:\n";
		computeLNeutralVarAdvanced(obsA,deltaDays, Nt,ap,qBmean, qBvar, contribsPHset, true);

		cout << "-------------------------------------------------------\n";
	}

	//Likelihood is nan, so print out
	if(::isnan(logPxAph) != 0 && print == false) {
		cout << "ERROR in Misc::computeLNeutralVarAdvanced: Likelihood for xA is nan. LogPxAph: " << logPxAph <<"\n";

		cout << "Recomputing mean and variance whilst printing out stats:\n";
		computeLNeutralVarAdvanced(obsA, deltaDays, Nt,ap,qBmean, qBvar, contribsPHset, true);

		cout << "-------------------------------------------------------\n";
	}

	//Clean up before exit
	gsl_vector_free(qBmeanGSL);
	gsl_vector_free(TqBmeanGSL);

	gsl_matrix_free(T);
	gsl_matrix_free(TqBvar);
	gsl_matrix_free(TqBvarTTrans);
	gsl_matrix_free(NaAlphaPlusNaMinusAlphaGammaMTqBmean);
	gsl_matrix_free(NaNaMinusAlphaDeltaTqBvarTTrans);

	return logPxAph;
}

//Reduce dimensionality up front
//First reduce full haps by dimension k and subsequently whichever
//dimension in the partial haps that contains the full hap k
double computeLNeutralVarReduced(vector<int> &obsA, int Nt, AnaParam *ap, vector<double> &qBmean, gsl_matrix* qBvar, vector<vector<int> >* contribsPHset, bool print) {

	double C = ap->getC();
	double logPxAph = -numeric_limits<double>::max(); //Set to neg infinity as precaution

	if(print == true) {
	
		cout << "Printing intermediate steps in computeLNeutralSelGVar.\n";
		cout << "First, the input variables:\n";
		cout << "obsA: "; printIntVector(obsA);
		cout << "Nt: " << Nt << "\n";
		cout << "C: " << C << "\n";
		cout << "qBmean: "; printDoubleVector(qBmean);
		cout << "qBvar: "; printMatrixMathematica(qBvar);
		cout << "Contribs: \n";
		for(unsigned int i=0; i<contribsPHset->size(); i++) {
			cout << "PH " << i << ": "; printIntVector((*contribsPHset)[i]);
		}
	}

	//Dimensions
	int dimFH = (int) qBmean.size();
	int dimPH = (int) obsA.size();
	int dimPHminusOne = -1;
	if(dimPH==1) { dimPHminusOne = 1; } //If dimPH==1, this corresponds to at least dimPH==2 with an unobserved second ph
	else if(dimPH>1) { dimPHminusOne = dimPH-1; }
	

	//Total observations
	int Na = sumOfVector(obsA);
	
	if(print == true) {
		cout << "dimFH: " << dimFH << "\n";
		cout << "dimPH: " << dimPH << "\n";
		cout << "Na: " << Na << "\n";
	}


	//Construct matrix T (dim j-1 x k-1)
	gsl_matrix *T = gsl_matrix_alloc(dimPHminusOne,dimFH-1);

	//Find the full haplotype dimension that needs to be omitted fo reduced dimensionality purposes
	int FHdimOmitted = -1;
	if(dimPH == 1) { //Not straightforward in this case

		//WLOG, the first non-negative integer that isn't a contrib of PH0 will be a FH contrib of an unobserved PH1
		for(int i=0; i<9999; i++) {
		
			if(i<(int)(*contribsPHset)[0].size()) {
				if((*contribsPHset)[0][i] != i) { //Assumes contribs ordered. E.g. if contribs of PH0: {0,1,2,4,6} then FHdimOmitted = 3.

					FHdimOmitted = i;
					break;
				}

			} else {
				FHdimOmitted = i; //Situation where e.g. contribs is {0,1,2,3,4}, then FHdimOmitted = 5
				break;
			}
		}
		
		
	} else { //Simpler, just omit the last PH dimension (j). Then omit the first of the contribs corresponding to this PH dim
		FHdimOmitted = (*contribsPHset)[dimPH-1][0];
	}

	if(print == true) {
		cout << "FHdimOmitted: " << FHdimOmitted << "\n";
	}
	
	for(int i=0; i<dimPHminusOne; i++) { //Loop over rows of T, omit last row

		int counter=0; //Contribs counter
		for(int j=0; j<dimFH; j++) { //Loop over columns of T

			if(j<FHdimOmitted) {
				if((*contribsPHset)[i][counter]==j && counter < (int) (*contribsPHset)[i].size()) { //This assumes the contribs are ordered
					gsl_matrix_set(T,i,j,1);
					counter++;
				} else {
					gsl_matrix_set(T,i,j,0);
				}

			} else if(j>FHdimOmitted) {
				if((*contribsPHset)[i][counter]==j && counter < (int) (*contribsPHset)[i].size()) { //This assumes the contribs are ordered
					gsl_matrix_set(T,i,j-1,1);
					counter++;
				} else {
					gsl_matrix_set(T,i,j-1,0);
				}
			}
		}
	}

	if(print == true) {
		cout << "Matrix T: "; printMatrixMathematica(T);
	}

	//Create reduced quantities
	vector<double> qBmeanReduced;
	for(unsigned int i=0; i<qBmean.size(); i++) {
		if(((int) i)!=FHdimOmitted) {
			qBmeanReduced.push_back(qBmean[i]);
		}
	} 
	vector<double> obsAreduced;
	for(int i=0; i<dimPHminusOne; i++) {
		obsAreduced.push_back((double) obsA[i]);
	}
		


	/*
	* * 
	* * After computation
	* *
	*/


	/*
 	*  First we create the mean of xAPH = Na*T*qBmean
 	*/ 
	gsl_vector *qBmeanGSL = gsl_vector_alloc(dimFH-1);
	for(unsigned int i=0;i<qBmeanReduced.size();i++) {
		gsl_vector_set(qBmeanGSL,i,qBmeanReduced[i]);
	}

	//Create TqBmean through matrix vector multiplication
	gsl_vector *TqBmeanGSL = gsl_vector_alloc(dimPHminusOne);
	gsl_blas_dgemv(CblasNoTrans, 1.0, T, qBmeanGSL, 0.0, TqBmeanGSL);  

	//Convert back to C++ vector
	vector<double> meanXaPH;
	vector<double> TqBmean; //To be used in computing variance
	for(int i=0; i<dimPHminusOne; i++) {
		TqBmean.push_back(gsl_vector_get(TqBmeanGSL, i));
		meanXaPH.push_back(Na*TqBmean[i]);
	}

	if(print == true) {
		cout << "meanXaPH: "; printDoubleVector(meanXaPH);
	}

	
	/*
 	* Next we compute the variance of xAPH = Na(alpha+(Na-alpha)beta)M(TqBmean) + Na(Na-alpha)gammaTqBvarTTrans 
	*/
	//Constructing the first term 
	gsl_matrix *NaAlphaPlusNaMinusAlphaBetaMTqBmean = gsl_matrix_alloc(dimPHminusOne,dimPHminusOne);
	constructMatrixM(TqBmean, NaAlphaPlusNaMinusAlphaBetaMTqBmean);
	
	//Multiplying by the factor Na(alpha+(Na-alpha)beta):
	int growthFactor = ap->getGrowthFactor();
	int Ng=growthFactor*Nt;
        double NtNgDouble = ((double)Nt)*((double)Ng); //This computation gets out of integer range if Nt>4634, so turn into double
	double alpha = (Na+C)/((double)(1+C));
	double beta = (Nt+Ng-1)/NtNgDouble;
        double NaAlphaPlusNaMinusAlphaBeta = Na*(alpha+(Na-alpha)*beta);
	gsl_matrix_scale(NaAlphaPlusNaMinusAlphaBetaMTqBmean, NaAlphaPlusNaMinusAlphaBeta);

	if(print==true) {
                cout << "Na(alpha+(Na-alpha)beta)MTqB: "; printMatrixMathematica(NaAlphaPlusNaMinusAlphaBetaMTqBmean);
        }	
	

	//Creating TqBvarTTrans
	double gamma = (NtNgDouble-Nt-Ng+1)/NtNgDouble;
	gsl_matrix* gammaqBvar = gsl_matrix_alloc(dimFH-1,dimFH-1);
	
	//Create reduced covariance matrix
	for(int i=0; i<dimFH; i++) {
		for(int j=0; j<dimFH; j++) {
	
			int deltai = 0;
			int deltaj = 0;
			if(i>FHdimOmitted) { deltai = -1; }
			if(j>FHdimOmitted) { deltaj = -1; }

			if(i!= FHdimOmitted && j!=FHdimOmitted) {
			
				gsl_matrix_set(gammaqBvar,i+deltai,j+deltaj, gsl_matrix_get(qBvar,i,j));
			}	
		}
	}
	gsl_matrix_scale(gammaqBvar, gamma);
	if(print == true) {
		cout << "gammqBvar: "; printMatrixMathematica(gammaqBvar);
	}



	//2B)Compute T*qBvar product
	gsl_matrix* TqBvar = gsl_matrix_alloc(dimPHminusOne,dimFH-1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, gammaqBvar, 0.0, TqBvar);


	//3) Computing TqBvarT
	gsl_matrix* three = gsl_matrix_alloc(dimPHminusOne,dimPHminusOne);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, TqBvar, T, 0.0, three);

	if(print == true) {
		cout << "3: "; printMatrixMathematica(three);
	}

	//
	gsl_matrix_scale(three, Na*(Na-alpha));

	
	//Add together
	gsl_matrix_add(NaAlphaPlusNaMinusAlphaBetaMTqBmean, three); //Final variance for after case. Output stored in first matrix.

	if(print==true) {

		cout << "Final quantities:\n";
		cout << "Reduced mean: "; printDoubleVector(meanXaPH);
		cout << "Reduced var: "; printMatrixMathematica(NaAlphaPlusNaMinusAlphaBetaMTqBmean);
		cout << "Reduced obsA: "; printDoubleVector(obsAreduced);
		cout << "Na: " << Na << "\n";

	}

	
	//Computing likelihoods depend on dimensionality
	if(dimPH > 2) { //Multivariate approach	

		//Compute likelihood
		logPxAph = logMultivariateNormalPDF(obsAreduced,meanXaPH,NaAlphaPlusNaMinusAlphaBetaMTqBmean);
		if(print==true) {
			logMultivariateNormalPDFPrint(obsAreduced,meanXaPH,NaAlphaPlusNaMinusAlphaBetaMTqBmean);
		}
		
		if(ap->getUseIntegralApproach() == true) {
			double* xmin = (double *) malloc(dimPHminusOne * sizeof(double));
			double* xmax = (double *) malloc(dimPHminusOne * sizeof(double));
			double* err = (double *) malloc(intOutputDim * sizeof(double));
			err[0] = 9999;
			for(int i=0; i<dimPHminusOne; i++) {
		
				xmin[i] = ((double) obsAreduced[i]) - 0.5; //Define hypercube around observation
				xmax[i] = ((double) obsAreduced[i]) + 0.5; //Boundaries may be < 0 or > Na
			}

			parameters params;
			params.muVec = meanXaPH;
			params.Sigma = NaAlphaPlusNaMinusAlphaBetaMTqBmean;
			double PxAph = 0; //Not log parameter
		
			//Perform integral. For dim <= 3, pcubature is best. hcubature is better for larger dimensions, although best for <7.
			int intOutcome = -1;
			if(dimPHminusOne <= 3) {
				intOutcome = pcubature(intOutputDim, integrandMultiNormal, &params, dimPHminusOne, xmin, xmax, maxEval, tol, tol, ERROR_INDIVIDUAL, &PxAph, err); 
			} else {
				intOutcome = hcubature(intOutputDim, integrandMultiNormal, &params, dimPHminusOne, xmin, xmax, maxEval, tol, tol, ERROR_INDIVIDUAL, &PxAph, err); 
			}
			bool error = false;
			if(intOutcome != 0) { //Error in integral

				cout << "Error computing numerical integral. Using original (compound) result instead.\n";
				error = true;
			}
			if(err[0] > 0.001) {
				
				cout << "Error in numerical integral is too large. The error is: " << err[0] << ". Using original (compound) result instead.\n";
				error = true;
			}
			double logPxAphInt = log(PxAph);
			if(::isinf(logPxAphInt) != 0) { //integral evaluation is infinity

				cout << "Error when computing numerical integral: log likelihood is infinity. Using original (compound) result instead.\n";
				error = true;
			}

			if(::isnan(logPxAphInt)!=0) { //integral evaluation is nan

				cout << "Error when computing numerical integral: log likelihood is nan. Using original (compound) result instead.\n";
				error = true;
			}
			if(error == false) { //Using integral approach value if no errors were encountered
				logPxAph = logPxAphInt;
				cout << "Using integral result.\n";
			}
			cout << "logPxAph: " << logPxAph << "\n";

		}

	}else { //I.e. dimPH==2 or dimPH == 1 (same computation)
		
		//Compute multivariate normal likelihood
		//In this case we retrieve the univariate case
		double meanXaPHuni = meanXaPH[0];
		double varXaPHuni = gsl_matrix_get(NaAlphaPlusNaMinusAlphaBetaMTqBmean,0,0);
		double stDevXaPHuni = sqrt(varXaPHuni);
		logPxAph = -(obsA[0]-meanXaPHuni)*(obsA[0]-meanXaPHuni)/(2*varXaPHuni) - log(stDevXaPHuni*sqrt(2*M_PI)); //log L
		

		if(ap->getUseIntegralApproach() == true) {

			double* xmin = (double *) malloc(1 * sizeof(double));
			double* xmax = (double *) malloc(1 * sizeof(double));
			xmin[0] = ((double) obsAreduced[0]) - 0.5; //Define +/- 0.5 interval around observation
			xmax[0] = ((double) obsAreduced[0]) + 0.5; //Boundaries may be < 0 or > Na
			double* err = (double *) malloc(intOutputDim * sizeof(double));
			err[0] = 9999;

			parameters params;
			params.mu = meanXaPHuni;
			params.var = varXaPHuni;
			params.stDev = stDevXaPHuni;
			double PxAph = 0; //Not log parameter
		
			//Perform integral using pcubature (best for dim <=3)
			int intOutcome = pcubature(intOutputDim, integrandNormal, &params, 1, xmin, xmax, maxEval, tol, tol, ERROR_INDIVIDUAL, &PxAph, err); 
			bool error = false;
			if(intOutcome != 0) { //Error in integral

				cout << "Error computing numerical integral. Using original (compound) result instead.\n";
				error = true;
			}
			if(err[0] > 0.001) {
				
				cout << "Error in numerical integral is too large. The error is: " << err[0] << ". Using original (compound) result instead.\n";
				error = true;
			}
			double logPxAphInt = log(PxAph);
			if(::isinf(logPxAphInt) != 0) { //integral evaluation is infinity

				cout << "Error when computing numerical integral: log likelihood is infinity. Using original (compound) result instead.\n";
				error = true;
			}

			if(::isnan(logPxAphInt)!=0) { //integral evaluation is nan

				cout << "Error when computing numerical integral: log likelihood is nan. Using original (compound) result instead.\n";
				error = true;
			}
			if(error == false) { //Using integral approach value if no errors were encountered
				logPxAph = logPxAphInt;
				cout << "Using integral result.\n";
			}
			cout << "logPxAph: " << logPxAph << "\n";

		}
	}

	//Check of likelihoods are still negative infinity, which should not happen
	if(logPxAph == -numeric_limits<double>::max() && print == false) {
		cout << "ERROR in Misc::computeLNeutralVarReduced: Likelihood for xA is negative infinity. LogPxAph: " << logPxAph <<"\n";

		cout << "Recomputing mean and variance whilst printing out stats:\n";
		computeLNeutralVarReduced(obsA,Nt,ap,qBmean, qBvar, contribsPHset, true);

		cout << "-------------------------------------------------------\n";
	}

	//Likelihood is nan, so print out
	if(::isnan(logPxAph) != 0 && print == false) {
		cout << "ERROR in Misc::computeLNeutralVarReduced: Likelihood for xA is nan. LogPxAph: " << logPxAph <<"\n";

		cout << "Recomputing mean and variance whilst printing out stats:\n";
		computeLNeutralVarReduced(obsA,Nt,ap,qBmean, qBvar, contribsPHset, true);

		cout << "-------------------------------------------------------\n";
	}
	
	//Clean up before exit
	gsl_vector_free(qBmeanGSL);
	gsl_vector_free(TqBmeanGSL);

	gsl_matrix_free(T);
	gsl_matrix_free(gammaqBvar);
	gsl_matrix_free(three);


	if(print == true) {
		cout << "Returning likelihood " << logPxAph << "\n"; 
	}

	return logPxAph;
}




//Reduce dimensionality up front
//First reduce full haps by dimension k and subsequently whichever
//dimension in the partial haps that contains the full hap k
double computeLNeutralSelGVar(vector<int> &obsA, int Nt, AnaParam *ap, vector<double> &qBmean, gsl_matrix* qBvar, vector<double>& hapFitG, vector<vector<int> >* contribsPHset, bool print) {

	double C = ap->getC();
	double logPxAph = -numeric_limits<double>::max(); //Set to neg infinity as precaution

	if(print == true) {
	
		cout << "Printing intermediate steps in computeLNeutralSelGVar.\n";
		cout << "First, the input variables:\n";
		cout << "obsA: "; printIntVector(obsA);
		cout << "Nt: " << Nt << "\n";
		cout << "C: " << C << "\n";
		cout << "qBmean: "; printDoubleVector(qBmean);
		cout << "qBvar: "; printMatrixMathematica(qBvar);
		cout << "hapFitG: "; printDoubleVector(hapFitG);
		cout << "Contribs: \n";
		for(unsigned int i=0; i<contribsPHset->size(); i++) {
			cout << "PH " << i << ": "; printIntVector((*contribsPHset)[i]);
		}
	}

	//Dimensions
	int dimFH = (int) qBmean.size();
	int dimPH = (int) obsA.size();
	int dimPHminusOne = -1;
	if(dimPH==1) { dimPHminusOne = 1; } //If dimPH==1, this corresponds to at least dimPH==2 with an unobserved second ph
	else if(dimPH>1) { dimPHminusOne = dimPH-1; }
	

	//Total observations
	int Na = sumOfVector(obsA);
	
	if(print == true) {
		cout << "dimFH: " << dimFH << "\n";
		cout << "dimPH: " << dimPH << "\n";
		cout << "Na: " << Na << "\n";
	}


	//Construct matrix T (dim j-1 x k-1)
	gsl_matrix *T = gsl_matrix_alloc(dimPHminusOne,dimFH-1);

	//Find the full haplotype dimension that needs to be omitted fo reduced dimensionality purposes
	int FHdimOmitted = -1;
	if(dimPH == 1) { //Not straightforward in this case

		//WLOG, the first non-negative integer that isn't a contrib of PH0 will be a FH contrib of an unobserved PH1
		for(int i=0; i<9999; i++) {
		
			if(i<(int)(*contribsPHset)[0].size()) {
				if((*contribsPHset)[0][i] != i) { //Assumes contribs ordered. E.g. if contribs of PH0: {0,1,2,4,6} then FHdimOmitted = 3.

					FHdimOmitted = i;
					break;
				}

			} else {
				FHdimOmitted = i; //Situation where e.g. contribs is {0,1,2,3,4}, then FHdimOmitted = 5
				break;
			}
		}
		
		
	} else { //Simpler, just omit the last PH dimension (j). Then omit the first of the contribs corresponding to this PH dim
		FHdimOmitted = (*contribsPHset)[dimPH-1][0];
	}

	if(print == true) {
		cout << "FHdimOmitted: " << FHdimOmitted << "\n";
	}
	
	for(int i=0; i<dimPHminusOne; i++) { //Loop over rows of T, omit last row

		int counter=0; //Contribs counter
		for(int j=0; j<dimFH; j++) { //Loop over columns of T

			if(j<FHdimOmitted) {
				if((*contribsPHset)[i][counter]==j && counter < (int) (*contribsPHset)[i].size()) { //This assumes the contribs are ordered
					gsl_matrix_set(T,i,j,1);
					counter++;
				} else {
					gsl_matrix_set(T,i,j,0);
				}

			} else if(j>FHdimOmitted) {
				if((*contribsPHset)[i][counter]==j && counter < (int) (*contribsPHset)[i].size()) { //This assumes the contribs are ordered
					gsl_matrix_set(T,i,j-1,1);
					counter++;
				} else {
					gsl_matrix_set(T,i,j-1,0);
				}
			}
		}
	}

	if(print == true) {
		cout << "Matrix T: "; printMatrixMathematica(T);
	}

	//Create reduced quantities
	vector<double> qBmeanReduced;
	for(unsigned int i=0; i<qBmean.size(); i++) {
		if(((int) i)!=FHdimOmitted) {
			qBmeanReduced.push_back(qBmean[i]);
		}
	} 
	vector<double> obsAreduced;
	for(int i=0; i<dimPHminusOne; i++) {
		obsAreduced.push_back((double) obsA[i]);
	}
		


	/*
	* * 
	* * After computation
	* *
	*/


	/*
 	* First we create the mean of xAPH = Na*T*sigmaG[qBmean]
 	*/ 
	vector<double> sigmaGqBmean;
	double totalG = 0;
	for(int i=0; i<dimFH; i++) {
		totalG += qBmean[i]*exp(hapFitG[i]);
	}
	if(print == true) {
		cout << "totalG: " << totalG << "\n";
	}
	for(int i=0; i<dimFH; i++) {
		if(i!=FHdimOmitted) {
			sigmaGqBmean.push_back(qBmean[i]*exp(hapFitG[i])/totalG);
//			sigmaGqBmean.push_back(qBmean[i]);
		}
	}

	//Rescale sigmaGqBmean to avoid issues with frequencies < 10^-10
	rescale(sigmaGqBmean); 

	if(print == true) {
		cout << "sigmaGqBmean: "; printDoubleVector(sigmaGqBmean);
	}

	gsl_vector *sigmaGqBmeanGSL = gsl_vector_alloc(dimFH-1);
	for(unsigned int i=0;i<sigmaGqBmean.size();i++) {
		gsl_vector_set(sigmaGqBmeanGSL,i,sigmaGqBmean[i]);
	}

	//Create TsigmaGqBmean through matrix vector multiplication
	gsl_vector *TsigmaGqBmeanGSL = gsl_vector_alloc(dimPHminusOne);
	gsl_blas_dgemv(CblasNoTrans, 1.0, T, sigmaGqBmeanGSL, 0.0, TsigmaGqBmeanGSL);  

	//Convert back to C++ vector
	vector<double> meanXaPH;
	vector<double> TsigmaGqBmean; //To be used in computing variance
	for(int i=0; i<dimPHminusOne; i++) {
		meanXaPH.push_back(Na*(gsl_vector_get(TsigmaGqBmeanGSL, i)));
		TsigmaGqBmean.push_back(gsl_vector_get(TsigmaGqBmeanGSL, i));
	}

	if(print == true) {
		cout << "meanXaPH: "; printDoubleVector(meanXaPH);
		cout << "TsigmaGqBmean: "; printDoubleVector(TsigmaGqBmean);
	}

	
	/*
	* Next we compute the variance of xAPH = alpha*Na*M(TsigmaG[qBmean]) + Na*(Na-alpha)*T*D[sigmaG]_qBmean * (beta*M(qBmean) + gamma* qBvar) * (D[sigmaG]_qBmean)^T * T^T
	*/
	//Constructing the first term alpha*Na*M(TsigmaG[qBmean])
	gsl_matrix *alphaNaMTsigmaGqBmean = gsl_matrix_alloc(dimPHminusOne,dimPHminusOne);
	constructMatrixM(TsigmaGqBmean, alphaNaMTsigmaGqBmean);
	if(print==true) {
		cout << "MTsigmaGqBmean: "; printMatrixMathematica(alphaNaMTsigmaGqBmean);
	}
	

	//Multiplying by the factor alpha*Na:
	double alpha = (Na+C)/((double)(1+C));
	gsl_matrix_scale(alphaNaMTsigmaGqBmean, alpha*Na);

	if(print==true) {
		cout << "alpha: " << alpha << "\n";
		cout << "alpha*Na*MTsigmaGqBmean: "; printMatrixMathematica(alphaNaMTsigmaGqBmean);
	}
	
	

	//For constructing the second term Na*(Na-alpha)*T*D[sigmaG]_qBmean * (beta*M(qBmean) + gamma * qBvar) * (D[sigmaG]_qBmean)^T * T^T
	//we break it up into 4 parts:
	//1) Compute beta*M(qBmean) + gamma* qBvar
	//2) Compute T*D[sigmaG]_qBmean
	//3) Compute 1*2^T
	//4) Compute Na(Na-alpha)*2*3

	//Doing 1) first
	//1a) Compute betaMqBmean
	gsl_matrix *betaMqBmean = gsl_matrix_alloc(dimFH-1, dimFH-1);
	constructMatrixM(qBmeanReduced, betaMqBmean);
	int growthFactor = ap->getGrowthFactor();
	int Ng=growthFactor*Nt;
	double NtNgDouble = ((double)Nt)*((double)Ng); //This computation gets out of integer range if Nt>4634, so turn into double
	double beta = (Nt+Ng-1)/NtNgDouble;
	gsl_matrix_scale(betaMqBmean, beta);

	if(print == true) {
		cout << "betaMqBmean: "; printMatrixMathematica(betaMqBmean);
	}

	//1b) Compute gamma*qBvar
//	double gamma = (NtNgDouble-Nt-Ng+1)/NtNgDouble;
//	gsl_matrix_view qBvarReducedView = gsl_matrix_submatrix(qBvar, 0,0, dimFH-1, dimFH-1); //Remove row k and column k
//	gsl_matrix* gammaqBvar = gsl_matrix_alloc(dimFH-1,dimFH-1);
//	gsl_matrix_memcpy(gammaqBvar, &(qBvarReducedView.matrix)); //Dest,source
//	gsl_matrix_scale(gammaqBvar, gamma);
	double gamma = (NtNgDouble-Nt-Ng+1)/NtNgDouble;
	gsl_matrix* gammaqBvar = gsl_matrix_alloc(dimFH-1,dimFH-1);
	
	//Create reduced covariance matrix
	for(int i=0; i<dimFH; i++) {
		for(int j=0; j<dimFH; j++) {
	
			int deltai = 0;
			int deltaj = 0;
			if(i>FHdimOmitted) { deltai = -1; }
			if(j>FHdimOmitted) { deltaj = -1; }

			if(i!= FHdimOmitted && j!=FHdimOmitted) {
			
				gsl_matrix_set(gammaqBvar,i+deltai,j+deltaj, gsl_matrix_get(qBvar,i,j));
			}	
		}
	}
	gsl_matrix_scale(gammaqBvar, gamma);
	if(print == true) {
		cout << "gammqBvar: "; printMatrixMathematica(gammaqBvar);
	}

	//If qBvar is used, add together
	if(ap->getNoVar() == false && ap->getMeanOnly() == false) {

	
		//1c) Perform 1a+1b (resulting variable is 1a)	
		gsl_matrix_add(betaMqBmean, gammaqBvar);
		if(print == true) {
			cout << "1c: "; printMatrixMathematica(betaMqBmean);
		}
	} else {

		//Don't add together, i.e. betaMqBmean remains unchanged
	}



	//Computing term 2: T*D[sigmaG]_qBmean
	//2A) Computing D[sigmaG]_qBmean
	gsl_matrix* DsigmaG_qBmean = gsl_matrix_alloc(dimFH-1,dimFH-1);
	for(int i=0; i<dimFH; i++) {
		for(int j=0; j<dimFH; j++) {
	
			int deltai=0;
			int deltaj=0;
			if(i>FHdimOmitted) { deltai = -1; }
			if(j>FHdimOmitted) { deltaj = -1; }

			if(i!=FHdimOmitted && j!=FHdimOmitted) {
				if(i==j) {
					gsl_matrix_set(DsigmaG_qBmean, i+deltai, j+deltaj, exp(hapFitG[i])/totalG -qBmean[i]*exp(hapFitG[i]) * (exp(hapFitG[i]) -exp(hapFitG[FHdimOmitted]))/(totalG*totalG));
//cout << "value for " << i+deltai << "," << j+deltaj << ": " << gsl_matrix_get(DsigmaG_qBmean,i+deltai, j+deltaj) << "\n"; 
//gsl_matrix_set(DsigmaG_qBmean,i+deltai,j+deltaj,1);
				} else {
				gsl_matrix_set(DsigmaG_qBmean, i+deltai, j+deltaj, -qBmean[i]*exp(hapFitG[i]) * (exp(hapFitG[j]) - exp(hapFitG[FHdimOmitted]))/(totalG*totalG));
//gsl_matrix_set(DsigmaG_qBmean,i+deltai,j+deltaj,0);
				}
			}
		}
	}
	
	if(print == true) {
		cout << "DsigmaG_qBmean: "; printMatrixMathematica(DsigmaG_qBmean);
	}

	//2B)Compute T*D product
	gsl_matrix* TDsigmaG_qBmean = gsl_matrix_alloc(dimPHminusOne,dimFH-1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, DsigmaG_qBmean, 0.0, TDsigmaG_qBmean);


	//3) Computing term 3 = 1*2^T
	gsl_matrix* three = gsl_matrix_alloc(dimFH-1,dimPHminusOne);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, betaMqBmean, TDsigmaG_qBmean, 0.0, three);

	if(print == true) {
		cout << "3: "; printMatrixMathematica(three);
	}
	//4) Computing term four = Na(Na-alpha)*2*3
	gsl_matrix* four = gsl_matrix_alloc(dimPHminusOne,dimPHminusOne);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, Na*(Na-alpha), TDsigmaG_qBmean, three, 0.0, four);

	if(print == true) {
		cout << "4: "; printMatrixMathematica(four);
	}

	if(print==true) {
		cout << "Na*(Na-alpha)*T*D[sigmaG]_qBmean * (beta*M(qBmean) + gamma * qBvar) * (D[sigmaG]_qBmean)^T * T^T:  "; printMatrixMathematica(four);
	}
	
	//Add together
	gsl_matrix_add(alphaNaMTsigmaGqBmean, four); //Final variance for after case. Output stored in first matrix.

	if(print==true) {

		cout << "Final quantities:\n";
		cout << "Reduced mean: "; printDoubleVector(meanXaPH);
		cout << "Reduced var: "; printMatrixMathematica(alphaNaMTsigmaGqBmean);
		cout << "Reduced obsA: "; printDoubleVector(obsAreduced);
		cout << "Na: " << Na << "\n";

	}

	
	//Computing likelihoods depend on dimensionality
	if(dimPH > 2) { //Multivariate approach	

		//Compute likelihood
		logPxAph = logMultivariateNormalPDF(obsAreduced,meanXaPH,alphaNaMTsigmaGqBmean);
		if(print==true) {
			logMultivariateNormalPDFPrint(obsAreduced,meanXaPH,alphaNaMTsigmaGqBmean);
		}
		
		if(ap->getUseIntegralApproach() == true) {
			double* xmin = (double *) malloc(dimPHminusOne * sizeof(double));
			double* xmax = (double *) malloc(dimPHminusOne * sizeof(double));
			double* err = (double *) malloc(intOutputDim * sizeof(double));
			err[0] = 9999;
			for(int i=0; i<dimPHminusOne; i++) {
		
				xmin[i] = ((double) obsAreduced[i]) - 0.5; //Define hypercube around observation
				xmax[i] = ((double) obsAreduced[i]) + 0.5; //Boundaries may be < 0 or > Na
			}

			parameters params;
			params.muVec = meanXaPH;
			params.Sigma = alphaNaMTsigmaGqBmean;
			double PxAph = 0; //Not log parameter
		
			//Perform integral. For dim <= 3, pcubature is best. hcubature is better for larger dimensions, although best for <7.
			int intOutcome = -1;
			if(dimPHminusOne <= 3) {
				intOutcome = pcubature(intOutputDim, integrandMultiNormal, &params, dimPHminusOne, xmin, xmax, maxEval, tol, tol, ERROR_INDIVIDUAL, &PxAph, err); 
			} else {
				intOutcome = hcubature(intOutputDim, integrandMultiNormal, &params, dimPHminusOne, xmin, xmax, maxEval, tol, tol, ERROR_INDIVIDUAL, &PxAph, err); 
			}
			bool error = false;
			if(intOutcome != 0) { //Error in integral

				cout << "Error computing numerical integral. Using original (compound) result instead.\n";
				error = true;
			}
			if(err[0] > 0.001) {
				
				cout << "Error in numerical integral is too large. The error is: " << err[0] << ". Using original (compound) result instead.\n";
				error = true;
			}
			double logPxAphInt = log(PxAph);
			if(::isinf(logPxAphInt) != 0) { //integral evaluation is infinity

				cout << "Error when computing numerical integral: log likelihood is infinity. Using original (compound) result instead.\n";
				error = true;
			}

			if(::isnan(logPxAphInt)!=0) { //integral evaluation is nan

				cout << "Error when computing numerical integral: log likelihood is nan. Using original (compound) result instead.\n";
				error = true;
			}
			if(error == false) { //Using integral approach value if no errors were encountered
				logPxAph = logPxAphInt;
				cout << "Using integral result.\n";
			}
			cout << "logPxAph: " << logPxAph << "\n";

		}

	}else { //I.e. dimPH==2 or dimPH == 1 (same computation)
		
		//Compute multivariate normal likelihood
		//In this case we retrieve the univariate case
		double meanXaPHuni = meanXaPH[0];
		double varXaPHuni = gsl_matrix_get(alphaNaMTsigmaGqBmean,0,0);
		double stDevXaPHuni = sqrt(varXaPHuni);
		logPxAph = -(obsA[0]-meanXaPHuni)*(obsA[0]-meanXaPHuni)/(2*varXaPHuni) - log(stDevXaPHuni*sqrt(2*M_PI)); //log L
		

		if(ap->getUseIntegralApproach() == true) {

			double* xmin = (double *) malloc(1 * sizeof(double));
			double* xmax = (double *) malloc(1 * sizeof(double));
			xmin[0] = ((double) obsAreduced[0]) - 0.5; //Define +/- 0.5 interval around observation
			xmax[0] = ((double) obsAreduced[0]) + 0.5; //Boundaries may be < 0 or > Na
			double* err = (double *) malloc(intOutputDim * sizeof(double));
			err[0] = 9999;

			parameters params;
			params.mu = meanXaPHuni;
			params.var = varXaPHuni;
			params.stDev = stDevXaPHuni;
			double PxAph = 0; //Not log parameter
		
			//Perform integral using pcubature (best for dim <=3)
			int intOutcome = pcubature(intOutputDim, integrandNormal, &params, 1, xmin, xmax, maxEval, tol, tol, ERROR_INDIVIDUAL, &PxAph, err); 
			bool error = false;
			if(intOutcome != 0) { //Error in integral

				cout << "Error computing numerical integral. Using original (compound) result instead.\n";
				error = true;
			}
			if(err[0] > 0.001) {
				
				cout << "Error in numerical integral is too large. The error is: " << err[0] << ". Using original (compound) result instead.\n";
				error = true;
			}
			double logPxAphInt = log(PxAph);
			if(::isinf(logPxAphInt) != 0) { //integral evaluation is infinity

				cout << "Error when computing numerical integral: log likelihood is infinity. Using original (compound) result instead.\n";
				error = true;
			}

			if(::isnan(logPxAphInt)!=0) { //integral evaluation is nan

				cout << "Error when computing numerical integral: log likelihood is nan. Using original (compound) result instead.\n";
				error = true;
			}
			if(error == false) { //Using integral approach value if no errors were encountered
				logPxAph = logPxAphInt;
				cout << "Using integral result.\n";
			}
			cout << "logPxAph: " << logPxAph << "\n";

		}
	}

	//Check of likelihoods are still negative infinity, which should not happen
	if(logPxAph == -numeric_limits<double>::max() && print == false) {
		cout << "ERROR in Misc::computeLNeutralSelGVar: Likelihood for xA is negative infinity. LogPxAph: " << logPxAph <<"\n";

		cout << "Recomputing mean and variance whilst printing out stats:\n";
		computeLNeutralSelGVar(obsA,Nt,ap,qBmean, qBvar, hapFitG, contribsPHset, true);

		cout << "-------------------------------------------------------\n";
	}

	//Likelihood is nan, so print out
	if(::isnan(logPxAph) != 0 && print == false) {
		cout << "ERROR in Misc::computeLNeutralSelGVar: Likelihood for xA is nan. LogPxAph: " << logPxAph <<"\n";

		cout << "Recomputing mean and variance whilst printing out stats:\n";
		computeLNeutralSelGVar(obsA,Nt,ap,qBmean, qBvar, hapFitG, contribsPHset, true);

		cout << "-------------------------------------------------------\n";
	}
	
	//Clean up before exit
	gsl_vector_free(sigmaGqBmeanGSL);
	gsl_vector_free(TsigmaGqBmeanGSL);

	gsl_matrix_free(T);
	gsl_matrix_free(alphaNaMTsigmaGqBmean);
	gsl_matrix_free(betaMqBmean);
	gsl_matrix_free(gammaqBvar);
	gsl_matrix_free(DsigmaG_qBmean);
	gsl_matrix_free(TDsigmaG_qBmean);
	gsl_matrix_free(three);
	gsl_matrix_free(four);


	if(print == true) {
		cout << "Returning likelihood " << logPxAph << "\n"; 
	}

	return logPxAph;
}




double computeLNeutralSelGVarOld(vector<int> &obsA, int Nt, AnaParam *ap, vector<double> &qBmean, gsl_matrix* qBvar, vector<double>& hapFitG, vector<vector<int> >* contribsPHset, bool print) {

	double C = ap->getC();

	if(print == true) {
	
		cout << "Printing intermediate steps in computeLselection.\n";
		cout << "First, the input variables:\n";
		cout << "obsA: "; printIntVector(obsA);
		cout << "Nt: " << Nt << "\n";
		cout << "C: " << C << "\n";
		cout << "qBmean: "; printDoubleVector(qBmean);
		cout << "Contribs: \n";
		for(unsigned int i=0; i<contribsPHset->size(); i++) {
			cout << "PH " << i << ": "; printIntVector((*contribsPHset)[i]);
		}


	}

	//Dimensions
	int dimFH = (int) qBmean.size();
	int dimPH = (int) obsA.size();
	

	//Total observations
	int Na = sumOfVector(obsA);
	
	if(print == true) {
		cout << "dimFH: " << dimFH << "\n";
		cout << "dimPH: " << dimPH << "\n";
		cout << "Na: " << Na << "\n";
	}

	//Construct matrix T
	gsl_matrix *T = gsl_matrix_alloc(dimPH,dimFH);
	for(int i=0; i<dimPH; i++) { //Loop over rows

		int counter=0;
		for(int j=0; j<dimFH; j++) { //Loop over columns

			if((*contribsPHset)[i][counter]==j && counter < (int) (*contribsPHset)[i].size()) { //This assumes the contribs are ordered
				gsl_matrix_set(T,i,j,1);
				counter++;
			} else {
				gsl_matrix_set(T,i,j,0);
			}
		}
	}

	if(print == true) {
		cout << "Matrix T: "; printMatrixMathematica(T);
	}


	/*
	* * 
	* * After computation
	* *
	*/


	/*
 	* First we create the mean of xAPH = Na*T*sigmaG[qBmean]
 	*/ 
	vector<double> sigmaGqBmean;
	double totalG = 0;
	for(int i=0; i<dimFH; i++) {
		totalG += qBmean[i]*exp(hapFitG[i]);
	}
	if(print == true) {
		cout << "totalG: " << totalG << "\n";
	}
	for(int i=0; i<dimFH; i++) {
		sigmaGqBmean.push_back(qBmean[i]*exp(hapFitG[i])/totalG);
	}

	if(print == true) {
		cout << "sigmaGqBmean: "; printDoubleVector(sigmaGqBmean);
	}

	gsl_vector *sigmaGqBmeanGSL = gsl_vector_alloc(dimFH);
	for(int i=0;i<dimFH;i++) {
		gsl_vector_set(sigmaGqBmeanGSL,i,sigmaGqBmean[i]);
	}

	//Create TsigmaGsigmaTqBmean through matrix vector multiplication
	gsl_vector *TsigmaGqBmeanGSL = gsl_vector_alloc(dimPH);
	gsl_blas_dgemv(CblasNoTrans, 1.0, T, sigmaGqBmeanGSL, 0.0, TsigmaGqBmeanGSL);  

	//Convert back to C++ vector
	vector<double> meanXaPH;
	vector<double> TsigmaGqBmean; //To be used in computing variance
	for(int i=0; i<dimPH; i++) {
		meanXaPH.push_back(Na*(gsl_vector_get(TsigmaGqBmeanGSL, i)));
		TsigmaGqBmean.push_back(gsl_vector_get(TsigmaGqBmeanGSL, i));
	}

	if(print == true) {
		cout << "meanXaPH: "; printDoubleVector(meanXaPH);
	}


	
	/*
	* Next we compute the variance of xAPH = alpha*Na*M(TsigmaG[qBmean]) + Na*(Na-alpha)*T*D[sigmaG]_qBmean * (beta*M(qBmean) + gamma* qBvar) * (D[sigmaG]_qBmean)^T * T^T
	*/
	//Constructing the first term alpha*Na*M(TsigmaG[qBmean])
	gsl_matrix *alphaNaMTsigmaGqBmean = gsl_matrix_alloc(dimPH,dimPH);
	constructMatrixM(TsigmaGqBmean, alphaNaMTsigmaGqBmean);
	
	//Multiplying by the factor alpha*Na:
	double alpha = (Na+C)/((double)(1+C));
	gsl_matrix_scale(alphaNaMTsigmaGqBmean, alpha*Na);

	if(print==true) {
		cout << "alpha*Na*MTsigmaGqBmean: "; printMatrixMathematica(alphaNaMTsigmaGqBmean);
	}
	


	//For constructing the second term Na*(Na-alpha)*T*D[sigmaG]_qBmean * (beta*M(qBmean) + gamma * qBvar) * (D[sigmaG]_qBmean)^T * T^T
	//we break it up into 4 parts:
	//1) Compute beta*M(qBmean) + gamma* qBvar
	//2) Compute T*D[sigmaG]_qBmean
	//3) Compute 1*2^T
	//4) Compute Na(Na-alpha)*2*3

	//Doing 1) first
	//1a) Compute betaMqBmean
	gsl_matrix *betaMqBmean = gsl_matrix_alloc(dimFH, dimFH);
	constructMatrixM(qBmean, betaMqBmean);
	int growthFactor = ap->getGrowthFactor();
	int Ng=growthFactor*Nt;
	double NtNgDouble = ((double)Nt)*((double)Ng); //This computation gets out of integer range if Nt>4634, so turn into double
	double beta = (Nt+Ng-1)/NtNgDouble;
	gsl_matrix_scale(betaMqBmean, beta);

	if(print == true) {
		cout << "betaMqBmean: "; printMatrixMathematica(betaMqBmean);
	}

	//1b) Compute gamma*qBvar
	double gamma = (NtNgDouble-Nt-Ng+1)/NtNgDouble;
	gsl_matrix* gammaqBvar = gsl_matrix_alloc(dimFH,dimFH);
	gsl_matrix_memcpy(gammaqBvar, qBvar); //Dest,source
	gsl_matrix_scale(gammaqBvar, gamma);

	//1c) Perform 1a+1b (resulting variable is 1a)	
	gsl_matrix_add(betaMqBmean, gammaqBvar);

	//Computing term 2: T*D[sigmaG]_qBmean
	//2A) Computing D[sigmaG]_qBmean
	gsl_matrix* DsigmaG_qBmean = gsl_matrix_alloc(dimFH,dimFH);
	for(int i=0; i<dimFH; i++) {
		for(int j=0; j<dimFH; j++) {
	
			if(i==j) {
				gsl_matrix_set(DsigmaG_qBmean, i, j, exp(hapFitG[i])/totalG - qBmean[i]*exp(2*hapFitG[i])/(totalG*totalG));
			} else {
				gsl_matrix_set(DsigmaG_qBmean, i, j, -qBmean[i]*exp(hapFitG[i]+hapFitG[j])/(totalG*totalG));
			}
		}
	}	

	//2B)Compute T*D product
	gsl_matrix* TDsigmaG_qBmean = gsl_matrix_alloc(dimPH,dimFH);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, DsigmaG_qBmean, 0.0, TDsigmaG_qBmean);

	//3) Computing term 3 = 1*2^T
	gsl_matrix* three = gsl_matrix_alloc(dimFH,dimPH);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, betaMqBmean, TDsigmaG_qBmean, 0.0, three);

	//4) Computing term four = Na(Na-alpha)*2*3
	gsl_matrix* four = gsl_matrix_alloc(dimPH,dimPH);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, Na*(Na-alpha), TDsigmaG_qBmean, three, 0.0, four);


	if(print==true) {
		cout << "Na*(Na-alpha)*T*D[sigmaG]_qBmean * (beta*M(qBmean) + gamma * qBvar) * (D[sigmaG]_qBmean)^T * T^T:  "; printMatrixMathematica(four);
	}
	
	//Add together
	gsl_matrix_add(alphaNaMTsigmaGqBmean, four); //Final variance for after case. Output stored in first matrix.

	
	//Computing likelihoods depend on dimensionality
	double logPxAph = -numeric_limits<double>::max(); //Set to neg infinity as precaution

	if(dimPH > 1) { //Most common case, i.e. multiple partial haps in a set

		/*
		 * Reduce dimensionality by 1 to ensure non-degeneracy
		 */
		//After case
		meanXaPH.pop_back(); //WLOG remove the last PH
		gsl_matrix_view varXaPHreducedView = gsl_matrix_submatrix(alphaNaMTsigmaGqBmean, 0,0, dimPH-1, dimPH-1); //Remove row k and column k
		gsl_matrix * varXaPH = gsl_matrix_alloc(dimPH-1, dimPH-1);
		gsl_matrix_memcpy(varXaPH, &(varXaPHreducedView.matrix)); //Copy the reduced matrix view into a "normal" matrix varXaPH
		vector<double> obsAreduced;
		for(int i=0; i<dimPH-1; i++) {
			obsAreduced.push_back(obsA[i]);
		}

		if(print == true) {
			cout << "Reduced quantities: \n";
			cout << "ObsAreduced: "; printDoubleVector(obsAreduced);
			cout << "meanXaPHreduced: "; printDoubleVector(meanXaPH);
			cout << "varXaPHreduced: "; printMatrixMathematica(varXaPH);
		}

		if(dimPH > 2) { //Multivariate approach	

			//Compute likelihood
			logPxAph = logMultivariateNormalPDF(obsAreduced,meanXaPH,varXaPH);
		
        
		}else { //I.e. dimPH==2
	
			//In this case we retrieve the univariate case when reducing dimensionality by one
			double meanXaPHuni = meanXaPH[0];
			double varXaPHuni = gsl_matrix_get(varXaPH,0,0);
			double stDevXaPHuni = sqrt(varXaPHuni);

			logPxAph = -(obsA[0]-meanXaPHuni)*(obsA[0]-meanXaPHuni)/(2*varXaPHuni) - log(stDevXaPHuni*sqrt(2*M_PI)); //log L
		}

		//Deallocate variables created in this scope
		gsl_matrix_free(varXaPH);

	//In theory, can combine dimPH==2 from above with this one, would be slightly faster, although more clear this way
	} else if (dimPH==1) {
		//Case where there is a single partial haplotype in a set.
		//This can be thought of as a case with a single observed partial haplotype
		//and one (or several) unobserved haplotypes.
		//If all loci in the partial haplotype set are variants (i.e. not monomorphic)
		//then there MUST exists at least one other haplotype containing the unobserved alleles,
		//however, this (or these) haplotype(s) are not observed.
		//As such, WLOG we can consider this a case of dimPH=2, which gets reduced to a one dimensional system
		//under reduction, i.e. similar to the dimPH==2 case of above.
				
		double meanXaPHuni = meanXaPH[0];
		double varXaPHuni = gsl_matrix_get(alphaNaMTsigmaGqBmean,0,0);
		double stDevXaPHuni = sqrt(varXaPHuni);	

		logPxAph = -(obsA[0]-meanXaPHuni)*(obsA[0]-meanXaPHuni)/(2*varXaPHuni) - log(stDevXaPHuni*sqrt(2*M_PI)); //log L

	} else { //error
	
		cout << "ERROR in computeLNeutralSelG: dimPH is <1\n";
	}

	//Check of likelihoods are still negative infinity, which should not happen
	if(logPxAph == -numeric_limits<double>::max() && print == false) {
		cout << "ERROR in Misc::computeLNeutralSelGVar: Likelihood for xA is negative infinity. LogPxAph: " << logPxAph <<"\n";

		cout << "Recomputing mean and variance whilst printing out stats:\n";
		computeLNeutralSelGVar(obsA,Nt,ap,qBmean, qBvar, hapFitG, contribsPHset, true);

		cout << "-------------------------------------------------------\n";
	}


	//Clean up before exit
	gsl_vector_free(sigmaGqBmeanGSL);
	gsl_vector_free(TsigmaGqBmeanGSL);

	gsl_matrix_free(T);
	gsl_matrix_free(alphaNaMTsigmaGqBmean);
	gsl_matrix_free(betaMqBmean);
	gsl_matrix_free(gammaqBvar);
	gsl_matrix_free(DsigmaG_qBmean);
	gsl_matrix_free(TDsigmaG_qBmean);
	gsl_matrix_free(three);
	gsl_matrix_free(four);

	if(print == true) {
		cout << "Returning likelihood " << logPxAph << "\n"; 
	}

	return logPxAph;
}





//Takes an already subsetted vector x !
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




//Takes an already subsetted vector x !
void constructMatrixMoriginal(vector<double> x, gsl_matrix *M) {

    //Create subset of vec
    int dim = (int) x.size();
   
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


//Compute the fitness of the full haplotypes fullHaps given the selection model selModel with coefficients selCoefsNew
//SHOULD THIS TAKE INTO ACCOUNT EPISTASIS? Is this method still in use?
vector<double> computeHapFit(vector<Sequence> &fullHaps, Sequence &selModel, vector<double> &selCoefsNew) {
    
    vector<double> hapFit(fullHaps.size(),0);
    for(unsigned int l=0;l<fullHaps.size();l++) { //Loop over haplotypes
        for(int m=0;m<fullHaps[l].getLength();m++) { //For each locus in haplotype
            //Compare nucleotide to selected nucleotide
            if(selModel.getBase(m) != '-') { //Locus under selection
                if(fullHaps[l].getBase(m)==selModel.getBase(m)) { //Haplotype under selection
                    hapFit[l]=hapFit[l]+selCoefsNew[m];
                }
            }
        }
    }
    
    return hapFit;
}

vector<double> computeSimulationHapFitT(vector<Sequence> &fullHaps, SimParamPH *spPH) {

	vector<double> hapFit(fullHaps.size(), 0);
	
	Sequence selSeq = *(spPH->getSel());
	vector<double> selMagVec = *(spPH->getSelMagVec());	

	//Selection coefficents
    for(unsigned int l=0;l<fullHaps.size();l++) { //Loop over haplotypes
        for(int m=0;m<fullHaps[l].getLength();m++) { //For each locus in haplotype
            //Compare nucleotide to selected nucleotide
            if(selSeq.getBase(m) != '-') { //Locus under selection
                if(fullHaps[l].getBase(m)==selSeq.getBase(m)) { //Haplotype under selection
                    hapFit[l]=hapFit[l]+selMagVec[m];
                }
            }
        }
    }

	vector<vector<int> > epiPosVec = *(spPH->getEpiPosVec());
	vector<double> epiMagVec = *(spPH->getEpiMagVec());
	

	//Epistasis
	for(unsigned int i=0; i<fullHaps.size(); i++) { //Loop over haplotypes
	
		cout << "EpiPosVec.size() " << epiPosVec.size() << "\n";
		for(unsigned int j=0; j<epiPosVec.size(); j++) { //Loop over epistatic effecs

			vector<int> currentPos = epiPosVec[j];
			bool epistasisPresent = true; //Assume tru, then check
			for(unsigned int k=0; k<currentPos.size(); k++) { //Loop over positions in epi model
		
				if(fullHaps[i].getBase(currentPos[k]) != selSeq.getBase(currentPos[k])) {

					epistasisPresent = false;
					break;
				}
			}

			if(epistasisPresent == true) {
				hapFit[i] += epiMagVec[j];
			}
			cout << "HapFit[i] = " << hapFit[i] << "\n";
		}
	}
    
    return hapFit;
}

vector<double> computeSimulationHapFitG(vector<Sequence> &fullHaps, SimParamPH *spPH) {

	vector<double> hapFitG(fullHaps.size(), 0);
	
	Sequence selSeqG = *(spPH->getSelG());
	vector<double> selMagVecG = *(spPH->getSelMagVecG());	

	//Selection coefficents
    for(unsigned int l=0;l<fullHaps.size();l++) { //Loop over haplotypes
        for(int m=0;m<fullHaps[l].getLength();m++) { //For each locus in haplotype
            //Compare nucleotide to selected nucleotide
            if(selSeqG.getBase(m) != '-') { //Locus under selection
                if(fullHaps[l].getBase(m)==selSeqG.getBase(m)) { //Haplotype under selection
                    hapFitG[l]=hapFitG[l]+selMagVecG[m];
                }
            }
        }
    }

	vector<vector<int> > epiPosVecG = *(spPH->getEpiPosVecG());
	vector<double> epiMagVecG = *(spPH->getEpiMagVecG());
	

	//Epistasis
	for(unsigned int i=0; i<fullHaps.size(); i++) { //Loop over haplotypes
	
		cout << "EpiPosVecG.size() " << epiPosVecG.size() << "\n";
		for(unsigned int j=0; j<epiPosVecG.size(); j++) { //Loop over epistatic effecs

			vector<int> currentPos = epiPosVecG[j];
			bool epistasisPresent = true; //Assume tru, then check
			for(unsigned int k=0; k<currentPos.size(); k++) { //Loop over positions in epi model
		
				if(fullHaps[i].getBase(currentPos[k]) != selSeqG.getBase(currentPos[k])) {

					epistasisPresent = false;
					break;
				}
			}

			if(epistasisPresent == true) {
				hapFitG[i] += epiMagVecG[j];
			}
			cout << "HapFitG[i] = " << hapFitG[i] << "\n";
		}
	}
    
    return hapFitG;
}

vector<DiploidSequence> createOneLocusCollapsed(vector<vector<DiploidSequence> > &oneLocusModels) {

	vector<DiploidSequence> oneLocusCollapsed;
    for(unsigned int i=0; i<oneLocusModels.size(); i++) { //Loop over selection models
        
        bool locusFound = false;
     
        for(unsigned int j=0; j<oneLocusModels[i].size();j++) { //Loop over genes within sel model
            if(i==0) { //First time around, add an empty sequence (i.e. one for each gene)
                oneLocusCollapsed.push_back(DiploidSequence(oneLocusModels[i][j].getMajor().getLength()));
            }
            
            for(int k=0; k<oneLocusModels[i][j].getMajor().getLength(); k++) { //Loop over length of sequence motif
             
                if(oneLocusModels[i][j].getMajor().getBase(k) != '-') { //If locus is not empty, add it to collapsed
                    
                    if(oneLocusCollapsed[j].getMajor().getBase(k) == '-') { //If collapsed is empty, add it
                        oneLocusCollapsed[j].setMajor(k,oneLocusModels[i][j].getMajor().getBase(k));
                        oneLocusCollapsed[j].setMinor(k,oneLocusModels[i][j].getMinor().getBase(k));
                    } else { //Check that no selection models have different nucleotides at the same positions (this should never happen)
                        
                        if(oneLocusModels[i][j].getMajor().getBase(k) != oneLocusCollapsed[j].getMajor().getBase(k)) {
                            
                            cout << "ERROR in createOneLocusCollapsed! Alleles of one locus models don't match.\n";
                        }
                    }
                    
                    locusFound = true;
                    break; //There should only be one locus per model in oneLocusModels
                }
            }
        }
    }

	return oneLocusCollapsed;   

}

vector<Model> generateModelsFromPreviousBest(vector<DiploidSequence> &collapsedFullHaps, Model &bestModelBefore) {

	vector<Model> result;

	//Loop over each gene and try adding a further selection coefficient or epistasis coefficient to the best model currently (i.e. the one from before)
	for(unsigned int i=0; i<collapsedFullHaps.size(); i++) {

		if(collapsedFullHaps[i].getLength() > 0) { //i.e. there are data for this gene - if not, no point in trying to find new combinations from this gene

			//First, try adding a new selection coefficient
			//Procedure: Find all locus without selection, and add one in turn
			Model::modelSingleGene MSG = bestModelBefore.getModelSingleGene(i);
			DiploidSequence model = MSG.getModel(); //This is the sites currently under selection

			for(int j=0; j<model.getLength(); j++) { //Loop over loci in this gene

				if(model.getMajor().getBase(j) == '-') { //If no selection currently present, add this loci to a model

					DiploidSequence newModel = model;
					newModel.setMajor(j,collapsedFullHaps[i].getMajor().getBase(j));
					newModel.setMinor(j,collapsedFullHaps[i].getMinor().getBase(j));
				
					Model::modelSingleGene newMSG = MSG;
					newMSG.setModel(newModel);
					newMSG.setSelectionPresent(true); //Selection is present at this gene (default==false comes from passing a neutral bestModelBefore)

					Model currentM = bestModelBefore;
					currentM.setModelSingleGene(i,newMSG);
					currentM.countNumParamsToBeFitted();

					result.push_back(currentM);

				}
			}

			//Next, try adding a new epistasis coefficient
			//First find the positions under selection
			vector<int> positionsUnderSel;
			for(int j=0; j<model.getLength(); j++) {

				if(model.getMajor().getBase(j) != '-') { //Selection present
					positionsUnderSel.push_back(j);
				}
			}

			//Next, find all j-way epistasis
			for(unsigned int j=2; j<=positionsUnderSel.size(); j++) { //Requires at least two selection parameters to have epistasis

				//Finds all combos of length j of the elements of positionsUnderSel
				//E.g. if positionsUnderSel has length 5 and j=2, then get (0,1), (0,2), (0,3), (0,4), (1,2), (1,3), (1,4), (2,3), (2,4), (3,4)
				//and if j=3 then get (0,1,2), (0,1,3), (0,1,4), (0,2,3), (0,2,4), (0,3,4), (1,2,3), (1,2,4), (1,3,4), (2,3,4)
				vector<vector<int> > combs = getAllCombs((int) positionsUnderSel.size(), j);

				for(unsigned int k=0; k<combs.size(); k++) { //Loop over combinations

		
					vector<Model::epistasis> epiModelCurrent = MSG.getEpiModel();
					bool combinationFound = false;
					for(unsigned int l=0; l<epiModelCurrent.size(); l++) { //Loop over epistasis coefficients

						if(epiModelCurrent[l].getDim() == (int) j) {


							//Convert a combination into a position e.g. a combination (1,3)
							//refers to a combination of the 1st and 3rd index of positionsUnderSel.
							//So if positionUnderSel = (0,2,3,5,8) then combPos = (2,5)
							vector<int> combPos;
							for(unsigned int m=0; m<combs[k].size(); m++) {
								combPos.push_back(positionsUnderSel[combs[k][m]]);
							}
	
							//Check if current combinations is identical with existing epistasis
							vector<int> epiCurrentPos = epiModelCurrent[l].getPositions();
							if(identicalIntVectors(combPos,epiCurrentPos) == true) { //Combination already included

								combinationFound = true; break;
							}
						}
					}

					if(combinationFound == false) { //Add combination

						Model::epistasis epiCurrent;
						vector<int> epiLociCurrent;
						for(unsigned int l=0; l<combs[k].size(); l++) {
							epiLociCurrent.push_back(positionsUnderSel[combs[k][l]]);
						}
						epiCurrent.setPositions(epiLociCurrent);
						epiCurrent.setDim(j);
						epiModelCurrent.push_back(epiCurrent);
				
						Model::modelSingleGene newMSG = MSG;
						newMSG.setEpiModel(epiModelCurrent);

						Model currentM = bestModelBefore;
						currentM.setModelSingleGene(i,newMSG);
						currentM.countNumParamsToBeFitted();

						result.push_back(currentM);
					}		
				}
			}
		} 
	}


	return result;

}

vector<ModelMR> generateModelsFromPreviousBest(vector<DiploidSequence> &collapsedFullHaps, ModelMR &bestModelBefore) {

	vector<ModelMR> result;
	
	//Loop over each gene and try adding a further selection coefficient or epistasis coefficient to the best model currently (i.e. the one from before)
	for(unsigned int i=0; i<collapsedFullHaps.size(); i++) {

		//First, try adding a new selection coefficient
		//Procedure: Find all locus without selection, and add one in turn
		Model::modelSingleGene MSG = bestModelBefore.getModelSingleGene(i);
		DiploidSequence selModel = MSG.getModel(); //This is the sites currently under selection

		for(int j=0; j<selModel.getLength(); j++) { //Loop over loci in this gene

			if(selModel.getMajor().getBase(j) == '-') { //If no selection currently present, add this loci to a model

				DiploidSequence newModel = selModel;
				newModel.setMajor(j,collapsedFullHaps[i].getMajor().getBase(j));
				newModel.setMinor(j,collapsedFullHaps[i].getMinor().getBase(j));
				
				Model::modelSingleGene newMSG = MSG;
				newMSG.setModel(newModel);
				newMSG.setSelectionPresent(true); //Selection is present at this gene (default==false comes from passing a neutral bestModelBefore)

				ModelMR currentM = bestModelBefore;
				currentM.setModelSingleGene(i,newMSG);
				currentM.countNumParamsToBeFitted();

				result.push_back(currentM);

			}
		}

		//Next, try adding a new epistasis coefficient
		//First find the positions under selection
		vector<int> positionsUnderSel;
		for(int j=0; j<selModel.getLength(); j++) {

			if(selModel.getMajor().getBase(j) != '-') { //Selection present
				positionsUnderSel.push_back(j);
			}
		}

		//Next, find all j-way epistasis
		for(unsigned int j=2; j<=positionsUnderSel.size(); j++) { //Requires at least two selection parameters to have epistasis

			//Finds all combos of length j of the elements of positionsUnderSel
			//E.g. if positionsUnderSel has length 5 and j=2, then get (0,1), (0,2), (0,3), (0,4), (1,2), (1,3), (1,4), (2,3), (2,4), (3,4)
			//and if j=3 then get (0,1,2), (0,1,3), (0,1,4), (0,2,3), (0,2,4), (0,3,4), (1,2,3), (1,2,4), (1,3,4), (2,3,4)
			vector<vector<int> > combs = getAllCombs((int) positionsUnderSel.size(), j);

			for(unsigned int k=0; k<combs.size(); k++) { //Loop over combinations

				vector<Model::epistasis> epiModelCurrent = MSG.getEpiModel();
				bool combinationFound = false;
				for(unsigned int l=0; l<epiModelCurrent.size(); l++) { //Loop over epistasis coefficients
	
					if(epiModelCurrent[l].getDim() == (int) j) {

						vector<int> epiCurrentPos = epiModelCurrent[l].getPositions();
						if(identicalIntVectors(combs[k],epiCurrentPos) == true) { //Combination already included

							combinationFound = true; break;
						}
					}
				}

				if(combinationFound == false) { //Add combination

					Model::epistasis epiCurrent;
					vector<int> epiLociCurrent;
					for(unsigned int l=0; l<combs[k].size(); l++) {
						epiLociCurrent.push_back(positionsUnderSel[combs[k][l]]);
					}
					epiCurrent.setPositions(epiLociCurrent);
					epiCurrent.setDim(j);
					epiModelCurrent.push_back(epiCurrent);
				
					Model::modelSingleGene newMSG = MSG;
					newMSG.setEpiModel(epiModelCurrent);

					ModelMR currentM = bestModelBefore;
					currentM.setModelSingleGene(i,newMSG);
					currentM.countNumParamsToBeFitted();

					result.push_back(currentM);
				}		
			}
		}
	}


	return result;

}



int binomialCoeff(int n, int k) {
    
    int result = 1;
    
    for(int i=1; i<=k; i++) {
        
        result*= ((n+1-i)/(double)i);
    }
    
    return result;
}

int logBinCoeffStirling(int n, int k) {
    
    return n*log(n) - k*log(k) -(n-k)*log(n-k);
}


//Finds the next number from start to max where a set alreadyPicked cannot be chosen from
//The number start is the value of the last entry of alreadyPicked +1
//If alreadyPicked[last entry] == max, then -1 is returned
int getNextLeft(vector<int> alreadyPicked, int max) {

	int start = alreadyPicked[alreadyPicked.size()-1] +1; //Current value of last entry +1 is the starting point
	if(start == max+1) { 
		return -1;
	} else {
		for(int i=start; i<=max; i++) {
			bool toBeIncluded = true; 
			for(unsigned int j=0; j<alreadyPicked.size(); j++) {
				if(i==alreadyPicked[j]) {
					toBeIncluded = false;
				}
			}
			if(toBeIncluded == true) {
				return i;
			}
		}
	}


	return -1; //This should never happen
}

//Find all the sequences in set but NOT in subset
vector<Sequence> findRemainingSeqs(vector<Sequence> subset, vector<Sequence> set) {

	vector<Sequence> result;
	for(unsigned int i=0; i<set.size(); i++) {

		bool found = false;
		for(unsigned int j=0; j<subset.size(); j++) {
			if(identicalSeqs(subset[j],set[i]) == true) {
				found = true;
				break;
			}
		}

		if(found == false) {
			result.push_back(set[i]);
		}
	}
	return result;
}


//Takes a vector of full haplotypes and constructs a single diploidSequence where the major and minor alleles correspond
//to the allels with the highest frequencies
//E.g. if the full haplotypes are ACT, AGT, CCA, AGA with freqs 0.10, 0.20, 0.30 0.40
//then the diploidSequence is major=AGA  minor=CCT
DiploidSequence constructDiploidFromFullHap(vector<Sequence> &fhs, vector<double> &freqs) {

	DiploidSequence result = DiploidSequence(fhs[0].getLength()); //Initialise empty diploid sequence
	for(int i=0; i<result.getLength(); i++) { //Loop over bases
	
		char allele1 = fhs[0].getBase(i); //WLOG set the first allele to be the base at position i from the first haplotype
		char allele2 = '-'; //Currently unknown
		double freqAllele1 = 0;
		double freqAllele2 = 0;
		for(unsigned int j=0; j<freqs.size(); j++) { //Loop over the full haplotype frequencies

			if(fhs[j].getBase(i) != allele1) {
				allele2 = fhs[j].getBase(i);
				freqAllele2 += freqs[j];
			} else {
				freqAllele1 += freqs[j];
			}		
		}

		//Check which is major and which is minor
		if(freqAllele1 >= freqAllele2) {
			result.setMajor(i,allele1);
			result.setMinor(i,allele2);
		} else {	
			result.setMajor(i,allele2);
			result.setMinor(i,allele1);
		}
	}

	return result;
}

//Takes a vector of full haplotypes and constructs a single diploidSequence with major = fhs[0]
DiploidSequence constructDiploidFromFullHap(vector<Sequence> &fhs) {

	if(fhs.size() > 0) { //If there are any full haplotypes

		DiploidSequence result = DiploidSequence(fhs[0].getLength()); //Initialise empty diploid sequence
		for(int i=0; i<result.getLength(); i++) { //Loop over bases
	
			char allele1 = fhs[0].getBase(i);
			char allele2 = '-'; //Currently unknown
			for(unsigned int j=0; j<fhs.size(); j++) { //Loop over the full haplotypes

				if(fhs[j].getBase(i) != allele1) {
					allele2 = fhs[j].getBase(i);
					break;
				}		
			}
		
			result.setMajor(i,allele1);
			result.setMinor(i,allele2);
		}
		
		return result;
	
	} else { //This gene is empty, so return empty diploid sequence

		DiploidSequence empty;
		return empty;
	}
}





//Check that none of the entres of the vector v are the same
////Exception: Zero allowed multiple times.
bool noEntryTheSameExceptZero(vector<double> &v) {


	for(unsigned int i=0; i<v.size()-1; i++) {

		double currentVal = v[i];
		if(currentVal != 0) {
			for(unsigned int j=i+1; j<v.size(); j++) {
		
				if(currentVal == v[j]) { return	false; }
			}	
		}
	}

	return true;
}



//Combinatorics code based on code from Rosetta Code http://rosettacode.org/wiki/Combinations#C.2B.2B
////If N=5 and k=2, then get (0,1), (0,2), (0,3), (0,4), (1,2), (1,3), (1,4), (2,3), (2,4), (3,4)
////and if k=3 then get (0,1,2), (0,1,3), (0,1,4), (0,2,3), (0,2,4), (0,3,4), (1,2,3), (1,2,4), (1,3,4), (2,3,4)
vector<vector<int> > getAllCombs(int N, int k) {

	vector<vector<int> > result;

	string bitmask(k, 1); // K leading 1's
	bitmask.resize(N, 0); // N-K trailing 0's
	//Store integers and permute bitmask
	do {
		vector<int> currentComb;
        	for (int i = 0; i < N; ++i) { // [0..N-1] integers
			if (bitmask[i]) {
				currentComb.push_back(i);
			}
		}
		result.push_back(currentComb);
	} while (prev_permutation(bitmask.begin(), bitmask.end()));

	return result;
}



vector<vector<int> > getLinkedPairs(vector<Sequence> &haps) {

	vector<vector<int> > result;

	int numLoci = haps[0].getLength();
	for(int i=0; i<numLoci-1; i++) {
		for(int j=i+1; j<numLoci; j++) {

			bool linked = true;
			char i0 = haps[0].getBase(i);
			char j0 = haps[0].getBase(j);
			
			//Go through all haps and see if linked
			for(unsigned int k=1; k<haps.size(); k++) {

				if(haps[k].getBase(i) == i0 && haps[k].getBase(j) != j0) {
					linked = false;
					break;

				} else if(haps[k].getBase(i) != i0 && haps[k].getBase(j) == j0) {
					linked = false;
					break;
				}
			}


			if(linked == true) { //Position i and j is a linked pair

				vector<int> pair = {i,j};
				result.push_back(pair);
			}
		}
	}


	return result;
}

vector<string> split(const string& s, char delim) {
	stringstream ss(s);
	string item;
	vector<string> tokens;
	while(getline(ss,item,delim)) {
		tokens.push_back(item);
	}
	return tokens;
}


//Related to gaining information on memory usage
int parseLine(char* line){
	// This assumes that a digit will be found and the line ends in " Kb".

	 int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

int getValue(){ //Note: this value is in KB!

FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmSize:", 7) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}


bool stringToBool(string& s) {


	if(s.compare("true") == 0 || s.compare("True") == 0 || s.compare("TRUE") == 0 || s.compare("1") == 0) { return true; }
	else if(s.compare("false") == 0 || s.compare("False") == 0 || s.compare("FALSE") == 0 || s.compare("0") == 0) { return false; }
	else{
		cout << "ERROR! argument: " << s << " in Misc:stringToBool(string& s) is not of the right form. Returning false by default.";
		return false;
	}
}

//Compares two vectors by size and returns true if b>a
bool compareBySize(const vector<int>& a, const vector<int>& b) {

	return (a.size() > b.size());

}


//Gamma and delta functions for compound distributions

//Computes 
//\begin{equation}
//\gamma_n=\frac{(-1)^n \lambda ^{-\frac{n^2}{2}-\frac{n}{2}} \left(N^T\right)^{-n-1} (\lambda 
//   N^T;\lambda )_n \left(\lambda  \left(N^T\right)^2 \left(\sum\limits_{j=0}^{n-1}
//      \frac{(\lambda  N^T-1) (-1)^{1-j} \lambda^{\frac{j^2}{2}+\frac{j}{2}-1} \left(N^T\right)^{j-1}}{(\lambda  N^T;\lambda)_{j+1}}\right)+\lambda  N^T-1\right)}{\lambda  N^T-1}
//      \end{equation}
//
// where lambda = growth factor
double gamma_n(int n, double lambda, int Nt) {

	double gamma = pow(-1,n)*pow(lambda,-(n*n)/2.0-n/2.0)*pow(Nt,-n-1)*qPochhammer(lambda*Nt,lambda,n)/(lambda*Nt-1); //Consider if division by integers is a problem for large lambda,Nt

	double sum = 0;
	for(int j=0; j<=n-1; j++) {

		sum += (lambda*Nt-1)*pow(-1,1-j)*pow(lambda,j*j/2.0+j/2.0-1)*pow(Nt,j-1)/qPochhammer(lambda*Nt,lambda,j+1);
	}

	gamma *= (lambda*Nt*Nt*sum + lambda*Nt -1);


	return gamma;
}

//Computes 
//\begin{equation}
//\delta_n=(-1)^n \lambda ^{-\frac{n^2}{2}-\frac{3 n}{2}} \left(N^T\right)^{-n-1} \left(N^T \lambda^n-1\right) (\lambda  N^T;\lambda )_n
//\end{equation}
//
// where lambda = growth factor
double delta_n(int n, double lambda, int Nt) {

	double delta = pow(-1,n)*pow(lambda,-(n*n)/2.0-n/2.0)*pow(Nt,-n-1)*(Nt-1)*qPochhammer(lambda*Nt,lambda,n);

	return delta;
}

//q-Pochhammer symbol: (a,q)_n = prod_0^{n-1} (1-a*q^k), with (a,q)_0 = 1
double qPochhammer(double a, double q, int n) {

	if(n == 0) {

		return 1;

	} else if(n > 0) {

		double prod = 1;
		for(int k=0; k<=n-1; k++) {

			prod *= (1-a*pow(q,k));
		}
		return prod;
	} else {

		cout << "Error in qPochhammer function: n is < 0: " << n << "\n";
		exit(1);
	}

}


















