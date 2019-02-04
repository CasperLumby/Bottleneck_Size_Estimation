#include "misc.h"
#include <math.h>
#include <iostream>
//#include "../cubature-master/cubature.h"
#include "cubature.h"

#if defined(PCUBATURE)
#  define cubature pcubature
#else
#  define cubature hcubature
#endif

#ifdef _OPENMP
        #include <omp.h>
#else
        #define omp_get_thread_num() 0
#endif

using namespace std;

struct parameters {
		vector<double> mu;
		gsl_matrix* Sigma;
};
int count = 0;

//Integration function
int f(unsigned dim, const double *x, void *fdata, unsigned int fdim, double *fval) {

	count++; //Increase number of evaluations
	parameters* params = (parameters*) (fdata); //Cast fval to a parameters object			
	vector<double> xvals(dim);
	for(unsigned int i=0; i<dim; i++) { xvals[i] = x[i]; }

	*fval = exp(logMultivariateNormalPDF(xvals, params->mu , params->Sigma));
	
	return  0;

}

int main (int argc, char* argv[]) {

	//Get start time
	double timeStart = omp_get_wtime();


	cout << "Program started.\n";

	vector<double> mu = {10,20,30,30,20,10};
	gsl_matrix* Sigma = gsl_matrix_alloc(mu.size(),mu.size());
	gsl_matrix_set_zero(Sigma);
	for(unsigned int i=0; i<mu.size(); i++) {
		gsl_matrix_set(Sigma,i,i,300);
	}

	
	parameters params;
	params.mu = mu;
	params.Sigma = Sigma;
	void *paramPtr = &params;
	unsigned int dim = mu.size();
	unsigned int fdim = 1;

	//Integration ranges
	double* xmin = (double *) malloc(dim * sizeof(double));
     	double* xmax = (double *) malloc(dim * sizeof(double));
     	for (int i = 0; i < ((int) dim); ++i) {
	  xmin[i] = 0;
	  xmax[i] = 50;
     	}

	//Other definitions
	double tol =1e-3;
	size_t maxEval = 0;
	double *val, *err;
  val = (double *) malloc(sizeof(double) * fdim);
     err = (double *) malloc(sizeof(double) * fdim);
     

	printf("%u-dim integral, tolerance = %g\n", dim, tol);
	cubature(fdim, f, paramPtr, dim, xmin, xmax, maxEval, tol, tol, ERROR_INDIVIDUAL, val, err);
//	hcubature(fdim, f, paramPtr, dim, xmin, xmax, maxEval, 0, 0, ERROR_INDIVIDUAL, val, err);
   
	for (int i = 0; i < ((int) fdim); ++i) {
		printf("integral = %0.11g and error: %0.11g \n", val[i], err[i]);
		printf("log(integral): %0.11g\n", log(val[i]));
//		printf("integrand %d: integral = %0.11g, est err = %g, true err = %g\n", 
//		which_integrand[i], val[i], err[i], 
//		fabs(val[i] - exact_integral(which_integrand[i], dim, xmax)));
	}
	printf("#evals = %d\n", count);


	//Print elapsed time to file 
	double timeFinish = omp_get_wtime();
	cout << "Elapsed time: " << timeFinish - timeStart << "\n";

	return 0;
}

