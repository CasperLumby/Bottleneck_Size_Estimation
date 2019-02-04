#include "misc.h"
#include <math.h>
#include <iostream>


#ifdef _OPENMP
        #include <omp.h>
#else
        #define omp_get_thread_num() 0
#endif

using namespace std;


int main (int argc, char* argv[]) {

	//Test 1
	vector<double> v1 = { 0.223435 };
	int dim1 = (int) v1.size();
	gsl_matrix* m1a = gsl_matrix_alloc(dim1,dim1); 
	constructMatrixMoriginal(v1,m1a);
	gsl_matrix* m1b = gsl_matrix_alloc(dim1,dim1); 
	constructMatrixM(v1,m1b);
	cout << "Given vector "; printDoubleVector(v1);
	cout << "the M matrix is (original method):\n"; printMatrixMathematica(m1a);
	cout << "the M matrix is (hard-coded):\n"; printMatrixMathematica(m1b);
	cout << "\n";
	

	//Testing timings
	//Get start time
	double timeStart1A = omp_get_wtime();

	for(int i=0; i<1000000; i++) { //1million rounds

		vector<double> vi = { 0.223435 };
		int dimi = (int) vi.size();
		gsl_matrix* mia = gsl_matrix_alloc(dimi,dimi); 
		constructMatrixMoriginal(vi,mia);
	}
	
	//Print elapsed time to file 
	double timeFinish1A = omp_get_wtime();
	cout << "Elapsed time 1A: " << timeFinish1A - timeStart1A << "\n";

	//Get start time
	double timeStart1B = omp_get_wtime();

	for(int i=0; i<1000000; i++) { //1million rounds

		vector<double> vi = { 0.223435 };
		int dimi = (int) vi.size();
		gsl_matrix* mia = gsl_matrix_alloc(dimi,dimi); 
		constructMatrixM(vi,mia);
	}
	
	//Print elapsed time to file 
	double timeFinish1B = omp_get_wtime();
	cout << "Elapsed time 1B: " << timeFinish1B - timeStart1B << "\n\n";



	//Test 2
	vector<double> v2 = { 0.99999999918647342056488014350179582834243774414062, 0.000001123 };
	int dim2 = (int) v2.size();
	gsl_matrix* m2a = gsl_matrix_alloc(dim2,dim2); 
	constructMatrixMoriginal(v2,m2a);
	gsl_matrix* m2b = gsl_matrix_alloc(dim2,dim2); 
	constructMatrixM(v2,m2b);
	cout << "Given vector "; printDoubleVector(v2);
	cout << "the M matrix is (original method):\n"; printMatrixMathematica(m2a);
	cout << "the M matrix is (hard-coded):\n"; printMatrixMathematica(m2b);
	cout << "\n";
	

	//Testing timings
	//Get start time
	double timeStart2A = omp_get_wtime();

	for(int i=0; i<1000000; i++) { //1million rounds

		vector<double> vi = { 0.99999999918647342056488014350179582834243774414062, 0.000001123 };
		int dimi = (int) vi.size();
		gsl_matrix* mia = gsl_matrix_alloc(dimi,dimi); 
		constructMatrixMoriginal(vi,mia);
	}
	
	//Print elapsed time to file 
	double timeFinish2A = omp_get_wtime();
	cout << "Elapsed time 2A: " << timeFinish2A - timeStart2A << "\n";

	//Get start time
	double timeStart2B = omp_get_wtime();

	for(int i=0; i<1000000; i++) { //1million rounds

		vector<double> vi = { 0.99999999918647342056488014350179582834243774414062, 0.000001123 };
		int dimi = (int) vi.size();
		gsl_matrix* mia = gsl_matrix_alloc(dimi,dimi); 
		constructMatrixM(vi,mia);
	}
	
	//Print elapsed time to file 
	double timeFinish2B = omp_get_wtime();
	cout << "Elapsed time 2B: " << timeFinish2B - timeStart2B << "\n\n";





	return 0;
}

