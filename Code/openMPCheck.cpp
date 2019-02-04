#include <iostream>
#include <vector>
#include <omp.h>
#include <unistd.h>

using namespace std;

int main () {

	
	vector<int> a(20,-1);


#pragma omp parallel default(none) shared(a, cout)
{

	#pragma omp for	schedule(runtime)
	for(unsigned int i=0; i<a.size(); i++) {
		int TID = omp_get_thread_num();
		a[i] = i*i;

		
		#pragma omp critical
		{
		cout << "Value of i: " << i << "\tCurrent TID: " << TID << "\n";
		}

		for(int j=0; j<((int) i); j++) {
			sleep(1);
		}

	}

}
	cout << "Final vector is: \n";
	for(unsigned int i=0; i<a.size(); i++) {
		cout << a[i] << " ";
	}
	cout << "\n";

	return 0;
}
