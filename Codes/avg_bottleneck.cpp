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
#include <time.h> 
#include <limits>

using namespace std;


int main(int argc, char* argv[])
{
	freopen(argv[3],"w",stdout); //print out the progress of the optimisation process in a file called update.txt
	vector<string> table;
	ifstream bottleneck;
	bottleneck.open(argv[1]); // a catted list of all likelihoods for one Transmission event with seven segments
	while (!bottleneck.eof())
	{
		string line;
		getline(bottleneck, line);
		table.push_back(line);
	}
	bottleneck.close();
	
	vector <double> row;
	for (unsigned int i = 0; i < table.size(); i++)
	{
		istringstream iss(table[i]);
		string term;
		while (iss >> term)
		{
			double num_1 = atof(term.c_str());
			row.push_back(num_1);
		}
	}
	
	vector<double> likelihood(1000, 1);
	double counter = 0;	
	for (unsigned int i = 0; i < row.size(); i++)
	{
		if((i+1)%1000 != 0) //every 1000 row belong to the likelihood of one segment.
		{
			likelihood[counter] += row[i];
			counter++;
		}
		
		else if((i+1)%1000 == 0)
		{
			likelihood[counter] += row[i];
			counter = 0;
		}
	}
	
	double Max_finder = -1e9;
	double Max_finder_count = 0;
	for (unsigned int i = 0; i < likelihood.size(); i++)
	{
		if (likelihood[i] > Max_finder)
		{
			Max_finder_count = i+1;
			Max_finder = likelihood[i];
		}
	}
	
	if (Max_finder_count!=0)
	{	
		cout << Max_finder_count << endl;
	}
	
	ofstream myfile;
	myfile.open(argv[2]); //
	for (unsigned int i = 0; i < likelihood.size(); i++)
	{
		myfile << likelihood[i] << endl;
	}
	myfile.close();


	return 0;
} 