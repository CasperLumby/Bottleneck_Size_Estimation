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
	//first check if the input format is right
	if(argc != 4) { 
		cout << "Error in the input format: " << endl << "Please provide the names for your concatenated likelihood file, full likelihood distribution, and the bottleneck size with maximum likelihood value, respectively." << endl;
		return false;
	}
	
	if(argc == 4) { 
		cout << "avg_bottleneck is now initiated..." << endl;
	}
	
	freopen(argv[3],"w",stdout); //This writes the stdout output to Transmission1_maximum_likelihood.txt
	vector<string> likelihoodsString; //Likelihoods in string form
	vector<double> likelihoods; //Likelihoods as doubles
	ifstream bottleneck;
	bottleneck.open(argv[1]); // a concatenated list of all likelihoods for one Transmission event with seven segments
	while (!bottleneck.eof())
	{
		string line;
		getline(bottleneck, line);
		likelihoodsString.push_back(line);
		likelihoods.push_back(atof(line.c_str()));
	}
	bottleneck.close();
	
	vector<double> likelihood(1000, 1);
	double index = 0;	
	for (unsigned int i = 0; i < likelihoods.size(); i++)
	{
		if((i+1)%1000 != 0) //every 1000 row belong to the likelihood of one segment.
		{
			likelihood[index] += likelihoods[i];
			index++;
		}
		
		else if((i+1)%1000 == 0)
		{
			likelihood[index] += likelihoods[i];
			index = 0;
		}
	}
	
	double Max_finder = std::numeric_limits<double>::lowest();
	double Max_finder_index = 0;
	for (unsigned int i = 0; i < likelihood.size(); i++)
	{
		if (likelihood[i] > Max_finder)
		{
			Max_finder_index = i+1;
			Max_finder = likelihood[i];
		}
	}
	
	if (Max_finder_index!=0)
	{	
		cout << Max_finder_index << endl; //This is the mean bottleneck size of the transmission pair and the information is saved in argv[3].txt
	}
	
	ofstream myfile;
	myfile.open(argv[2]); //This is the full likelihood distribution of the bottleneck estimates (from 1 to 1000) for this transmission pair
	for (unsigned int i = 0; i < likelihood.size(); i++)
	{
		myfile << likelihood[i] << endl;
	}
	myfile.close();


	return 0;
} 