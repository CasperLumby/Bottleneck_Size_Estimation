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
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
using namespace std;

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

int main(int argc, char* argv[])
{
	string qStar = argv[1];
	freopen("maximum_likelihood.txt","w",stdout); //print out the progress of the optimisation process in a file called update.txt

	if (FILE *file1 = fopen(qStar.c_str(), "r")) 
	{
        fclose(file1);
		vector<string> table;
		ifstream test;
		test.open(argv[1]); //get the reconstructed haplotypes in qStar
		while (!test.eof())
		{
			string line;
			getline(test, line);
			table.push_back(line);
		}
		if (table.empty())
		{
			return false;
		}
		test.close();
		vector < vector<string> > rows;
		for (unsigned int i = 0; i < table.size(); i++)
		{
			vector<string> row;
			istringstream iss(table[i]);
			string term;
			while (iss >> term)
			{
				row.push_back(term);
			}
			if (!row.empty())
			{
				rows.push_back(row);
			}
		}
		
		vector<double> q_star_b;
		vector<double> q_star_a;
		vector<string> temp;
		vector<string> matching_seuquence; // we record the haplotype number in outcome_1.txt to then match them with their corresponding qStarStar frequency 
		for (unsigned int i = 0; i < rows.size(); i++)
		{
			temp = (rows[i]);
			unsigned int sz = temp.size();
			string freq_b = rows[i][sz-2];
			string freq_a = rows[i][sz-1];
			double frequency_a = atof(freq_a.c_str());
			double frequency_b = atof(freq_b.c_str());
			q_star_b.push_back(frequency_b);
			q_star_a.push_back(frequency_a);
			matching_seuquence.push_back(rows[i][sz-4]);
		}
		
		vector<string> table_2;
		ifstream test_2;
		test_2.open(argv[2]); //get the full list of catted test_i.txt files that are currently stored in Transmission1_bottleneck/segment_i/
		while (!test_2.eof())
		{
			string line_2;
			getline(test_2, line_2);
			table_2.push_back(line_2);
		}
		test_2.close();
		vector < vector<string> > rows_2;
		for (unsigned int i = 0; i < table_2.size(); i++)
		{
			vector<string> row_2;
			istringstream iss(table_2[i]);
			string term_2;
			while (iss >> term_2)
			{
				row_2.push_back(term_2);
			}
			if (!row_2.empty())
			{
				rows_2.push_back(row_2);
			}
		}
		
		
		double qStarStar_size = (rows_2.size())/(rows.size()); //this gives the total number of q**s generated for each haplotype in outcome_1.txt
		//cout << "Size of qStarStar = " << qStarStar_size << endl;
		
		vector<double> q_star_star_b(matching_seuquence.size(), 0); //this file calculates the variance for each inferred haplotype frequency q* before transmission (i.e. first timepoint)
		vector<double> q_star_star_a(matching_seuquence.size(), 0); //this file calculates the variance for each inferred haplotype frequency q* after transmission (i.e. second timepoint)
		int count = 0;
		for (unsigned int i = 0; i < rows_2.size(); i++)
		{
			count ++;
			temp = (rows_2[i]);
			unsigned int sz = temp.size();
			for (unsigned int j = 0; j < matching_seuquence.size(); j++) //for every element of qStarStar, we go over all the qStar haplotypes and match the corresponding haplotypes to each other.
			{
				int number = stoi(matching_seuquence[j]);
				number -= 1; //outcome_1.txt starts numbering haplotypes from 1 whereas test_i.txt starts from 0 --> the indexing for outcome_1.txt should be shifted down by 1 unit.
				string hap_num = to_string(number);
				if (temp[0] == hap_num) // if the two haplotype numbers are the same
				{
					string temp_1 = temp[sz-2];
					double num_1 = atof(temp_1.c_str());
					string temp_2 = temp[sz-1];
					double num_2 = atof(temp_2.c_str());
					q_star_star_b[j] += (pow((num_1 - q_star_b[j]), 2))*1/qStarStar_size;
					q_star_star_a[j] += (pow((num_2 - q_star_a[j]), 2))*1/qStarStar_size;
				}
			}
		}
		
		vector<double> likelihood(1000, 1); // calculate the log-likelihood given the bottleneck size could vary anywhere between 1 and 1000
		for (unsigned int  n = 0; n < 1000; n++)
		{
			for (unsigned int i = 0; i < q_star_a.size(); i++)
			{
				double nn = n + 1;
				double x = q_star_a[i];
				double m = q_star_b[i]; //mean of the Gaussian
				double s = q_star_star_a[i] + (1-1/(22*nn))*((1/nn)*(m*(1 - m)) + (1 - 1/nn)*q_star_star_b[i]) + (1/(22*nn))*(m*(1 - m)); //variance of the Gaussian
				s = sqrt(s);
				likelihood[n] += log((1/(sqrt(2*3.141592653*s*s)))*exp(-0.5*((x - m)/s)*(x - m)/s)); // the likelihood is approximated as a Gaussian and the log-likelihood of each haplotype is added together
			}
		}
		
		double Max_finder = -1e20;
		double Max_finder_count = 0; //shows at what bottleneck size the likelihood function is maximised
		for (unsigned int i = 0; i < likelihood.size(); i++)
		{
			if (likelihood[i] > Max_finder)
			{
				Max_finder_count = i+1;
				Max_finder = likelihood[i];
			}
		}
		
		if (Max_finder_count != 0)
		{	
			cout << "Bottleneck size with maximum likelihood is = " << Max_finder_count << endl;
		}
		
		
		if (Max_finder != 1 && Max_finder != 0)
		{
			ofstream myfile;
			myfile.open(argv[3]); // save the likelihood function
			for (unsigned int i = 0; i < likelihood.size(); i++)
			{
				if (likelihood[i] > -827) //this is the lowest the log-likelihood can get. If it goes any lower, cpp shows -inf.
				{
					myfile << likelihood[i] << endl;
				}
				else
				{
					myfile << -827 << endl;
				}
			}
			myfile.close();
		}
		
		return 0;
	}
	
	else
	{
		cout << "This segment contains no data. Either outcome_1.txt and/or test_*.txt are empty" << endl;
		return false;
	}

} 