#include <algorithm>
#include <iostream>
#include <list>
#include <string>
#include <array>
#include <fstream>
#include <math.h> 
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

struct container_is_empty_t {
	template<class C>
	bool operator()(C const& c)const {
		return c.empty();
	}
	template<class T, size_t N>
	bool operator()(T(&)[N])const {
		return N > 0;
	}
};
static container_is_empty_t const container_is_empty;

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

//Calculate the log Dirichlet Multinomial given c=noise parameter, n=reads, x=frequencies of haplotypes
double Likelihood(int c, vector<double> n, vector<double> x)
{
	double sum_1 = 0;
	for (unsigned int i = 0; i < n.size(); i++)
	{
		sum_1 += n[i];
	}
	double sum_2 = 0;
	for (unsigned int i = 0; i < n.size(); i++)
	{
		sum_2 += gsl_sf_lngamma(n[i] + 1);
	}
	double sum_3 = 0;
	for (unsigned int i = 0; i < x.size(); i++)
	{
		sum_3 += (c*x[i]);
	}
	double sum_4 = 0;
	for (unsigned int i = 0; i < x.size(); i++)
	{
		sum_4 += (n[i] + c * x[i]);
	}
	double sum_5 = 0;
	for (unsigned int i = 0; i < n.size(); i++)
	{
		sum_5 += gsl_sf_lngamma(n[i] + c * x[i]);
	}
	double sum_6 = 0;
	for (unsigned int i = 0; i < n.size(); i++)
	{
		sum_6 += gsl_sf_lngamma(c*x[i]);
	}
	return gsl_sf_lngamma(sum_1 + 1) - sum_2 + gsl_sf_lngamma(sum_3) - gsl_sf_lngamma(sum_4) + sum_5 - sum_6;
}

/*Given a timepoint (0 for the donor and 1 for the recipient), 
this function finds the likelihood of a set of haplotypes with 
their frequencies (freq) and their contribution to the reads of 
a given partial haplotype set (CONTRIBS) with full partitioning 
of the partial sets (partitions) and inferred noise parameter c*/
double Dirichlet_m(double timepoint, vector<string> haplotypes, vector<double> freq, vector<vector<vector<int>>> &CONTRIBS, vector<vector<vector<string>>> &partitions, double c)
{
	//In Multi_locus_trajectories.out, go two columns back from the last column on the right --> first timepoint 
	if (timepoint == 0)
	{
		timepoint = 3;
	}
	//In Multi_locus_trajectories.out, go zero columns back from the last column on the right --> last timepoint 
	else if (timepoint == 1)
	{
		timepoint = 1;
	}
	vector<string> candidates = haplotypes;
	vector<double> q = freq;
	double ULTIMA_THULE = 0; //For a candidate haplotype set, this parameter calculates the sum of the likelihoods for each of the i partial haplotype sets
	int candidates_size = candidates.size();
	for (unsigned int i = 0; i<CONTRIBS.size(); i++)
	{
		vector<double> inf; //contains the sum of the frequencies of all contributing candidate haplotypes for each of the j members of the ith partial haplotype set
		vector<double> nn; //contains the corresponding total number of reads for each member of inf vector
		vector<double> temp = q;
		temp.erase(remove(temp.begin(), temp.end() - 1, q[candidates_size]));//erase the unkknown haplotype qx out of the frequency list 
		unsigned int check = 0;
		vector<double> majorCheckvec;
		double NUM1 = 0;
		double sum_nn = 0;
		unsigned int majorCheck1 = 0;
		vector<double> vaccum_holder;
		for (unsigned int j = 0; j<(CONTRIBS[i]).size(); j++)
		{
			double sum = 0;
			int majorCheck = 0;
			for (unsigned int k = 0; k<(CONTRIBS[i][j]).size(); k++)
			{
				//add the frequency of all the haplotypes that share the same partial haplotype of interest
				if (!(CONTRIBS[i][j]).empty())
				{
					sum += q[(CONTRIBS[i][j][k])];//add the frequency of all haplotypes that contribute to the ith partial haplotype set (CONTRIBS[i][j][k] = ith partial haplotype set, jth element, kth contributing haplotype)
					majorCheck++;
					majorCheck1++;
					majorCheckvec.push_back(q[(CONTRIBS[i][j][k])]);
				}
			}
			//if none of the candidate haplotypes match with the partial haplotype of interest (i.e. do not 'contribute' to that partial haplotype set), the corresponding reads for that partial haplotype goes to the unknown haplotype qx)
			if ((CONTRIBS[i][j]).empty())
			{
				check++;
				double SIZE1_1 = (partitions[i][j]).size();
				NUM1 = atoi((partitions[i][j]).at(SIZE1_1 - timepoint).c_str());
				sum_nn += NUM1;
				vaccum_holder.push_back(NUM1);
			}
			//if, at least, there were one candidate haplotype that contributed to the ith partial haplotype set, uptade nn and inf
			if (majorCheck != 0)
			{
				inf.push_back(sum);
				double SIZE = (partitions[i][j]).size();
				double NUM = atoi((partitions[i][j]).at(SIZE - timepoint).c_str());
				nn.push_back(NUM);
			}
		}
		//if the candidate haplotypes covered all the partial haplotypes in the ith set, then there are zero reads left to be attributed to qx
		if (check == 0 && majorCheck1 == temp.size())
		{
			inf.push_back(q[candidates_size]);
			nn.push_back(0.0);
		}
		unsigned int check2 = 0;
		
		//if some (but not all) of the candidate haplotypes covered the entire partial haplotype set i, the remaining candidates would correspond to zero reads (so as the qx haplotype)
		if (check == 0 && majorCheck1 != temp.size())
		{
			for (unsigned int m = 0; m < majorCheckvec.size(); m++)
			{
				temp.erase(remove(temp.begin(), temp.end() - 1, majorCheckvec[m]));
			}
			double temp_sum = 0;
			for (unsigned int g = 0; g < temp.size(); g++)
			{
				temp_sum += temp[g];
			}
			inf.push_back(temp_sum);
			nn.push_back(0.0);
			inf.push_back(q[candidates_size]);
			nn.push_back(0.0);
			check2++;
		}
		
		//if some (but not all) of the candidate haplotypes did not contribute to the ith set, then they correspond to zero reads and the remaining unidentified partial haplotypes go to qx
		unsigned int check3 = 0;
		if (check != 0 && majorCheck1 != temp.size() && majorCheck1 != 0 && check2 == 0)
		{
			for (unsigned int m = 0; m < majorCheckvec.size(); m++)
			{
				temp.erase(remove(temp.begin(), temp.end() - 1, majorCheckvec[m]));
			}
			double temp_sum = 0;
			for (unsigned int g = 0; g < temp.size(); g++)
			{
				temp_sum += temp[g];
			}
			double total = accumulate(vaccum_holder.begin(), vaccum_holder.end(), 0);
			inf.push_back(temp_sum);
			nn.push_back(0.0);
			inf.push_back(q[candidates_size]);
			nn.push_back(total);
			check3++;
		}
		//if non of the candidate haplotypes contributed to ith set, then they all correspond to zero reads and everything else would be classified as qx
		if (check != 0 && majorCheck1 != temp.size() && majorCheck1 == 0 && check2 == 0 && check3 == 0)
		{
			double temp_sum = 0;
			for (unsigned int g = 0; g < temp.size(); g++)
			{
				temp_sum += temp[g];
			}
			double total = accumulate(vaccum_holder.begin(), vaccum_holder.end(), 0);
			inf.push_back(q[candidates_size]);
			nn.push_back(total);
			inf.push_back(temp_sum);
			nn.push_back(0.0);
		}

		//if all the candidate haplotypes covered some of the partial haplotypes, then all the remaining reads correspond to qx.
		if (check != 0 && majorCheck1 == temp.size() && check2 == 0 && check3 == 0)
		{
			inf.push_back(q[candidates_size]);
			nn.push_back(sum_nn);
		}
		ULTIMA_THULE += Likelihood(c, nn, inf);
		temp = q;
	}
	return ULTIMA_THULE;
}

//For a candidate haplotype set, this function finds an optimal frequency given the short-read data collected from Multi_locus_trajectories.out file and returns the log-likelihood value for that set
double OptimumFreqs(double timepoint, vector<string> &candidates, vector<double> &init_freqs, vector<vector<vector<int>>> &CONTRIBS, vector<vector<vector<string>>> &partitions, double c)
{
	gsl_rng_env_setup();
	gsl_rng *rgen = gsl_rng_alloc(gsl_rng_taus);
	int seed = (int)time(NULL);
	gsl_rng_set(rgen, seed);
	double init_size = init_freqs.size();
	vector<double> init_freqs_store = init_freqs;
	double L = -1e7; //initial (extremely low) value for the likelihood
	unsigned int check = 1;
	unsigned int max_it = 8e3; //maximum number of attempts for optimising the frequencies
	double changex = 1e-2; //magnitude of the incremental (random) change in frequency
	double L_store1 = -1e7; //improved likelihood after frequency optimisation

	for (unsigned int it = 0; it < max_it; it++)
	{
		double r = 2 * (gsl_rng_uniform(rgen) - 0.5);//random direction between -1 to +1
		int j = floor(init_freqs.size()*gsl_rng_uniform(rgen));
		init_freqs[j] = init_freqs[j] + (r*changex);

		for (unsigned int i = 0; i < init_size; i++)
		{
			if (init_freqs[i] < 1e-11)
			{
				//Warning: if you change this threshold below from 1e-11 to 0, there will be a problem with calculating the gamma function
				//we assume that the lowest possible frequency is 10^-11
				init_freqs[i] = 1e-11;
			}
		}

		NormaliseFreqs(init_freqs); //frequencies should add up to one after being randomly changed
		if (init_freqs[init_size - 1] > 1e-2)
		{
			init_freqs[init_size - 1] = 1e-2; //the frequency of the 'X' haplotype, qx, could only go up to 1%
			double TOTAL = 0;
			for (unsigned int i = 0; i < init_size-1; i++)
			{
				TOTAL += init_freqs[i];
			}
			for (unsigned int i = 0; i < init_size-1; i++)
			{
				init_freqs[i] = init_freqs[i]*(0.99/TOTAL);//in this case, the frequency of the rest of the haplotypes should add up to 99%
			}
		}
		L = Dirichlet_m(timepoint, candidates, init_freqs, CONTRIBS, partitions, c);//calculate the Dirichlet likelihood after frequency optimisation step 

		if (L > L_store1)//if the likelihood improved, store it and go from there in the next step
		{
			L_store1 = L;
			init_freqs_store = init_freqs;
			check = 1;
		}
		if (L <= L_store1)//if there was no improvement for over 50 consecutive attempts, then we have reached a maximum likelihood peak 
		{
			init_freqs = init_freqs_store;
			check++;
		}
		if (check >= 50)//50 is found heuristically from testing multiple simulations
		{
			break;
		}
	}
	return L_store1;
}

int main(int argc, char* argv[])
{
	string xStar = argv[1];
	string reconstructed_haplotypes = argv[2];
	double c = atoi(argv[5]); //fifth input is the inferred noise parameter C (same as MLHapRec)

	if (FILE *file1 = fopen(xStar.c_str(), "r")) 
	{
		if (FILE *file2 = fopen(reconstructed_haplotypes.c_str(), "r"))
			{
			fclose(file1);
			vector<string> table_2;
			ifstream test_2;
			//Change the directory below
			test_2.open(argv[2]);
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
			
			vector<double> f_before;
			vector<double> f_after;
			vector<string> temp;
			vector<string> candidates;
			for (unsigned int i = 0; i < rows_2.size(); i++)
			{
				temp = (rows_2[i]);
				unsigned int sz = temp.size();
				candidates.push_back(rows_2[i][sz-3]);
				string temp_1 = rows_2[i][sz-2];
				double num_1 = atof(temp_1.c_str());				
				f_before.push_back(num_1);
				string temp_2 = rows_2[i][sz-1];
				double num_2 = atof(temp_2.c_str());
				f_after.push_back(num_2);
			}
			
			vector<string> table;
			ifstream reads;
			reads.open(argv[1]);
			while (!reads.eof())//read the content of Multi_locus_trajectories.out and store it in table
			{
				string line;
				getline(reads, line);
				table.push_back(line);
			}
			reads.close();
			
			vector < vector<string> > rows; //save each row of Multi_locus_trajectories.out
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
			
			vector<string> loci; //collect the loci numbers from Multi_locus_trajectories.out
			vector<string> temp2;
			for (unsigned int i = 0; i < rows.size(); i++)
			{
				temp2 = (rows[i]);
				unsigned int sz = temp2.size();
				for (unsigned int j = 1; j < sz - 6; j++) //according to the current format of SAMFIRE, from the second column all the way to the seventh last column (from the right) contains loci numbers
				{
					loci.push_back(rows[i][j]);
				}
			}
			
			int TOTAL_N = 0; //calculating the total number of reads for all the partial haplotypes -- this quantity is required for BIC calculations 
			for (unsigned int i=0; i<rows.size(); i++)
			{
				unsigned int sz = rows[i].size();
				string temp_1 = rows[i][sz-1];
				string temp_2 = rows[i][sz-3];
				int num_1 = atoi(temp_1.c_str());
				int num_2 = atoi(temp_2.c_str());
				TOTAL_N += num_1 + num_2;
			}
			
			vector<int> positions;//The following ~30 lines of code sorts the order of loci from 0 to the last variant locus and store them in Fpositions 
			for (unsigned int i = 0; i < loci.size(); i++)
			{
				int num = atoi(loci.at(i).c_str());
				positions.push_back(num);
			}
			sort(positions.begin(), positions.end());
			positions.erase(unique(positions.begin(), positions.end()), positions.end());
			vector<string> Fpositions;
			for (int i : positions)
			{
				Fpositions.push_back(to_string(i));
			}
			vector<vector<string>> locus = rows;
			vector<string> temp1;
			for (unsigned int i = 0; i < rows.size(); i++)
			{
				temp1 = (rows[i]);
				unsigned int sz1 = temp1.size();
				for (unsigned int j = 1; j < sz1 - 6; j++)
				{
					int count = -1;
					for (string k : Fpositions)
					{
						count++;
						if (rows[i][j] == k)
						{
							locus[i][j] = " ";
							locus[i][j] = to_string(count);
							break;
						}
					}
				}
			}

			
			vector<vector<vector<string>>> partitions; //partitions the partial haplotypes from Multi_locus_trajectories.out into groups with matching loci numbers
			for (unsigned int k = 1; k <= Fpositions.size(); k++)
			{
				for (unsigned int i = 0; i < Fpositions.size(); i++)
				{
					vector<string> pick;
					vector<vector<string>> partialPar;
					for (unsigned int j = 0; j < locus.size(); j++)
					{
						pick = locus[j];
						string pos_num = to_string(k);
						unsigned int check_match = 0;
						if (pick[0] == pos_num)
						{
							for (unsigned int m = 0; m < k; m++)
							{
								if (pick[m + 1] == to_string(i + m))
								{
									check_match++;
								}
							}
							if (check_match == k)
							{
								partialPar.push_back(pick);
							}
						}
					}
					if (!partialPar.empty())
					{
						partitions.push_back(partialPar);
					}
				}
			}
			
			vector<string> tampa;//Here we sort each element of a partition set from largest to smallest with respect to the total number of reads
			for (unsigned int i = 0; i < (partitions).size(); i++)
			{
				double SIZE0 = (partitions[i]).size();
				for (unsigned int j = 0; j < SIZE0 - 1; j++)
				{
					double SIZE1 = (partitions[i][j]).size();
					double NUM1 = atoi((partitions[i][j]).at(SIZE1 - 1).c_str());
					double SIZE2 = (partitions[i][j + 1]).size();
					double NUM2 = atoi((partitions[i][j + 1]).at(SIZE2 - 1).c_str());
					if (NUM2 > NUM1)
					{
						tampa = partitions[i][j];
						partitions[i][j] = partitions[i][j + 1];
						partitions[i][j + 1] = tampa;
					}
				}
			}

			//Uncomment the following lines if you want to see how the jth member of the ith partial haplotype set looks like
			for (unsigned int i = 0; i<partitions.size(); i++)
			{
				for (unsigned int j = 0; j<(partitions[i]).size(); j++)
				{
					for (unsigned int k = 0; k<(partitions[i][j]).size(); k++)
					{
						//cout << partitions[i][j][k] << "  ";
					}
					//cout << endl << "--------------------------------------" << endl;
				}
				//cout << "**************************************" << endl;
			}
			
			
			vector<vector<vector<int>>> CONTRIBS;
			for (unsigned int i = 0; i < partitions.size(); i++)
			{
				vector<vector<int>> contribs;
				for (unsigned int j = 0; j < (partitions[i]).size(); j++)
				{
					vector<int> Contribs;
					vector<string> tempa;
					tempa = partitions[i][j];
					unsigned int sz = tempa.size();
					int num = sz - 6;
					for (unsigned int s = 0; s < candidates.size(); s++)
					{
						int counter = 0;
						for (int k = 1; k < num; k++)
						{
							int tem = atoi((tempa[k]).c_str());
							if ((candidates[s])[tem] == (tempa[num])[k - 1])
							{
								counter++;
							}
						}
						int dum = atoi((tempa[0]).c_str());
						if (counter == dum)
						{
							Contribs.push_back(s);
						}
					}
					contribs.push_back(Contribs);
				}
				CONTRIBS.push_back(contribs);
			}

			gsl_rng_env_setup();
			gsl_rng *rgen = gsl_rng_alloc(gsl_rng_taus);
			int seed = (int)time(NULL);
			gsl_rng_set(rgen, seed);
			
			vector<double> init_freqs1;//set the initial frequency of candidate haplotypes for the first timepoint to be a randomly distributed number between 0 and 1
			for (unsigned int i = 0; i < candidates.size() + 1; i++)
			{
				init_freqs1.push_back(gsl_rng_uniform(rgen));			
			}
			
			NormaliseFreqs(init_freqs1);
			vector<double> init_freqs1_store = init_freqs1;
			
			vector<double> init_freqs2;//set the initial frequency of candidate haplotypes for the second timepoint to be a randomly distributed number between 0 and 1
			for (unsigned int i = 0; i < candidates.size() + 1; i++)
			{
				init_freqs2.push_back(gsl_rng_uniform(rgen));			
			}
				
			NormaliseFreqs(init_freqs2);
				

			OptimumFreqs(0, candidates, init_freqs1, CONTRIBS, partitions, c);
			OptimumFreqs(1, candidates, init_freqs2, CONTRIBS, partitions, c);
			
			//L_store = L_store1 + L_store2;
			//cout << endl << L_store1 << " + " << L_store2 << " = " << L_store << endl;
			//cout << L_store << endl;
			
			/*for (unsigned int i = 0; i < f_after.size(); i++)
			{
				cout << init_freqs2[i]  << "\t" << f_after[i] << endl << init_freqs1[i] << "\t" << f_before[i] << endl;
			}
			
			cout << endl << "****************************************" << endl;*/
			
			for (unsigned int i = 0; i < candidates.size(); i++)
			{
				if (init_freqs1[i] < 1e-7)
				{
					init_freqs1[i] = 0;
				}
				if (init_freqs2[i] < 1e-7)
				{
					init_freqs2[i] = 0;
				}
			}
				
			chdir(argv[3]);
			ofstream myfile;
			myfile.open(argv[4]);
			for (unsigned int i = 0; i < f_after.size(); i++)
			{
				myfile << i << "\t" << init_freqs1[i] << "\t" << init_freqs2[i] << endl;
			}
			myfile.close();

			return 0;
			}
	}
	
	else
	{
		//cout << "No SimulatedData_Mahan_Gene.txt file found in " << xStar << endl;
		return false;
	}
}
