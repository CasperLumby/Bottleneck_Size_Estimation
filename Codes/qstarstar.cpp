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
		vector<double> qCopy = q; //copy of the input frequencies, freq
		qCopy.erase(remove(qCopy.begin(), qCopy.end() - 1, q[candidates_size]));//erase the unkknown haplotype qx out of the copy of the frequency list, qCopy. 
		unsigned int count_no_contrib = 0; //counts the number of haplotypes with no contribution to a partial haplotype read
		double no_contrib_reads = 0; //partial haplotype reads with no matching haplotype that contributes
		vector<double> no_contrib_haps; //list of all haplotype frequencies that do not contribute
		vector<double> contrib_haps; //list of all haplotype frequencies that contribute
		double contrib_reads = 0; //partial haplotype reads with, at least, one matching haplotype that contributes
		unsigned int count_contrib = 0; //counts the number of haplotypes contributing to a partial haplotype read
		double sum_nn = 0; //adds the number of reads
		//the following loop goes over all the haplotypes that contribute to the reads of a given partial haplotype set
		for (unsigned int j = 0; j<(CONTRIBS[i]).size(); j++) 
		{
			double sum = 0;
			int count_if_any_contrib = 0; //checks, for each partial haplotype, if any haplotype contributes
			for (unsigned int k = 0; k<(CONTRIBS[i][j]).size(); k++)//this loop adds the frequency of all the haplotypes that share the same partial haplotype of interest
			{
				//the if clause below insures that we add the frequency of all haplotypes that contribute to the ith partial haplotype set (CONTRIBS[i][j][k] = ith partial haplotype set, jth element, kth contributing haplotype) 
				sum += q[(CONTRIBS[i][j][k])];
				count_if_any_contrib++;
				count_contrib++;
				contrib_haps.push_back(q[(CONTRIBS[i][j][k])]);
			}
			//if none of the candidate haplotypes match with the partial haplotype of interest (i.e. do not 'contribute' to that partial haplotype set), the corresponding reads for that partial haplotype goes to the unknown haplotype qx)
			if ((CONTRIBS[i][j]).empty())
			{
				count_no_contrib++;
				double SIZE1_1 = (partitions[i][j]).size();
				no_contrib_reads = atoi((partitions[i][j]).at(SIZE1_1 - timepoint).c_str());
				sum_nn += no_contrib_reads;
				no_contrib_haps.push_back(no_contrib_reads);
			}
			//if, at least, there were one candidate haplotype that contributed to the ith partial haplotype set, uptade nn and inf
			if (count_if_any_contrib != 0)
			{
				inf.push_back(sum);
				double SIZE = (partitions[i][j]).size();
				contrib_reads = atoi((partitions[i][j]).at(SIZE - timepoint).c_str());
				nn.push_back(contrib_reads);
			}
		}
		//if the candidate haplotypes covered all the partial haplotypes in the ith set, then there are zero reads left to be attributed to qx
		if (count_no_contrib == 0 && count_contrib == qCopy.size())
		{
			inf.push_back(q[candidates_size]);
			nn.push_back(0.0);
		}
		unsigned int check2 = 0; 
		
		//if some (but not all) of the candidate haplotypes covered the entire partial haplotype set i, the remaining candidates would correspond to zero reads (so as the qx haplotype)
		if (count_no_contrib == 0 && count_contrib != qCopy.size())
		{
			for (unsigned int m = 0; m < contrib_haps.size(); m++)//this loop finds all the haplotypes that do not contribute and removes their frequencies from qCopy
			{
				qCopy.erase(remove(qCopy.begin(), qCopy.end() - 1, contrib_haps[m]));
			}
			double temp_sum = 0;
			for (unsigned int g = 0; g < qCopy.size(); g++)
			{
				temp_sum += qCopy[g];
			}
			inf.push_back(temp_sum);
			nn.push_back(0.0);
			inf.push_back(q[candidates_size]);
			nn.push_back(0.0);
			check2++;
		}
		
		//if some (but not all) of the candidate haplotypes did not contribute to the ith set, then they correspond to zero reads and the remaining unidentified partial haplotypes go to qx
		unsigned int check3 = 0;
		if (count_no_contrib != 0 && count_contrib != qCopy.size() && count_contrib != 0 && check2 == 0)
		{
			for (unsigned int m = 0; m < contrib_haps.size(); m++)
			{
				qCopy.erase(remove(qCopy.begin(), qCopy.end() - 1, contrib_haps[m]));
			}
			double temp_sum = 0;
			for (unsigned int g = 0; g < qCopy.size(); g++)
			{
				temp_sum += qCopy[g];
			}
			double total = accumulate(no_contrib_haps.begin(), no_contrib_haps.end(), 0);
			inf.push_back(temp_sum);
			nn.push_back(0.0);
			inf.push_back(q[candidates_size]);
			nn.push_back(total);
			check3++;
		}
		//if non of the candidate haplotypes contributed to ith set, then they all correspond to zero reads and everything else would be classified as qx
		if (count_no_contrib != 0 && count_contrib != qCopy.size() && count_contrib == 0 && check2 == 0 && check3 == 0)
		{
			double temp_sum = 0;
			for (unsigned int g = 0; g < qCopy.size(); g++)
			{
				temp_sum += qCopy[g];
			}
			double total = accumulate(no_contrib_haps.begin(), no_contrib_haps.end(), 0);
			inf.push_back(q[candidates_size]);
			nn.push_back(total);
			inf.push_back(temp_sum);
			nn.push_back(0.0);
		}

		//if all the candidate haplotypes covered some of the partial haplotypes, then all the remaining reads correspond to qx.
		if (count_no_contrib != 0 && count_contrib == qCopy.size() && check2 == 0 && check3 == 0)
		{
			inf.push_back(q[candidates_size]);
			nn.push_back(sum_nn);
		}
		ULTIMA_THULE += Likelihood(c, nn, inf);
		qCopy = q; //resetting the copy of the original frequencies
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

int main(int argc, char* argv[])
{
	
	if(argc < 6 || argc > 6) { 
		cout << "Error in the input format: " << endl << "Please use the following input order:" << endl;
		cout << "(1) Directory to simulated reads, i.e. x*" << endl;
		cout << "(2) Directory to reconstructed haplotypes, i.e. q*" << endl;
		cout << "(3) Directory to store re-inferred haplotype frequencies q**" << endl;
		cout << "(4) File name for the re-inferred haplotypes for a particular segment $t" << endl;
		cout << "(5) Inferred noise parameter C" << endl;
		return false;
	}
	
	else
	{
		cout << "qstarstar is now initiated..." << endl;
	}
	
	string xStar = argv[1];
	string reconstructed_haplotypes = argv[2];
	double c = atoi(argv[5]); 
	
	gsl_rng_env_setup();
	gsl_rng *rgen = gsl_rng_alloc(gsl_rng_taus);
	int seed = (int)time(NULL);
	gsl_rng_set(rgen, seed);
	
	if (fileExists(xStar) == true) 
	{
		if (fileExists(reconstructed_haplotypes) == true)
			{
			vector<string> table;
			ifstream test;
			test.open(reconstructed_haplotypes);
			while (!test.eof())
			{
				string line;
				getline(test, line);
				table.push_back(line);
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
			
			vector<double> f_before;
			vector<double> f_after;
			vector<string> candidates; //the list of all haplotypes (already constructed from q*)
			for (unsigned int i = 0; i < rows.size(); i++)
			{
				unsigned int sz = (rows[i]).size();
				candidates.push_back(rows[i][sz-3]);
				//string temp_1 = rows[i][sz-2];
				double num_1 = atof((rows[i][sz-2]).c_str());				
				f_before.push_back(num_1);
				//string temp_2 = rows[i][sz-1];
				double num_2 = atof((rows[i][sz-1]).c_str());
				f_after.push_back(num_2);
			}
			
			vector<string> table_2;
			ifstream reads;
			reads.open(xStar);
			while (!reads.eof())//read the content of Multi_locus_trajectories.out and store it in table
			{
				string line_2;
				getline(reads, line_2);
				table_2.push_back(line_2);
			}
			reads.close();
			
			vector < vector<string> > rows_2; //save each row of Multi_locus_trajectories.out
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
			
			int reads_total = 0; //calculating the total number of reads for all the partial haplotypes -- this quantity is required for BIC calculations 
			for (unsigned int i=0; i<rows_2.size(); i++)
			{
				unsigned int sz = (rows_2[i]).size();
				string reads_before_s = rows_2[i][sz-1];
				string reads_after_s = rows_2[i][sz-3];
				int reads_before = atoi(reads_before_s.c_str());
				int reads_after = atoi(reads_after_s.c_str());
				reads_total += reads_before + reads_after;
			}
			
			vector<int> positions; //collect the loci numbers from Multi_locus_trajectories.out
			for (unsigned int i = 0; i < rows_2.size(); i++)
			{
				unsigned int numLoci = atoi(rows_2[i][0].c_str()); //get the number of loci for the ith row
				for (unsigned int j = 1; j < 1 + numLoci; j++) //loop over the loci in the file
				{
					positions.push_back(atoi((rows_2[i][j]).c_str()));
				}
			}
			
			//The following ~30 lines of code sorts the order of loci from 0 to the last variant locus and store them in positions_string 
			sort(positions.begin(), positions.end());
			positions.erase(unique(positions.begin(), positions.end()), positions.end()); //up to this point, positions contain a list of all unique loci numbers sorted from smallest to largest, e.g. 121 322 1123 5694
			vector<string> positions_string; //convert positions to a vector of strings to do comparison on characters
			for (int i : positions)
			{
				positions_string.push_back(to_string(i));
			}
			vector<vector<string>> temprows = rows_2;//we want temprows to include all the information about Multi_locus_trajectories.out BUT with loci numbers that are sorted from 0 to n-1 for each partial haplotype of size n
			for (unsigned int i = 0; i < rows_2.size(); i++)//loop over all the rows in Multi_locus_trajectories.out
			{
				unsigned int numLoci = atoi(rows_2[i][0].c_str()); //get the number of loci for the ith row
				for (unsigned int j = 1; j < 1 + numLoci; j++)//loop over the loci for each partial haplotype (i.e. each row in Multi_locus_trajectories.out)
				{
					int count = -1;
					for (string k : positions_string)
					{
						count++;
						if (rows_2[i][j] == k)
						{
							temprows[i][j] = " ";
							temprows[i][j] = to_string(count); //what this does is that it would swap the actual loci number from SAMFIRE and re-label it as a sorted integer between [0,n-1] 
							break;
						}
					}
				}
			}
			//the updated rows of temprows are identical to original SAMFIRE Multi_locus_trajectories.out except that the loci are numbered from 0 to n-1

			
			vector<vector<vector<string>>> partitions; //partitions the partial haplotypes from Multi_locus_trajectories.out into groups with matching loci numbers
			for (unsigned int k = 1; k <= positions_string.size(); k++)//loop from 1 to n-1
			{
				for (unsigned int i = 0; i < positions_string.size(); i++)//loop over from 0 to n-1
				{
					vector<string> pick_row;
					vector<vector<string>> partialPar;
					for (unsigned int j = 0; j < temprows.size(); j++)//loop over all the rows in Multi_locus_trajectories.out
					{
						pick_row = temprows[j];
						string pos_num = to_string(k);
						unsigned int check_match = 0;
						if (pick_row[0] == pos_num)//if partial haplotype j has k number of loci
						{
							for (unsigned int m = 0; m < k; m++)//go over all of loci in partial haplotype j and see if it contains a unique sequence of loci
							{
								if (pick_row[m + 1] == to_string(i + m))
								{
									check_match++;
								}
							}
							if (check_match == k)//if there is a  perfect match, then it is an element of a given partition set
							{
								partialPar.push_back(pick_row);
							}
						}
					}
					if (!partialPar.empty())
					{
						partitions.push_back(partialPar);
					}
				}
			}
			
			
			bool changeOccurred = true;
			while(changeOccurred == true) 
			{
				changeOccurred = false;
				vector<string> sort_by_recipeint; //Here we sort each element of a partition set from largest to smallest with respect to the total number of reads in the recipient 
				// -- the pattern of sorting should not change if you do it with the donor reads (because it is always the case that a partial haplotype that has the 
				// highest number of reads in the donor would have the highest number of reads in the recipient)
				for (unsigned int i = 0; i < (partitions).size(); i++) //go over every element of the partition (defined above)
				{
					double SIZE0 = (partitions[i]).size();
					for (unsigned int j = 0; j < SIZE0 - 1; j++)// loop over the elements and sort partial haplotypes based on reads in the recipient
					{
						double SIZE1 = (partitions[i][j]).size();//select the jth partial haplotype
						double NUM1 = atoi((partitions[i][j]).at(SIZE1 - 1).c_str());//SIZE1 - 1 corresponds to the the number of reads in the recipient of every partial haplotype in Multi_locus_trajectories.out
						double SIZE2 = (partitions[i][j + 1]).size();//select the (j+1)th partial haplotype
						double NUM2 = atoi((partitions[i][j + 1]).at(SIZE2 - 1).c_str());
						if (NUM2 > NUM1)
						{
							sort_by_recipeint = partitions[i][j];
							partitions[i][j] = partitions[i][j + 1];
							partitions[i][j + 1] = sort_by_recipeint;
							changeOccurred = true;
						}
					}
				}
			}



			//Uncomment the following lines if you want to see how the jth element of the ith partial haplotype set looks like
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
			
			vector<vector<vector<int>>> CONTRIBS; //find the contribution of the candidate set given this change at a randomly selected locus (or loci) of one haplotype
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
			
			chdir(argv[3]); //directory to store re-inferred haplotype frequencies q**
			ofstream myfile;
			myfile.open(argv[4]); //file name for the re-inferred haplotypes for a particular segment $t
			for (unsigned int i = 0; i < f_after.size(); i++)
			{
				myfile << i << "\t" << init_freqs1[i] << "\t" << init_freqs2[i] << endl;
			}
			myfile.close();

			return 0;
			}
		else
		{
			//cout << "No outcome_1.txt file found in " << reconstructed_haplotypes << endl;
			return false;
		}
	}
	
	else
	{
		//cout << "No SimulatedData_Mahan_Gene.txt file found in " << xStar << endl;
		return false;
	}
}