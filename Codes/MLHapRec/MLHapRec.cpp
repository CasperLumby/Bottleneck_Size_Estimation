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
double OptimumFreqs(double timepoint, vector<string> &candidates, vector<double> &init_freqs, vector<vector<vector<int>>> &CONTRIBS, vector<vector<vector<string>>> &partitions, double c, int temp_seed)
{
	gsl_rng_env_setup();
	gsl_rng *rgen = gsl_rng_alloc(gsl_rng_taus);
	//int seed = (int)time(NULL);
	gsl_rng_set(rgen, temp_seed);
	double init_size = init_freqs.size();
	vector<double> init_freqs_store = init_freqs;
	double L = -1e7; //initial (extremely low) value for the likelihood
	unsigned int check = 1;
	unsigned int max_it = 8e3; //maximum number of attempts for optimising the frequencies
	double changex = 1e-2; //magnitude of the incremental (random) change in frequency
	double L_store1 = -1e7; //improved likelihood after frequency optimisation

	for (unsigned int it = 0; it < max_it; it++)
	{
		double s = gsl_rng_uniform_int(rgen,1000000)+1;
		s = s/(1000000);//a random number between zero and 1 discretised for 1 part per million 
		double r = 2 * (s - 0.5);//random direction between -1 to +1
		int j = floor(init_freqs.size()*s);
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
		if (check >= 50)//50 is found heuristically to be sufficient (from testing multiple simulations)
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
	//first check if the input format is right
	if(argc < 5 || argc > 5) { 
		cout << "Error in the input format: " << endl << "Please provide Multi_locus_trajectories.out file and noise parameter C" << endl;
		cout << "If you do not wish to provide a seed number or number of attempts, put 0 (or two 0s separated by a space if you do not wish to provide either) after inputting your file and noise parameter accordingly.";
		return false;
	}
	
	//taking in the input: first is the Multi_locus_trajectories.out file, second is the C-parameter for noise, third is the seed number, and forth is the number of attempts
	string multi_locus_file = argv[1];
	double c = atoi(argv[2]); 
	int seed = atoi(argv[3]); 
	unsigned int attempts_max = atoi(argv[4]); 

	//if only the ‘Multi_locus_trajectories.out’ and C-value is provided, MLHapRec would start with no initial seed number or a priori set of haplotype and reconstruct everything from the scratch.
	if(argc == 5) {
		if (seed != 0 && attempts_max == 0)
		{
			attempts_max = 20;
			cout << "You have inputted ‘Multi_locus_trajectories.out’ with a C-value = " << c <<", and seed number = " << seed << ".";
			cout << " A default number of attemps = " << attempts_max << " is used.";
			cout << endl << "MLHapRec is now initiated..." << endl;
		}
		else if (seed == 0 && attempts_max != 0)
		{
			seed = 1;
			cout << "You have inputted ‘Multi_locus_trajectories.out’ with a C-value = " << c <<", and number of attempts = " << attempts_max << ".";
			cout << " A default seed number = " << seed << " is used.";
			cout << endl << "MLHapRec is now initiated..." << endl;
		}
		else if (seed == 0 && attempts_max == 0)
		{
			seed = 1;
			attempts_max = 20;
			cout << "You have inputted ‘Multi_locus_trajectories.out’ with a C-value = " << c <<".";
			cout << " A default seed number = " << seed << " and number of attempts = " << attempts_max<< " are used.";
			cout << endl << "MLHapRec is now initiated..." << endl;
		}
		
		else
		{
			cout << "You have inputted ‘Multi_locus_trajectories.out’ with a C-value = " << c <<", seed number = " << seed << ", and number of attempts = " << attempts_max << ".";
			cout << endl << "MLHapRec is now initiated..." << endl;
		}
	}
	
	
	gsl_rng_env_setup();
	gsl_rng *rgen = gsl_rng_alloc(gsl_rng_taus);
	//seed = (int)time(NULL); uncomment this and it will make the MLHapRec reconstruction outcomes different each time you run the code.
	gsl_rng_set(rgen, seed); // this is the seed provided by the used. If no seed is provided, then seed=1 is used.
	
	freopen("update.txt","w",stdout); //print out the progress of the optimisation process in a file called update.txt
	if (fileExists(multi_locus_file) == true) //check whether the Multi_locus_trajectories.out file exists
	{
		vector<string> table;
		ifstream reads;
		reads.open(multi_locus_file);
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
		
		int reads_total = 0; //calculating the total number of reads for all the partial haplotypes -- this quantity is required for BIC calculations 
		for (unsigned int i=0; i<rows.size(); i++)
		{
			unsigned int sz = rows[i].size();
			string reads_before_s = rows[i][sz-1];
			string reads_after_s = rows[i][sz-3];
			int reads_before = atoi(reads_before_s.c_str());
			int reads_after = atoi(reads_after_s.c_str());
			reads_total += reads_before + reads_after;
		}
		
		vector<int> positions; //collect the loci numbers from Multi_locus_trajectories.out
		for (unsigned int i = 0; i < rows.size(); i++)
		{
			unsigned int numLoci = atoi(rows[i][0].c_str()); //get the number of loci for the ith row
			for (unsigned int j = 1; j < 1 + numLoci; j++) //loop over the loci in the file
			{
				positions.push_back(atoi((rows[i][j]).c_str()));
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
		vector<vector<string>> temprows = rows;//we want temprows to include all the information about Multi_locus_trajectories.out BUT with loci numbers that are sorted from 0 to n-1 for each partial haplotype of size n
		for (unsigned int i = 0; i < rows.size(); i++)//loop over all the rows in Multi_locus_trajectories.out
		{
			unsigned int numLoci = atoi(rows[i][0].c_str()); //get the number of loci for the ith row
			for (unsigned int j = 1; j < 1 + numLoci; j++)//loop over the loci for each partial haplotype (i.e. each row in Multi_locus_trajectories.out)
			{
				int count = -1;
				for (string k : positions_string)
				{
					count++;
					if (rows[i][j] == k)
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
		
		vector<vector<string>> shuffle_partitions; //this stores all the information of every partial haplotype in one list which is later going to be used as a shuffling set for optimising a set of candidate haplotypes
		for (unsigned int i = 0; i < partitions.size(); i++)
		{
			for (unsigned int j = 0; j < (partitions[i]).size(); j++)
			{
				shuffle_partitions.push_back(partitions[i][j]);
			}
		}

		
		//our 'default' candidate haplotype has no information in it, i.e. it starts with unknown nucleotide X at each position
		string unknown;
		for (unsigned int i = 0; i<positions.size(); i++)
		{
			unknown.append("X");
		}
					
		vector<string> candidates;//a list of candidate (or known) haplotypes; The code optimises and appends new ones until BIC value becomes negative	
		
		int row_numbers = rows.size();//number of rows in Multi_locus_trajectories.out
		double L_new = -1e7;
		double L_store = -1e7;
		double L_preserve_1 = -1e7 - 1;
		double L_preserve_2 = -1e7 - 1;
		double loop = 0; //keeps the track of how big is the candidate haplotype vector (i.e. same as candidates.size())
		vector<double> all_likelihoods;
		vector<double> haplotypes;
		vector<double> freqs_before;
		vector<double> freqs_after;
		double BIC_check = 1;
		while (BIC_check > 0)
		{
			loop++;
			candidates.push_back(unknown); //every round, we add the 'default' haplotype (...XXXX...) into candidate list and optimise it based on information from short-read data
			unsigned int candidates_size = candidates.size();
			L_preserve_1 = L_preserve_2;
			vector<double> temp_likelihoods;
			vector<vector<double>> freqs_1;
			vector<vector<double>> freqs_2;
			vector<vector<string>> temp_candidates;
			int candidate_size1 = candidates.size();
			//The following loop determines, for a given set of candidate haplotypes, how many attempts we make to optimise the set such that BIC > 0
			//In the default case, the code makes 20 attempts to optimize candidates.size() haplotypes.
			for (unsigned int attempts = 0; attempts < attempts_max; attempts++)
			{
				L_new = L_preserve_2 - 100 / (2*candidate_size1); //try to explore the possible haplotype space locally, i.e. start optimising likelihood from a value close to the previous step of the optimisation with one fewer candidate haplotype
				L_store = -1e5;
				candidates[candidate_size1 - 1] = unknown; //the most recent haplotype is unknown
				double candidate_size = candidates.size();
				double THRESHOLD = attempts_max * loop * row_numbers; //depending on how big is the Multi_locus_trajectories.out file and how many haplotypes are already added as potential candidates, the number of attempts to re-shuffling the candidate haplotype set must be adjusted 
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
				vector<double> init_freqs2_store = init_freqs2;
				
				vector<vector<vector<int>>> CONTRIBS; //finds which (if any) candidate haplotype matches with the jth element of the ith partial haplotype set in Multi_locus_trajectories.out (i.e. which haplotypes 'contribute' to the jth partial haplotype of the ith set)
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


				for (unsigned int s = 0; s < 50*(positions_string.size()); s++)//set the 'initial conditions' for the most recent (unknown) candidate haplotype using pieces from the partial haplotype set until there are no undetermined loci (X) left
				{
					double randomX = floor(shuffle_partitions.size()*gsl_rng_uniform(rgen));
					string replaceX;
					vector<double> posX;
					double loci_numberX = atoi((shuffle_partitions[randomX][0]).c_str());
					for (unsigned int i = 0; i < loci_numberX; i++)
					{
						double shuf_sizeX = shuffle_partitions[randomX].size();
						replaceX.push_back(shuffle_partitions[randomX][shuf_sizeX - 6][i]);
						posX.push_back(atoi((shuffle_partitions[randomX][i + 1]).c_str()));
					}
					string sub_hapX = candidates[candidate_size - 1].substr(posX[0], posX.size());
					candidates[candidate_size - 1].replace(posX[0], posX.size(), replaceX);
				}
				
				double L_store1;//maximum likelihood of the candidate haplotypes given the first timepoint (i.e. donor population)
				double L_store2;//maximum likelihood of the candidate haplotypes given the second timepoint (i.e. recipient population)
				double failure_counts = 0;//This defines the cut-off on when to stop making further attempts to optimize the haplotypes
				while (true)
				{
					double j = floor(candidates.size()*gsl_rng_uniform(rgen));//randomly select a candidate haplotype
					double random = floor(shuffle_partitions.size()*gsl_rng_uniform(rgen));//and replace a randomly selected position(s) with a partial haplotype that could exist at that position
					string REPLACE;
					vector<double> pos;
					double loci_number = atoi((shuffle_partitions[random][0]).c_str());
					for (unsigned int i = 0; i < loci_number; i++)
					{
						double shuf_size = shuffle_partitions[random].size();
						REPLACE.push_back(shuffle_partitions[random][shuf_size - 6][i]);
						pos.push_back(atoi((shuffle_partitions[random][i + 1]).c_str()));
					}
					string sub_hap = candidates[j].substr(pos[0], pos.size());
					candidates[j].replace(pos[0], pos.size(), REPLACE);

					vector<vector<vector<int>>> CONTRIBS; //find the contribution of the candidate set given this change at a randomly selected locus (or loci) of one haplotype
					for (unsigned int i = 0; i < partitions.size(); i++)
					{
						vector<vector<int>> contribs;
						for (unsigned int j = 0; j < (partitions[i]).size(); j++)
						{
							vector<int> Contribs;
							vector<string> tempa;
							tempa = partitions[i][j];
							int sz = tempa.size();
							int num = sz - 6;
							for (unsigned int s = 0; s < candidates.size(); s++)
							{
								unsigned int counter = 0;
								for (int k = 1; k < num; k++)
								{
									int tem = atoi((tempa[k]).c_str());
									if ((candidates[s])[tem] == (tempa[num])[k - 1])
									{
										counter++;
									}
								}
								unsigned int dum = atoi((tempa[0]).c_str());
								if (counter == dum)
								{
									Contribs.push_back(s);
								}
							}
							contribs.push_back(Contribs);
						}
						CONTRIBS.push_back(contribs);
					}

					L_store1 = OptimumFreqs(0, candidates, init_freqs1, CONTRIBS, partitions, c, seed);//optimise the frequencies of + calculate the log-likelihood of the candidate haplotype set with noise parameter c for the first timepoint (i.e. before transmission in the donor population) 
					L_store2 = OptimumFreqs(1, candidates, init_freqs2, CONTRIBS, partitions, c, seed);//optimise the frequencies of + calculate the log-likelihood of the candidate haplotype set with noise parameter c for the second timepoint (i.e. after transmission in the recipient population)
					L_store = L_store1 + L_store2;//the log-likelihood of the first and second timepoint are independent of each other

					if (L_store > L_new)//if this random change was beneficial (i.e. increased the likelihood of the candidate set)
					{
						//cout << endl << L_store1 << " + " << L_store2 << " = " << L_store << "  " << L_new << endl;
						L_new = L_store;
						//cout << endl << "---->IMPROVE!" << endl;
						failure_counts = 1;
						for (unsigned int i = 0; i < init_freqs1.size(); i++)
						{
							init_freqs1_store[i] = init_freqs1[i];
						}
						for (unsigned int i = 0; i < init_freqs2.size(); i++)
						{
							init_freqs2_store[i] = init_freqs2[i];
						}
					}
					else if (L_store <= L_new)//if it was a bad change, reverse the random change and try again until the failure_counts=THRESHOLD
					{
						//cout << endl << L_store1 << " + " << L_store2 << " = " << L_store << "  " << L_new << endl;
						candidates[j].replace(pos[0], pos.size(), sub_hap);
						//cout << endl << "---->FAIL!" << endl;
						failure_counts++;
					}
					if (failure_counts > THRESHOLD || shuffle_partitions.size() == 0)//if we tried 'enough' number of random changes (set by the THRESHOLD) and no improvement in the likelihood value was gained, it means that we have found the optimal set
					{
						break;
					}
				}
				
				
				//Uncomment the following lines if you like to see the output of each attempt in the optimisation process
				/*cout << "candidate haplotype(s) is (are):" << endl;
				for (unsigned int i = 0; i < candidates.size(); i++)
				{
				cout << candidates[i] << "   ";
				}
				cout << endl;
				cout << "likelihood found in attempt " << attempts << ": " << L_new << endl;
				for (unsigned int i = 0; i < init_freqs1.size(); i++)
				{
				cout << init_freqs1[i] << "  ";
				}
				cout << endl;
				for (unsigned int i = 0; i < init_freqs1.size(); i++)
				{
				cout << init_freqs2[i] << "  ";
				}
				cout << endl << endl;*/
				
				
				if (loop > 1)
				{
					int candidates_check = 0;
					for (unsigned int i = 0; i < candidate_size - 1; i++)
					{	
						if (candidates[candidate_size-1] != candidates[i]) // check to see if the most recently added haplotype is a new haplotype (i.e. not a repeated haplotype that already exists in the candidates list)
						{
							candidates_check++; //if there are no repeats, this should be the same as the total number of haplotypes - 1
						}
					}
					
					int bad_hap = 0;
					for (unsigned int i = 0; i < candidate_size; i++) //if the frequency of any haplotype is effectively zero before and after the transmission, then it is really not a good candidate and the code has got stuck in a local maximum
					{	
						if (init_freqs1_store[i] < 1e-9 && init_freqs2_store[i] < 1e-9)
						{
							bad_hap++;
						}
					}
					
					if (candidates_check == candidate_size - 1 && bad_hap == 0)//if there is no duplicate haplotype and their frequencies BOTH before and after tranmission is above 10^-9, optimisation is normal
					{
						temp_likelihoods.push_back(L_new);
						freqs_1.push_back(init_freqs1_store);
						freqs_2.push_back(init_freqs2_store);
						temp_candidates.push_back(candidates);
					}
					else if ((candidates_check != candidate_size - 1) || (candidates_check == candidate_size - 1 && bad_hap != 0))//otherwise, make another attempt to find a better set of haplotypes
					{
						attempts = attempts - 1;
					}
				}
				else if (loop == 1) //this is the first haplotype so it is fine.
				{
					temp_likelihoods.push_back(L_new);
					freqs_1.push_back(init_freqs1_store);
					freqs_2.push_back(init_freqs2_store);
					temp_candidates.push_back(candidates);
				}
			}
			
			//find which of the 20 attempts (default) had the highest likelihood value
			double max = temp_likelihoods[0];
			double max_it = 0;
			for (unsigned int ml = 1; ml < temp_likelihoods.size(); ml++)
			{
				if (temp_likelihoods[ml] > max)
				{
					max = temp_likelihoods[ml];
					max_it = ml;
				}
			}
			
			cout << "------------------------------------------------------------------" << endl;
			cout << "***ROUND " << candidates_size << "***" << endl;
			cout << "number of haplotypes are currently = " << candidates_size << endl;
			cout << "maximum log-likelihood value = " << temp_likelihoods[max_it] << endl << endl;
			cout << "below is the current list of candidate haplotypes:" << endl;
			cout << "<optimised haplotype> | <donor frequency> | <recipient frequency> | " << endl;
			for (unsigned int i = 0; i < temp_candidates[max_it].size(); i++)
			{
				cout << temp_candidates[max_it][i] << '\t' << freqs_1[max_it][i] << '\t' << freqs_2[max_it][i] << endl;
			}
			cout << "------------------------------------------------------------------" << endl;
			
			all_likelihoods.push_back(temp_likelihoods[max_it]);
			L_preserve_2 = temp_likelihoods[max_it];
			//the following line calculates the difference between the BIC of the current step and one step before; BIC = -2log(L) + klog(n)
			BIC_check = -((-2*L_preserve_2)+candidates_size*log (reads_total)) + ((-2*L_preserve_1)+(candidates_size-1)*log (reads_total));
			if (BIC_check > 0) // if Delta_BIC > 0 --> the current set of N haplotypes are better than the N-1 haplotypes in the previous step
			{
				ofstream myfile;
				myfile.open("raw_haplotypes.txt");
				for (unsigned int i = 0; i < temp_candidates[max_it].size(); i++)
				{
					myfile << i + 1 << "\t" << temp_candidates[max_it][i] << "\t" << freqs_1[max_it][i] << "\t" << freqs_2[max_it][i] << endl;
				}
				myfile.close();
				
				for (unsigned int i = 0; i < temp_candidates[max_it].size(); i++)
				{
					if (freqs_1[max_it][i] < 1e-7)
					{
						freqs_1[max_it][i] = 0;
					}
					if (freqs_2[max_it][i] < 1e-7)
					{
						freqs_2[max_it][i] = 0;
					}
				}
				
				ofstream myfile2;
				myfile2.open("outcome_1.txt");
				for (unsigned int i = 0; i < temp_candidates[max_it].size(); i++)
				{
					if (freqs_1[max_it][i] == 0 && freqs_2[max_it][i] > 0)
					{
						//This is a mutation!
					}
					if (freqs_1[max_it][i] != 0 && freqs_1[max_it][i] != freqs_2[max_it][i])
					{
						myfile2 << i + 1 << "\t" << temp_candidates[max_it][i] << "\t" << freqs_1[max_it][i] << "\t" << freqs_2[max_it][i] << endl;
					}
				}
				myfile2.close();
			}
			
			if (BIC_check < 0)
			{
				cout << "The optimisation process stops here." << endl;
				cout << "The BIC in ROUND " << candidates_size << " is negative" << endl;
			}
		}
		
		return 0;
	}
	
	else
	{
		cout << "No Multi_locus_trajectories.out file found in " << multi_locus_file << endl;
		return false;
	}
}