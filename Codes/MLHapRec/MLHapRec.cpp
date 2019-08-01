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
	//First check if the input format is right, if not, exit
	if(argc != 5) { 
		cout << "Error in the input format: " << endl << "Please provide Multi_locus_trajectories.out file and noise parameter C" << endl;
		cout << "If you do not wish to provide a seed number or number of attempts, put 0 (or two 0s separated by a space if you do not wish to provide either) after inputting your file and noise parameter accordingly.";
		return false;
	}
	
	//Get input: first is the Multi_locus_trajectories.out file, second is the (noise) C-parameter, third is the random number generator seed, and fourth is the number of attempts
	string multiLocusFilePath = argv[1];
	double c = atoi(argv[2]); 
	int seed = atoi(argv[3]); 
	unsigned int attempts_max = atoi(argv[4]); 

	//Define default parameters for seed and attempts_max if these are not provided by the user
	if (seed != 0 && attempts_max == 0) { //Seed not given

		attempts_max = 20;
		cout << "You have inputted ‘Multi_locus_trajectories.out’ with a C-value = " << c <<", and seed number = " << seed << ".";
		cout << " A default number of attemps = " << attempts_max << " is used.";
		cout << endl << "MLHapRec is now initiated..." << endl;
	
	} else if (seed == 0 && attempts_max != 0) { //Attempts_max not given

		seed = 1;
		cout << "You have inputted ‘Multi_locus_trajectories.out’ with a C-value = " << c <<", and number of attempts = " << attempts_max << ".";
		cout << " A default seed number = " << seed << " is used.";
		cout << endl << "MLHapRec is now initiated..." << endl;
	
	} else if (seed == 0 && attempts_max == 0) { //Both seed and attempts_max not given

		seed = 1;
		attempts_max = 20;
		cout << "You have inputted ‘Multi_locus_trajectories.out’ with a C-value = " << c <<".";
		cout << " A default seed number = " << seed << " and number of attempts = " << attempts_max<< " are used.";
		cout << endl << "MLHapRec is now initiated..." << endl;

	} else { //User defined everything

		cout << "You have inputted ‘Multi_locus_trajectories.out’ with a C-value = " << c <<", seed number = " << seed << ", and number of attempts = " << attempts_max << ".";
		cout << endl << "MLHapRec is now initiated..." << endl;
		
	}
	
	
	gsl_rng_env_setup();
	gsl_rng *rgen = gsl_rng_alloc(gsl_rng_taus);
	//seed = (int)time(NULL); uncomment this and it will make the MLHapRec reconstruction outcomes different each time you run the code.
	gsl_rng_set(rgen, seed); // this is the seed provided by the used. If no seed is provided, then seed=1 is used.
	
//	freopen("update.txt","w",stdout); //print out the progress of the optimisation process in a file called update.txt

	//Progress with inference only if multi-locus file exists, otherwise print error
	if (fileExists(multiLocusFilePath) == true) {

		/*
		 * Load multi-locus file. It looks like the following:
		 *
		 *
		 * 2 126 245 TT 2 0 1094 1 0 
		 * 2 126 245 AT 2 0 3086 1 5637 
		 * 2 126 245 TC 2 0 293 1 0 
		 * 2 126 245 AC 2 0 139 1 0 
		 * 1 126 T 2 0 18020 1 0 
		 * 1 126 A 2 0 32940 1 75400 
		 * 1 126 G 2 0 18 1 19 
		 * 1 126 C 2 0 16 1 0 
		 * 1 245 T 2 0 83172 1 141401 
		 * 1 245 C 2 0 14081 1 108 
		 * 1 1346 A 2 0 106412 1 129827 
		 * 1 1346 G 2 0 84 1 5226 
		 * 1 1346 T 2 0 0 1 15 
		 *
		 *
		 * i.e. one line for each partial haplotype entry.
		 */


		//Save the data as a vector of partial haplotypes
		vector<string> multiLocusData;
		ifstream multiLocusFile;
		multiLocusFile.open(multiLocusFilePath);
		while(!multiLocusFile.eof()) { //Read one line at a time

			string line;
			getline(multiLocusFile, line);
			multiLocusData.push_back(line);

		}
		if(multiLocusData.empty()) {

			cout << "Multi_locus_trajectories.out file was empty. Exiting.\n";
			exit(1);
		}
		multiLocusFile.close();


		//Now convert each entry in the q* data into a vector of strings, i.e. split by whitespace
		vector<vector<string> > multiLocusDataSplit;
		for(unsigned int i=0; i<multiLocusData.size(); i++) { //Loop over each entry

			vector<string> row;
			istringstream iss(multiLocusData[i]);
			string term;
			while (iss >> term)
			{
				row.push_back(term);
			}
			if (!row.empty())
			{
				//Check if row corresponds to a partial haplotype in the category X, which we don't include.
				//This may happen depending on how the SAMFIRE output is generated.
				int numLoci = atoi(row[0].c_str());
				if(row[numLoci+1] != "X") {
					multiLocusDataSplit.push_back(row);
				}
			}
		}
		
		
		//Next we calculcate the total number of reads across all partial haplotypes.
		//This quantity is used in the calculation of BIC statistics
		int reads_total = 0; 
		for (unsigned int i=0; i<multiLocusDataSplit.size(); i++) { //Loop over partial haplotype entries

			unsigned int sz = multiLocusDataSplit[i].size();
			string reads_before_string = multiLocusDataSplit[i][sz-3];
			string reads_after_string = multiLocusDataSplit[i][sz-1];
			int reads_before = atoi(reads_before_string.c_str());
			int reads_after = atoi(reads_after_string.c_str());
			reads_total += reads_before + reads_after;
		}
		
		//Next we collect the loci numbers from Multi_locus_trajectories.out
		vector<int> positions; //List of loci, not unique, i.e. repeated values may be found
		for (unsigned int i=0; i<multiLocusDataSplit.size(); i++) { //Loop over partial haplotype entries

			unsigned int numLoci = atoi(multiLocusDataSplit[i][0].c_str()); //get the number of loci for the ith row
			for (unsigned int j = 1; j < 1 + numLoci; j++) //loop over the loci in the file
			{
				positions.push_back(atoi((multiLocusDataSplit[i][j]).c_str()));
			}
		}
		
		//The following ~10 lines of code sorts the order of loci from 0 to the last variant locus and then store them in positions_string 
		//First sort the positions (duplicates exist)
		sort(positions.begin(), positions.end());

		//Next, remove duplicates and decrease size of positions to match unique values
		positions.erase(unique(positions.begin(), positions.end()), positions.end()); //positions now contain a list of all unique loci numbers sorted from smallest to largest, e.g. 121 322 1123 5694

		//Convert positions to a vector of strings, as this allows us to do direct comparisons with the multi-locus file (multiLocusDataSplit)
		vector<string> positions_string;
		for (int pos : positions){ 

			positions_string.push_back(to_string(pos));
		}


		/*
		 * Next we want to create a version of multiLocusDataSplit which has loci numbers, i.e. values from [0,n) where
		 * n is the number of loci in the system, rather than positions.
		 *
		 * E.g. 
		 *
		 * 2 126 245 TT 2 0 1094 1 0 
		 * 2 126 245 AT 2 0 3086 1 5637 
		 * 2 126 245 TC 2 0 293 1 0 
		 * 2 126 245 AC 2 0 139 1 0 
		 * 1 126 T 2 0 18020 1 0 
		 * 1 126 A 2 0 32940 1 75400 
		 * 1 126 G 2 0 18 1 19 
		 * 1 126 C 2 0 16 1 0 
		 * 1 245 T 2 0 83172 1 141401 
		 * 1 245 C 2 0 14081 1 108 
		 * 1 1346 A 2 0 106412 1 129827 
		 * 1 1346 G 2 0 84 1 5226 
		 * 1 1346 T 2 0 0 1 15 
		 *
		 * 
		 * becomes
		 *
		 * 2 0 1 TT 2 0 1094 1 0 
		 * 2 0 1 AT 2 0 3086 1 5637 
		 * 2 0 1 TC 2 0 293 1 0 
		 * 2 0 1 AC 2 0 139 1 3 
		 * 1 0 T 1 0 18020 1 5 
		 * 1 0 A 1 0 32940 1 75400 
		 * 1 0 G 1 0 18 1 19 
		 * 1 0 C 1 0 16 1 1 
		 * 1 1 T 1 0 83172 1 141401 
		 * 1 1 C 1 0 14081 1 108 
		 * 1 2 A 1 0 106412 1 129827 
		 * 1 2 G 1 0 84 1 5226 
		 *
		 * as n=2.
		 */
		vector<vector<string> > multiLocusDataSplitLoci = multiLocusDataSplit;

		for (unsigned int i = 0; i < multiLocusDataSplit.size(); i++) {//Loop over partial haplotype entries

			int numLoci = atoi(multiLocusDataSplit[i][0].c_str()); //Get the number of loci for the ith entry
			for (int j = 1; j < 1 + numLoci; j++) { //Go through each locus in the partial haplotype


				for(unsigned int posIndex=0; posIndex<positions_string.size(); posIndex++) { //Loop over all n loci

					//Check if pos at posIndex matches the current position in entry i, locus j
					if(multiLocusDataSplit[i][j] == positions_string[posIndex]) {

						//Update the [i][j] entry to be the locus value (i.e. an integer on [0,n)) rather than the (nucleotide) position
						multiLocusDataSplitLoci[i][j] = to_string(posIndex);
						break;
					}
				}
			}
		}


		//Next we wish to partition the multi-locus data into sets of partial haplotypes covering the same loci.
		
		vector<vector<vector<string>>> partitions; //partitions the partial haplotypes from Multi_locus_trajectories.out into groups with matching loci numbers
		for (unsigned int k = 1; k <= positions_string.size(); k++) { //Loop over k, the number of loci in the partial haplotype, i.e. k is on [1,n]

			for (unsigned int i = 0; i < positions_string.size(); i++) { //Loop over the starting positions, i.e. i is on [0,n)

				vector<vector<string>> partition; //I.e. a vector of partial haplotype entries that share a common set of loci

				//Go through all the partial haplotype entries and see if they match the current loci combination, if so add them to the partition
				for (unsigned int j = 0; j < multiLocusDataSplitLoci.size(); j++) { 

					vector<string> partialHapEntry = multiLocusDataSplitLoci[j];
					string numLoci = to_string(k); //Number of loci covered by partial haplotypes currently considered
					if (partialHapEntry[0] == numLoci) { //Check if partial haplotype j has k number of loci


						//Loop over the specific loci in the partial hap j to see if they correspond to the loci currently considered, i.e. the loci {i,i+1,...,i+k}
						bool partialHapMatch = true; //Assume the partial hap matches the current set of loci, then check
						for (unsigned int m = 0; m < k; m++) {

							if (partialHapEntry[m + 1] != to_string(i + m)) {

								partialHapMatch = false;
								break;
							}
						}

						if (partialHapMatch == true) {  //Partial hap matches the set of loci currently explored, so add it to partition

							partition.push_back(partialHapEntry);
						}
					}
				}
				if (!partition.empty()) { //If this partition contains data, add it to partitions

					partitions.push_back(partition);
				}
			}
		}
		

		//Here we sort each element of a partition set from largest to smallest with respect to the total number of reads in the recipient.
		bool changeOccurred = true;
		while(changeOccurred == true) { //Keep sorting until no change occur 

			changeOccurred = false;
			for (unsigned int i = 0; i < (partitions).size(); i++) { //Loop over all partitions

				//Suppose k haplotypes in partition i. We here loop over j defined in [0,k-1].
				//We then compare the number of reads in the recipient for haplotypes j and j+1
				//and swap the haplotypes if hap j+1 has more reads than hap j.
				for (unsigned int j=0; j<partitions[i].size()-1; j++) { 

					//Partial haplotype j
					int lenJ = partitions[i][j].size(); //Number of entries in haplotype j, e.g. if have 2 0 1 AT 2 0 3086 1 5637, then lenJ = 9.
					int numReadsJ = atoi(partitions[i][j][lenJ-1].c_str()); //Number of reads in recipient for partial hap j
	
					//Partial haplotype j+1
					int lenJplusOne = partitions[i][j+1].size(); //Number of entries in haplotype j+1, e.g. if have 2 0 1 AT 2 0 3086 1 5637, then lenJplusOne = 9.
					int numReadsJplusOne = atoi(partitions[i][j+1][lenJplusOne-1].c_str()); //Number of reads in recipient for partial hap j+1				

					//Swap haplotypes if number of reads for j+1 > number of reads for j
					if (numReadsJplusOne > numReadsJ) {

						vector<string> haplotypeJtemporary = partitions[i][j];
						partitions[i][j] = partitions[i][j + 1];
						partitions[i][j + 1] = haplotypeJtemporary;
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
					cout << partitions[i][j][k] << "  ";
				}
				cout << endl << "--------------------------------------" << endl;
			}
			cout << "**************************************" << endl;
		}
		
		/*
		 * We now move on to haplotype reconstruction
		 */ 

		//In the process of reconstructing haplotypes, we generate candidate haplotypes from a starting
		//point of an uninformative haplotype, here represented by the unknown nucleotide X being found
		//at  every position.
		string unknown;
		for (unsigned int i = 0; i<positions.size(); i++)
		{
			unknown.append("X");
		}
		
		//We create a vector of candidate haplotypes, which we optimise/extend progressively.
		//The code optimises and appends new haplotypes until DeltaBIC (difference between two
		//models) is < 0.		
		vector<string> candidates;
		
		int numPartialHaps = (int) multiLocusDataSplit.size(); //Number of rows in Multi_locus_trajectories.out, i.e. number of partial haplotypes

		//Likelihoods used in optimisation process
		double L_best = -1e7;
		double L_new = -1e7;
		double L_preserve_1 = -1e7 - 1;
		double L_preserve_2 = -1e7 - 1;

		
		double numCandidateHaplotypes = 0; //keeps the track of how big is the candidate haplotype vector (i.e. same as candidates.size())
		vector<double> all_likelihoods;
		vector<double> haplotypes;
		vector<double> freqs_before;
		vector<double> freqs_after;
		double DeltaBIC = 1; //This is actually delta BIC. Larger than 0 indicates that k haplotypes is better than k-1 haplotypes. First assume this is the case, then check.
		while (DeltaBIC > 0) {

			numCandidateHaplotypes++; //New round, so we attempt to add another haplotype
			candidates.push_back(unknown); //Every round, we add the 'default' haplotype (...XXXX...) into candidate list and optimise it based on information from short-read data

			//Save likelihood 2 as likelihood 1?
			L_preserve_1 = L_preserve_2;

			vector<double> temp_likelihoods;
			vector<vector<double>> freqs_1;
			vector<vector<double>> freqs_2;
			vector<vector<string>> temp_candidates;


			//Given a number k of candidate haplotypes, we attempt to optimise these, i.e. reconstruct the haplotypes and their frequencies.
			//The goal is to do this such that DeltaBIC > 0, which indicates that  k haplotypes fit the data better than k-1 haplotypes.
			//However, in order to limit the computational overhead, we here only attempt to optimise the halotypes for a total of attempts_max
			//times, after which we give up (and accept the k-1 number of haplotype solution). As default, we have attempts_max = 20, but the
			//user may alter this number with a command line argument.
			for (unsigned int attempts = 0; attempts < attempts_max; attempts++) {


				//When it comes to optimising the k haplotypes/frequencies we want to ensure that our optimisation doesn't end up too far
				//from our optimisation for k-1 haplotypes. To this end, we restrict the likelihood for the k haplotypes to be in the range
				//of the k-1 haplotype likelihood, rather than having an unconstrained starting point. This ensures a faster and more local
				//optimisation.
				//The constrain is heuristic, here set to L_k = L_{k-1} - 100/(2*k), e.g. if k=4 and L_{k-1} = -5000, then L_k starts at -5000-100/8 = -5012.25.  
				L_best = L_preserve_2 - 100 / (2*numCandidateHaplotypes); 
				L_new = -1e5; //What is this?

				//All the first k-1 haplotypes remain unchanged, but the kth haplotype is set to the unknown (XXXXXX) haplotype, as we start a new optimisation attempt
				candidates[numCandidateHaplotypes - 1] = unknown;



				//When it comes to generating new candidate haplotypes, we do this by a process of reshuffling partial haplotypes.
				//The following variable, THRESHOLD, defines the number of attempts made at reshuffling partial haplotypes. This number
				//necessarily needs to be larger, the more data (partial haplotypes) we have available, in order to ensure sufficient
				//sampling. In general, it depends on:
				//
				// - the number of attempts we make at constructing haplotypes (why??)
				// - the number of candidate haplotypes currently explored (why??)
				// - the number of partial haplotypes to reshuffle
				double THRESHOLD = attempts_max * numCandidateHaplotypes * numPartialHaps; //depending on how big is the Multi_locus_trajectories.out file and how many haplotypes are already added as potential candidates, the number of attempts to re-shuffling the candidate haplotype set must be adjusted 

				//Here we initialise the frequencies for the haplotypes for the first time point.
				//These are initialised to random numbers between 0 and 1, and then normalised.
				vector<double> init_freqs1;
				for (unsigned int i = 0; i < candidates.size() + 1; i++)
				{
				
					init_freqs1.push_back(gsl_rng_uniform(rgen));			
				}
				NormaliseFreqs(init_freqs1);
				vector<double> init_freqs1_store = init_freqs1; //Here we store the initial frequencies, for some reason. Check why?
	
				//Here we initialise the frequencies for the haplotypes for the second time point.
				//These are initialised to random numbers between 0 and 1, and then normalised.			
				vector<double> init_freqs2;
				for (unsigned int i = 0; i < candidates.size() + 1; i++)
				{
				
					init_freqs2.push_back(gsl_rng_uniform(rgen));			
				}
				NormaliseFreqs(init_freqs2);
				vector<double> init_freqs2_store = init_freqs2;
				

				
				//The quantity contribsPartitions represents information about which of the candidate haplotypes contributes to which partial haplotypes.
				//In particular, in the below we identify which haplotypes contribute to the jth haplotype in the ith haplotype set.
				vector<vector<vector<int>>> contribsPartitions; 
				for (unsigned int i = 0; i < partitions.size(); i++) { //Loop over sets of partial haplotypes, denoted by i

					vector<vector<int>> contribsPartialHapSet; //Define contributions for this set
					for (unsigned int j = 0; j < (partitions[i]).size(); j++) { //Loop over partial haplotypes in set i

						vector<int> contribsPartialHap; //Define contributions for this haplotype (set i, partial hap j)
						vector<string> currentPartialHap = partitions[i][j];
						int numLociCurrentPartialHap = atoi(currentPartialHap[0].c_str());

						//Loop over the candidate haplotypes, then check if they contribute to the current partial haplotype
						for (unsigned int s = 0; s < candidates.size(); s++) {

							bool matchPartialHap = true; //We assume that the candidate haplotype s matches with the partial haplotype, then check if it is true
							for (int k = 0; k < numLociCurrentPartialHap; k++) { //Loop over the loci and check if the alleles matches that of candidate haplotype s

								//Here we define the loci associated with loci position k, e.g. if there are 5 loci in total
								//and the partial hap covers just loci 2 and 3, then:
								//k = 0 corresponds to currentLoci = 2
								//k = 1 corresonds to currentLoci = 3
								int currentLoci = atoi((currentPartialHap[1+k]).c_str());

								//We now compare the allele at currentLoci in the candidate haplotype s with the allele at
								//position k in the partial haplotype. If they don't match, haplotype s cannot contribute
								//to haplotype j
								if ((candidates[s])[currentLoci] != (currentPartialHap[1+numLociCurrentPartialHap])[k]) {

									matchPartialHap = false;
									break;
								}
							}

							//If the candidate haplotype s matches at all loci in the partial hap, then we say that haplotype s contributes to the partial hap
							if (matchPartialHap == true) {
								contribsPartialHap.push_back(s);

								//Print out to check performance:
								//cout << "The full haplotype " << s << " with sequence " << candidates[s] << " was found as the contributor to the partial haplotype: ";
								//for(unsigned int k=0; k<currentPartialHap.size(); k++) {
								//	cout << currentPartialHap[k] << " ";
								//}
								//cout << "\n\n";
							}
						}
						contribsPartialHapSet.push_back(contribsPartialHap);
					}
					contribsPartitions.push_back(contribsPartialHapSet);
				}


				//What's the rationale for 50??
				for (unsigned int s = 0; s < 50*(positions_string.size()); s++)//set the 'initial conditions' for the most recent (unknown) candidate haplotype using pieces from the partial haplotype set until there are no undetermined loci (X) left
				{

					//To populate the most recent candidate haplotype we need to add partial haplotypes. First we find a random haplotype
					int randomPartialHapIndex = floor(multiLocusDataSplitLoci.size()*gsl_rng_uniform(rgen)); //Random number between 0 and n, where n is number of partial haplotypes
					vector<string> randomPartialHap = multiLocusDataSplitLoci[randomPartialHapIndex];

					//Find number of loci in the random haplotype
					int numLociPartialHap = atoi(randomPartialHap[0].c_str());

					//Find the sequence of the partial hap
					string seqPartialHap = randomPartialHap[1+numLociPartialHap];
					
	
					//Find the loci corresponding to sequence
					vector<double> lociPartialHap;
					for (int i=0; i<numLociPartialHap; i++) { //Loop over loci in haplotype

						//Append loci positions
						lociPartialHap.push_back(atoi(randomPartialHap[i+1].c_str()));

					}

					//Add this sequence to the most recent candidate sequence
					candidates[numCandidateHaplotypes - 1].replace(lociPartialHap[0], lociPartialHap.size(), seqPartialHap); //Replaces string[loci[0],loci[0]+loci.size()] with seq
				}
				
				double failure_counts = 0;//This defines the cut-off on when to stop making further attempts to optimize the haplotypes
				while (true) { //Infinite loop

					int randomHapIndex = floor(candidates.size()*gsl_rng_uniform(rgen)); //Random number between 0 and n, where n is number of candidate haplotypes
					int randomPartialHapIndex = floor(multiLocusDataSplitLoci.size()*gsl_rng_uniform(rgen)); //Random number between 0 and n, where n is number of partial haplotypes
					vector<string> randomPartialHap = multiLocusDataSplitLoci[randomPartialHapIndex];
					int numLociPartialHap = atoi((randomPartialHap[0]).c_str());
					
					//Obtain the sequence of the partial hap
					string seqPartialHap = randomPartialHap[1+numLociPartialHap];

					//Find the loci corresponding to sequence
					vector<double> lociPartialHap;
					for (int i=0; i<numLociPartialHap; i++) { //Loop over loci in partial haplotype

						lociPartialHap.push_back(atoi((randomPartialHap[i + 1]).c_str())); //Add loci i to pos
					}

					//Save the corresponding sequence in the random (full) haplotype so that we can revert the change later
					string randomHapOriginalString = candidates[randomHapIndex].substr(lociPartialHap[0], lociPartialHap.size());
					//Replace this sequence with that of the partial hap at these loci
					candidates[randomHapIndex].replace(lociPartialHap[0], lociPartialHap.size(), seqPartialHap);


					//Having changed the set of candidate haplotypes we need to recompute the contribution of the full haplotypes
					//to the partial haplotypes in order to be able to evaluate likelihoods.
					vector<vector<vector<int>>> contribsPartitions; 
					for (unsigned int i = 0; i < partitions.size(); i++) { //Loop over sets of partial haplotypes, denoted by i

						vector<vector<int>> contribsPartialHapSet; //Define contributions for this set
						for (unsigned int j = 0; j < (partitions[i]).size(); j++) { //Loop over partial haplotypes in set i

							vector<int> contribsPartialHap; //Define contributions for this haplotype (set i, partial hap j)
							vector<string> currentPartialHap = partitions[i][j];
							int numLociCurrentPartialHap = atoi(currentPartialHap[0].c_str());

							//Loop over the candidate haplotypes, then check if they contribute to the current partial haplotype
							for (unsigned int s = 0; s < candidates.size(); s++) {

								bool matchPartialHap = true; //We assume that the candidate haplotype s matches with the partial haplotype, then check if it is true
								for (int k = 0; k < numLociCurrentPartialHap; k++) { //Loop over the loci and check if the alleles matches that of candidate haplotype s
	
									//Here we define the loci associated with loci position k, e.g. if there are 5 loci in total
									//and the partial hap covers just loci 2 and 3, then:
									//k = 0 corresponds to currentLoci = 2
									//k = 1 corresonds to currentLoci = 3
									int currentLoci = atoi((currentPartialHap[1+k]).c_str());
	
									//We now compare the allele at currentLoci in the candidate haplotype s with the allele at
									//position k in the partial haplotype. If they don't match, haplotype s cannot contribute
									//to haplotype j
									if ((candidates[s])[currentLoci] != (currentPartialHap[1+numLociCurrentPartialHap])[k]) {
	
										matchPartialHap = false;
										break;
									}
								}
	
								//If the candidate haplotype s matches at all loci in the partial hap, then we say that haplotype s contributes to the partial hap
								if (matchPartialHap == true) {
									contribsPartialHap.push_back(s);
	
									//Print out to check performance:
									//cout << "The full haplotype " << s << " with sequence " << candidates[s] << " was found as the contributor to the partial haplotype: ";
									//for(unsigned int k=0; k<currentPartialHap.size(); k++) {
									//	cout << currentPartialHap[k] << " ";
									//}
									//cout << "\n\n";
								}
							}
							contribsPartialHapSet.push_back(contribsPartialHap);
						}
						contribsPartitions.push_back(contribsPartialHapSet);
					}


					//We here optimise the frequencies of the candidate haplotypes for the first (before transmission) timepoint.
					//Optimised frequencies stored in init_freqs1.
					//We use a noise parameter C. This generates a log likelihood.
					double L_new_before  = OptimumFreqs(0, candidates, init_freqs1, contribsPartitions, partitions, c, seed);

					//We here optimise the frequencies of the candidate haplotypes for the second (after transmission) timepoint.
					//Optimised frequencies stored in init_freqs2.
					//We use a noise parameter C. This generates a log likelihood.
					double L_new_after = OptimumFreqs(1, candidates, init_freqs2, contribsPartitions, partitions, c, seed);
	

					//Define total log likelihood at sum of log likelihood before and after transmission
					L_new = L_new_before + L_new_after;//the log-likelihood of the first and second timepoint are independent of each other

					//If the current likelihood is better than the best likelihood, then the change to the set of candidate haplotypes
					//was beneficial. We then store the new frequencies.
					if (L_new > L_best) {
						//cout << endl << L_new_before << " + " << L_new_after << " = " << L_new << "  " << L_best << endl;
						L_best = L_new;
						//cout << endl << "---->IMPROVE!" << endl;
						failure_counts = 1; //Why not zero?

						//Save the new frequencies
						init_freqs1_store = init_freqs1;
						init_freqs2_store = init_freqs2;

					} else if (L_new <= L_best) { //This change was worse than the best case so far, so we revert the changes and increment the failure_count
						
						//cout << endl << L_new_before << " + " << L_new_after << " = " << L_new << "  " << L_best << endl;
						candidates[randomHapIndex].replace(lociPartialHap[0], lociPartialHap.size(), randomHapOriginalString);
						//cout << endl << "---->FAIL!" << endl;
						failure_counts++; //After THRESHOLD failures, we exit optimisation
					}

					//Check if the number of failures has exceed the threshold, i.e. we need to stop optimisation.
					if (failure_counts > THRESHOLD) {
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
				cout << "likelihood found in attempt " << attempts << ": " << L_best << endl;
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
				
				
				//If we have more than one candidate haplotype, we need to check if they are unique, or, rather,
				//if the last added haplotype is unique from all the other ones.
				if (numCandidateHaplotypes > 1)
				{
					///Check for uniqueness
					bool unique = true; //Assume unique, then check
					for (unsigned int i = 0; i < numCandidateHaplotypes - 1; i++) { //Loop over the first n-1 haplotypes

						//Compare the nth haplotype with the ith haplotype (i in [0,n-1])	
						if (candidates[numCandidateHaplotypes-1] != candidates[i]) {
							
							unique = false;
							break;
						}
					}
					
					int bad_hap = 0;
					for (unsigned int i = 0; i < numCandidateHaplotypes; i++) //if the frequency of any haplotype is effectively zero before and after the transmission, then it is really not a good candidate and the code has got stuck in a local maximum
					{	
						if (init_freqs1_store[i] < 1e-9 && init_freqs2_store[i] < 1e-9)
						{
							bad_hap++;
						}
					}
					
					if (unique == true && bad_hap == 0)//if there is no duplicate haplotype and their frequencies BOTH before and after tranmission is above 10^-9, optimisation is normal
					{
						temp_likelihoods.push_back(L_best);
						freqs_1.push_back(init_freqs1_store);
						freqs_2.push_back(init_freqs2_store);
						temp_candidates.push_back(candidates);
					}
					else if (unique == false || bad_hap != 0)//otherwise, make another attempt to find a better set of haplotypes
					{
						attempts = attempts - 1;
					}
				}
				else if (numCandidateHaplotypes == 1) //this is the first haplotype so it is fine.
				{
					temp_likelihoods.push_back(L_best);
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
			cout << "***ROUND " << numCandidateHaplotypes << "***" << endl;
			cout << "number of haplotypes are currently = " << numCandidateHaplotypes << endl;
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
			DeltaBIC = -((-2*L_preserve_2)+numCandidateHaplotypes*log (reads_total)) + ((-2*L_preserve_1)+(numCandidateHaplotypes-1)*log (reads_total));
			if (DeltaBIC > 0) // if Delta_BIC > 0 --> the current set of N haplotypes are better than the N-1 haplotypes in the previous step
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
			
			if (DeltaBIC < 0)
			{
				cout << "The optimisation process stops here." << endl;
				cout << "The BIC in ROUND " << numCandidateHaplotypes << " is negative" << endl;
			}
		}
		
		return 0;

	} else {

		cout << "No Multi_locus_trajectories.out file found in " << multiLocusFilePath << ". Exiting bottleneck inference.\n";
		return false;
	}
}
