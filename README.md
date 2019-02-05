# Estimating the bottleneck size using MLHapRec 
This repository contains a series of codes to infer the bottleneck size from short-read data using the haplotype reconstruction package `MLHapRec`.
Some of the codes in this repository are taken from the [Transmission_Project](https://bitbucket.org/casperlu/transmission_project) by C. K. Lumby and [SAMFIRE](https://github.com/cjri/samfire) by C. J. R. Illingworth.

The codes in this repository are created to analyse transmission events in influenza viruses with short-read data from two timepoints, i.e. one before the transmission (in the donor population) and one after transmission (in the recipient population).   

<div style="text-align:center"><img src="overview.png" width="380">

## Requirements
This code requires `Multi_locus_trajectories.out`, `Loci*.dat`, and `Hap_data*.dat` files with an inferred noise parameter `C` from `SAMFIRE`. The files should be available for all the eight (flu type A and B) or seven (C and D) segments of the genome. Create a folder with the name of the segments, i.e. HA, MP, NA, NP, NS, PA, PB1, PB2, and put the corresponding files (mentioned above) into each folder. If some of the segments do not contain a variant, i.e. `Multi_locus_trajectories.out` is empty, please create the corresponding folder and leave it empty.

## Usage
To explain the features of this package, we consider the following example:
Suppose we have a dataset with one transmission event (we will later discuss how to deal with larger datasets). We store all the required files in a folder called `Transmission1`. 

* **The first step is** to construct the underlying haplotypes in each segment, q*, by typing the following in the command line:
```bash
t=0
mkdir Transmission1_qStar #making a directory to deposit all the reconstructed haplotypes (make sure to check your local directory and where you create the files)
cd Transmission1_qStar
for i in HA MP NA NP NS PA PB1 PB2; do #these are the names of the folders for the 8 segments of influenza A and B and are sorted alphabetically
    mkdir test_$t #reconstructed haplotypes for the ith segment (following the alphabetical order above) will be saved in a folder named test_i 
	cd test_$t
	echo $i #this shows the segment for which the haplotype reconstruction is currently taking place (note that the progress here could be slow depending on how large is the Multi_locus_trajectories.out
	/path/to/directory/Codes/MLHapRec/./run_MLHapRec /path/to/directory/Transmission1/$i/Multi_locus_trajectories.out 660
    cd ..
	t=`expr $t + 1`
done 
```
This creates a folder named `Transmission1_qStar` *(make sure the folder is in your local directory where other folders, including Transmission1, are located)* which contains 8 subfolders named `test_0`, `test_1`, ..., and `test_7` with the corresponding reconstructed haplotype sets saved in the `outcome_1.txt` for each segment `HA`, `MP`, `NA`, `NP`, `NS`, `PA`, `PB1`, and `PB2`, respectively.

Note that this step could take several minutes (hours) depending on how large your Multi_locus_trajectories.out files are.

* **The second step is** to create simulated short-read data, x*, by typing:
```bash
C=660 #change this to your inferred noise parameter 
for s2 in `seq 1 100`; do
## Update path to .run file
/path/to/directory/Codes/Codes_for_xStar/SimRealObsCount.run -if /path/to/directory/Transmission1 -of /path/to/directory/Transmission1_xStar -simOnly true -importHapsMahan /path/to/directory/Transmission1_qStar -format Mahan -rmMonoSim false -C $C -ss $s2 -filterData False 
done 
```
This creates a folder named `Transmission1_xStar` which contains 100 replica of the `Multi_locus_trajectories.out` for each segment of `Transmission1` assuming a noise parameter `C=660` and a fixed random seed generator `-ss $s2` for each replicate sample -- You can change the total number of replica by changing the for loop `for s2 in seq 1 NUMBER; do` and the noise parameter `C=VALUE`. A detailed explanation of what each flag does can be found [here](https://bitbucket.org/casperlu/transmission_project/).  

* **The third step is** to infer the frequency of the reconstructed haplotypes for each replicate sample x*. We call this the q** set. Note that in this step, there is no *haplotype reconstruction* for the simulated read, but we rather (re)infer the *frequency* of our initially reconstructed haplotype set `Transmission1_qStar` by typing the following commands:
```bash
mkdir Transmission1_qStarStar  #making a directory to deposit all the re-inferred haplotype frequencies q**
for s in `seq 1 100`; do #this varies from 1 to the total number of replicate samples which, in this case, we set to be equal to 100
   echo $s
   cd /path/to/directory/Transmission1_qStarStar
   mkdir Seed_$s #making a directory to deposit the re-inferred frequencies of set $s
   cd ..
   for t in `seq 0 7`; do
   Codes/./run_qstarstar /path/to/directory/Transmission1_xStar/Seed_$s/SimulatedData_Mahan_Gene_$t.dat /path/to/directory/Transmission1_qStar/test_$t/outcome_1.txt /path/to/directory/Transmission1_qStarStar/Seed_$s qStarStar_$t.txt 660
   done
done 
```
This creates a directory `Transmission1_qStarStar` with all the inferred haplotype frequencies q** (taken from x*) stored in folders `Seed_$s` where `$s` goes from 1 to 100 to cover all the 100 replicate samples.

We then concatenate the inferred frequencies of each segment across the 100 replicates and store them in a new folder called `Transmission1_bottleneck` with 8 subfolders corresponding to each gene segment, `segment_$t`, and concatenated file name `test_$t.txt` (where `$t` varies from 0 to 7) by typing the following in the command line:
```bash
mkdir Transmission1_bottleneck
for t in `seq 0 7`; do
   cd Transmission1_bottleneck
   mkdir segment_$t
   cd ..
   cat /path/to/directory/Transmission1_qStarStar/*/qStarStar_$t.txt > /path/to/directory/Transmission1_bottleneck/segment_$t/test_$t.txt
done 
```
* **The fourth step is** to calculate the bottleneck size for each segment by typing:
```bash
for t in `seq 0 7`; do
   cd /path/to/directory/Transmission1_bottleneck/segment_$t/
   /path/to/directory/Codes/./run_bottleneck /path/to/directory/Transmission1_qStar/test_$t/outcome_1.txt /path/to/directory/Transmission1_bottleneck/segment_$t/test_$t.txt likelihood_distribution.txt
   cd ..
done 
```
This calculates the liklihood function 
