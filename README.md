# Estimating the bottleneck size using MLHapRec 
This repository contains a series of codes to infer the bottleneck size from short-read data using the haplotype reconstruction package `MLHapRec`.
Parts of codes are taken from the [Transmission_Project](https://bitbucket.org/casperlu/transmission_project) by C. K. Lumby and [SAMFIRE](https://github.com/cjri/samfire) by C. J. R. Illingworth.

These codes are developed for the purpose of studying transmission events in influenza virus with short-read data from two timepoints, i.e. one before the transmission in the donor population and one after in the recipient population.   

## Requirements
This code requires `Multi_locus_trajectories.out`, `Loci*.dat`, and `Hap_data*.dat` files with an inferred noise parameter `C` from `SAMFIER`. The files should be available for all the eight (flu type A and B) or seven (C and D) segments of the genome. Create a folder with the name of the segments, i.e. HA, MP, NA, NP, NS, PA, PB1, PB2, and put the corresponding files (mentioned above) into each folder. If some of the segments do not contain a variant, i.e. `Multi_locus_trajectories.out` is empty, please create the corresponding folder and leave it empty.

## Usage
To explain how to use this package and its features, we now consider the following example:
Suppose we have a dataset with one transmission event (we will later discuss how to deal with larger datasets). We store all the required files in a folder called `Transmission1`. 

**The first step is** to reconstruct the underlying haplotypes in each segment by typing the following in the command line:
```bash
t=0
mkdir Transmission1_qStar
cd Transmission1_qStar
for i in HA MP NA NP NS PA PB1 PB2; do
    mkdir test_$t
	cd test_$t
	echo $i
	/rds/user/mg878/hpc-work/package_lab/MLHapRec/./run_MLHapRec /rds/user/mg878/hpc-work/package_lab/Transmission1/$i/Multi_locus_trajectories.out 660
    cd ..
	t=`expr $t + 1`
done 
```
good?
