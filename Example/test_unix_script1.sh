t=0
mkdir Transmission1_qStar #making a directory to deposit all the reconstructed haplotypes (make sure to check your local directory and where you create the files)
cd Transmission1_qStar
for i in HA MP NA NP NS PA PB1 PB2; do #these are the names of the folders for the 8 segments of influenza A and B and are sorted alphabetically
    mkdir test_$t #reconstructed haplotypes for the ith segment (following the alphabetical order above) will be saved in a folder named test_i 
	cd test_$t
	echo $i #this shows the segment for which the haplotype reconstruction is currently taking place (note that the progress here could be slow depending on how large is the Multi_locus_trajectories.out
	/rds/user/mg878/hpc-work/package_lab/Codes/MLHapRec/./run_MLHapRec /rds/user/mg878/hpc-work/package_lab/Transmission1/$i/Multi_locus_trajectories.out 660
    cd ..
	t=`expr $t + 1`
done 