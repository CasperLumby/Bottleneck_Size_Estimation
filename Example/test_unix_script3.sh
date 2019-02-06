mkdir Transmission1_qStarStar  #making a directory to deposit all the re-inferred haplotype frequencies q**
for s in `seq 1 100`; do #this varies from 1 to the total number of replicate samples which, in this case, we set to be equal to 100
   echo $s
   cd /rds/user/mg878/hpc-work/package_lab/Transmission1_qStarStar
   mkdir Seed_$s #making a directory to deposit the re-inferred frequencies of set $s
   cd ..
   for t in `seq 0 7`; do
   Codes/./run_qstarstar /rds/user/mg878/hpc-work/package_lab/Transmission1_xStar/Seed_$s/SimulatedData_Mahan_Gene_$t.dat /rds/user/mg878/hpc-work/package_lab/Transmission1_qStar/test_$t/outcome_1.txt /rds/user/mg878/hpc-work/package_lab/Transmission1_qStarStar/Seed_$s qStarStar_$t.txt 660
   done
done 