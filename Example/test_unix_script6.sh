cd /rds/user/mg878/hpc-work/package_lab/Transmission1_bottleneck
cat */likelihood_distribution.txt > catted_likelihoods.txt
/rds/user/mg878/hpc-work/package_lab/Codes/./run_avg_bottleneck catted_likelihoods.txt Transmission1_overal_likelihood.txt Transmission1_maximum_likelihood.txt