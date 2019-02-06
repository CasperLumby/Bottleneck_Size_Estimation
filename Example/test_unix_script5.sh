for t in `seq 0 7`; do
   cd /rds/user/mg878/hpc-work/package_lab/Transmission1_bottleneck/segment_$t/
   /rds/user/mg878/hpc-work/package_lab/Codes/./run_bottleneck /rds/user/mg878/hpc-work/package_lab/Transmission1_qStar/test_$t/outcome_1.txt /rds/user/mg878/hpc-work/package_lab/Transmission1_bottleneck/segment_$t/test_$t.txt likelihood_distribution.txt
   cd ..
done 